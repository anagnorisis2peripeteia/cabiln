"""Bounded, detached assembly results shared by edits, site checks and drawing."""

from collections import OrderedDict
from hashlib import sha256
import json
from threading import RLock

from rdkit import Chem

from pyPept.interfaces.reaction_library import REACTION_FINGERPRINT
from pyPept.peptide import Endpoint
from pyPept.site_chemistry import CHEMISTRY_FINGERPRINT


_MAX_BYTES = 32 * 1024 * 1024
_MAX_ENTRIES = 64
_PROPERTIES = int(Chem.PropertyPickleOptions.AllProps) & ~int(Chem.PropertyPickleOptions.ComputedProps)
_cache = OrderedDict()
_bytes = 0
_lock = RLock()


def assembly_key(peptide):
    digest = sha256(REACTION_FINGERPRINT + CHEMISTRY_FINGERPRINT.encode())
    definitions = {}
    for node in peptide.occurrences:
        if node.definition is None:
            raise ValueError('Assembly requires resolved monomer definitions')
        identity = id(node.definition)
        if identity not in definitions:
            template = node.definition.copy_template()
            definitions[identity] = sha256(
                template.ToBinary(_PROPERTIES)
                + repr(node.definition.attachment_metadata).encode()
            ).digest()
        digest.update(str(node.id).encode() + b':' + definitions[identity])
    digest.update(repr(tuple(edge.endpoints for edge in peptide.connections)).encode())
    return digest.digest()


def get_assembly(key):
    with _lock:
        value = _cache.get(key)
        if value is None:
            return None
        _cache.move_to_end(key)
    product, ports, labels = value
    return Chem.Mol(product), Chem.Mol(ports), {
        Endpoint(identity, slot): label for identity, slot, label in json.loads(labels)
    }


def put_assembly(key, product, ports, labels):
    global _bytes
    value = (
        product.ToBinary(_PROPERTIES), ports.ToBinary(_PROPERTIES),
        json.dumps([(site.occurrence_id, site.slot, label) for site, label in labels.items()]).encode(),
    )
    size = sum(map(len, value)) + len(key)
    if size > _MAX_BYTES:
        return
    with _lock:
        if key in _cache:
            _bytes -= sum(map(len, _cache.pop(key))) + len(key)
        _cache[key] = value
        _bytes += size
        while _bytes > _MAX_BYTES or len(_cache) > _MAX_ENTRIES:
            oldest, removed = _cache.popitem(last=False)
            _bytes -= len(oldest) + sum(map(len, removed))


def clear_assembly_cache():
    global _bytes
    with _lock:
        _cache.clear()
        _bytes = 0
