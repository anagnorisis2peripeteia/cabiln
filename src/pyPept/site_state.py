"""Current attachment chemistry on an assembled graph with unconsumed ports.

The graph retains placeholders for each free leaving group. Numbering and
occupancy come from the resolved peptide; reaction context comes from assembly.
Several ports can share one anchor and its cached chemical assessment.
"""

from pyPept.attachments import _SUPPORTED_TYPES, reaction_for_types
from pyPept.peptide import Endpoint
from pyPept.site_chemistry import Perception


_CONSUMABLE_GROUPS = frozenset({
    'alkyne_c', 'cyclooctyne_c', 'azide_alpha_c', 'tetrazine_c', 'tco_c',
    'nhs_ester', 'maleimide_c', 'terminal_alkene', 'thia_michael_c',
    'aldehyde', 'formamide_c',
})
_ROLES = frozenset({'backbone_n', 'backbone_n_mod', 'backbone_c', 'backbone_o'})


class CurrentSites:
    def __init__(self, peptide, molecule, labels):
        self.peptide = peptide
        self.molecule = molecule
        self.labels = labels
        self.perception = Perception(molecule, targeted=True)
        self.origins = {
            (atom.GetIntProp('_residue_idx'), atom.GetIntProp('_template_atom')): atom.GetIdx()
            for atom in molecule.GetAtoms()
            if atom.HasProp('_residue_idx') and atom.HasProp('_template_atom')
        }
        self.dummies = {atom.GetIsotope() for atom in molecule.GetAtoms() if not atom.GetAtomicNum()}
        self.states = {}

    def for_occurrence(self, identity):
        node = self.peptide.occurrence(identity)
        return [self.site(Endpoint(identity, site.slot)) for site in node.sites]

    def site(self, endpoint):
        if endpoint not in self.states:
            node = self.peptide.occurrence(endpoint.occurrence_id)
            site = self.peptide.site(endpoint)
            metadata = next(item for item in node.definition.attachment_metadata if item[0] == endpoint.slot)
            _, original, leaving, declared = metadata
            index = self.origins.get((node.id, site.anchor))
            label = self.labels[endpoint]
            used = self.peptide.connection_at(endpoint) is not None
            role = original if original in _ROLES else 'sidechain'
            kind = 'unknown'
            if index is not None:
                kind = self.perception.classify(index, label, leaving, original).functionality
            reaction_type = kind
            if role in ('backbone_c', 'backbone_o') or (
                role in ('backbone_n', 'backbone_n_mod')
                and kind in ('amine_primary', 'amine_secondary')
            ):
                reaction_type = role
            consumed = index is None or (not used and label not in self.dummies) or (
                original in _CONSUMABLE_GROUPS and kind != original
            )
            supported = not consumed and reaction_type in _SUPPORTED_TYPES
            reason = 'This site is already connected' if used else (
                'This reactive group was consumed by another connection' if consumed else
                f'No supported connection for {kind.replace("_", " ")}' if not supported else ''
            )
            self.states[endpoint] = {
                'slot': endpoint.slot, 'chem_type': reaction_type,
                'functionality': kind, 'role': role, 'used': used,
                'supported': supported, 'reason': reason, 'leaving': leaving or '',
                'declared_chem_type': declared,
            }
        return dict(self.states[endpoint])


def connection_reaction(left, right):
    """Assess actual sites; the caller must still validate the resulting product."""
    for site in (left, right):
        if site['used'] or not site['supported']:
            raise ValueError(f"R{site['slot']}: {site['reason']}")
    reaction = reaction_for_types(left['chem_type'], right['chem_type'])
    if reaction is None:
        raise ValueError(f"No reaction for {left['chem_type']} + {right['chem_type']}")
    return reaction
