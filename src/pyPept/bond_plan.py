"""Batch the SMIRKS subset that only replaces two ports with one single bond.

All retained atoms and bonds must be identical in the reaction templates.
Context-sensitive queries, stereochemical edits and other transformations use
the general reaction executor. No independent list of reaction chemistry lives here.
"""

from functools import lru_cache

from rdkit import Chem

from pyPept.interfaces.reaction_library import _compiled_reaction


def _bonds(molecules):
    return {
        tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))): bond.GetBondType()
        for molecule in molecules for bond in molecule.GetBonds()
        if bond.GetBeginAtom().GetAtomMapNum() and bond.GetEndAtom().GetAtomMapNum()
    }


@lru_cache(maxsize=128)
def _join_queries(smirks):
    reaction = _compiled_reaction(smirks)
    if reaction.GetNumReactantTemplates() != 2 or reaction.GetNumProductTemplates() != 1:
        return None
    # Templates returned by RDKit borrow storage from the reaction cache.
    reactants = [Chem.Mol(reaction.GetReactantTemplate(index)) for index in range(2)]
    product = reaction.GetProductTemplate(0)
    atoms, ports = {}, []
    for molecule in reactants:
        dummies = [atom for atom in molecule.GetAtoms() if not atom.GetAtomicNum()]
        if len(dummies) != 1 or dummies[0].GetDegree() != 1 or dummies[0].GetAtomMapNum():
            return None
        ports.append((dummies[0].GetIdx(), dummies[0].GetNeighbors()[0].GetAtomMapNum()))
        for atom in molecule.GetAtoms():
            if not atom.GetAtomicNum():
                continue
            mapping = atom.GetAtomMapNum()
            symbol = atom.GetSymbol().lower() if atom.GetIsAromatic() else atom.GetSymbol()
            if not mapping or mapping in atoms or atom.GetSmarts() != f'[{symbol}:{mapping}]':
                return None
            atoms[mapping] = atom.GetSmarts()
    if product.GetNumAtoms() != len(atoms) or {
        atom.GetAtomMapNum(): atom.GetSmarts() for atom in product.GetAtoms()
    } != atoms:
        return None
    before, after = _bonds(reactants), _bonds([product])
    join = tuple(sorted(port[1] for port in ports))
    if join in before or after != {**before, join: Chem.BondType.SINGLE}:
        return None
    if any(bond.GetStereo() != Chem.BondStereo.STEREONONE
           for molecule in [*reactants, product] for bond in molecule.GetBonds()):
        return None
    return tuple((molecule, port[0]) for molecule, port in zip(reactants, ports))


def _matches(molecule, template, isotope):
    query, dummy = template
    targeted = Chem.RWMol(query)
    targeted.ReplaceAtom(dummy, Chem.AtomFromSmarts(f'[{isotope}#0]'))
    return molecule.HasSubstructMatch(targeted)


def batch_join(pool, connections, labels, reactions):
    """Return a new port graph, or None when any connection requires SMIRKS."""
    ports = {}
    for identity, molecule in pool.items():
        # Reaction products can discard disconnected spectator fragments.
        if len(Chem.GetMolFrags(molecule)) != 1:
            return None
        if any(molecule.GetAtomWithIdx(index).GetAtomicNum() == 0
               for bond in molecule.GetBonds() for index in bond.GetStereoAtoms()):
            return None
        for atom in molecule.GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetDegree() != 1:
                    return None
                anchor = atom.GetNeighbors()[0]
                if anchor.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                    return None
                ports[atom.GetIsotope()] = (identity, atom.GetIdx(), anchor.GetIdx())
    for edge, reaction in zip(connections, reactions):
        if len(reaction['steps']) != 1 or reaction.get('take_largest'):
            return None
        queries = _join_queries(reaction['steps'][0])
        if queries is None:
            return None
        if any(query.GetAtomWithIdx(dummy).GetIsotope() != reaction.get(key)
               for (query, dummy), key in zip(queries, ('slot_a', 'slot_b'))):
            return None
        left, right = edge.endpoints
        if not any(
            _matches(pool[left.occurrence_id], queries[a], labels[left])
            and _matches(pool[right.occurrence_id], queries[b], labels[right])
            for a, b in ((0, 1), (1, 0))
        ):
            return None
    product, offsets = Chem.RWMol(), {}
    for identity, molecule in pool.items():
        offsets[identity] = product.GetNumAtoms()
        product.InsertMol(molecule)
    consumed = []
    for index, edge in enumerate(connections):
        anchors = []
        for endpoint in edge.endpoints:
            identity, dummy, anchor = ports[labels[endpoint]]
            anchors.append(offsets[identity] + anchor)
            consumed.append(offsets[identity] + dummy)
        if anchors[0] == anchors[1] or product.GetBondBetweenAtoms(*anchors) is not None:
            return None
        product.AddBond(*anchors, Chem.BondType.SINGLE)
        product.GetBondBetweenAtoms(*anchors).SetIntProp('_connection_idx', index)
    product.BeginBatchEdit()
    for index in consumed:
        product.RemoveAtom(index)
    product.CommitBatchEdit()
    product.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(product)
    return product.GetMol()
