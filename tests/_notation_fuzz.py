"""Small constructive notation inputs; expected edges precede serialization."""

from dataclasses import dataclass, replace
from functools import lru_cache

from hypothesis import strategies as st

from pyPept.peptide import Connection, Endpoint, Peptide
from pyPept.sequence import Sequence


TOKENS = (
    "G", "Gly", "A", "Ala", "K", "Lys", "C", "Cys", "Ser", "Val",
    "<[1*]N[C@@H]([13CH2][O-])C([2*])=O>",
    "<[1*]N[C@H](C[NH3+])C([2*])=O>",
    "<[1*]N[C@@H]([2H])C([2*])=O>",
)


def edge_set(edges):
    return {tuple(sorted(((a, sa), (b, sb)))) for a, sa, b, sb in edges}


def observed_edges(peptide):
    return {
        tuple((end.occurrence_id, end.slot) for end in edge.endpoints)
        for edge in peptide.connections
    }


@lru_cache(maxsize=64)
def singleton(token):
    sequence = Sequence(token, track_source=True, warning_sink=lambda _: None)
    return Peptide.from_sequence(sequence).occurrences[0]


@dataclass(frozen=True)
class GraphCase:
    family: str
    tokens: tuple
    edges: tuple

    def peptide(self, order):
        nodes = tuple(
            replace(singleton(self.tokens[i]), id=i, source=None, order_key=())
            for i in order
        )
        edges = tuple(
            Connection(Endpoint(a, sa), Endpoint(b, sb))
            for a, sa, b, sb in self.edges
        )
        return Peptide(nodes, edges)

    def tagged(self, order, prefix, bare):
        annotations = {i: [] for i in range(len(self.tokens))}
        positions = {identity: index for index, identity in enumerate(order)}
        for number, (a, sa, b, sb) in enumerate(self.edges):
            if positions[a] > positions[b]:
                a, sa, b, sb = b, sb, a, sa
            label = f"!{prefix}{number}"
            annotations[a].append(f".{label}({sa},{sb})")
            annotations[b].append(
                f".{label}" if bare else f".{label}({sb},{sa})"
            )
        return "%".join(self.tokens[i] + "".join(annotations[i]) for i in order)


@st.composite
def graph_cases(draw):
    family = draw(st.sampled_from((
        "chain", "cycle", "crosslink", "parallel", "duplicates", "tree",
        "unusual_slot",
    )))
    count = draw(st.integers(2, 6))
    tokens = draw(st.lists(st.sampled_from(TOKENS), min_size=count, max_size=count))
    edges = [(i, 2, i + 1, 1) for i in range(count - 1)]
    if family == "cycle":
        edges.append((count - 1, 2, 0, 1))
    elif family in ("crosslink", "parallel"):
        partner = 1 if family == "parallel" else count - 1
        tokens[0] = tokens[partner] = draw(st.sampled_from(("C", "Cys")))
        edges.append((0, 4, partner, 4))
    elif family == "duplicates":
        tokens *= 2
        edges += [(a + count, sa, b + count, sb) for a, sa, b, sb in edges]
    elif family == "tree":
        tokens = [draw(st.sampled_from(("K", "Lys"))) for _ in range(count)]
        edges, available = [], [(0, 1), (0, 4)]
        for child in range(1, count):
            parent, slot = draw(st.sampled_from(available))
            available.remove((parent, slot))
            edges.append((parent, slot, child, 2))
            available.extend(((child, 1), (child, 4)))
    elif family == "unusual_slot":
        slot = draw(st.sampled_from((7, 17, 31, 257)))
        tokens[0] = f"<[1*]N([{slot}*])[C@@H]([13CH2][O-])C([2*])=O>"
        tokens.append("ac")
        edges.append((0, slot, count, 2))
    return GraphCase(family, tuple(tokens), tuple(edges))


@dataclass(frozen=True)
class BranchCase:
    tokens: tuple
    children: tuple
    protected: tuple

    @property
    def graph(self):
        edges = [(0, 2, 1, 1)]
        edges += [
            (host, slot, child, 2)
            for host, children in enumerate(self.children)
            for slot, child in children
        ]
        return GraphCase("recursive", self.tokens, tuple(edges))

    def spelling(self, mode):
        order, protected_ids = [0, 1], set()

        def body(identity, slot, inherited):
            order.append(identity)
            inherited |= self.protected[identity]
            if inherited:
                protected_ids.add(identity)
            value = f"{self.tokens[identity]}({slot},2)"
            children = self.children[identity]
            if mode == "hub" and children:
                value = "[" + value + "]"
            for child_slot, child in children:
                if (mode == "implied" and len(children) == 1
                        and not self.protected[child]):
                    value += "." + body(child, child_slot, inherited)
                else:
                    value += group(child, child_slot, inherited)
            return value

        def group(identity, slot, inherited, root=False):
            value = body(identity, slot, inherited)
            if self.protected[identity]:
                return ".{" + value + "}"
            if mode == "legacy" and not root:
                return "[." + value + "]"
            return ".[" + value + "]"

        source = self.tokens[0]
        for slot, child in self.children[0]:
            source += group(child, slot, False, root=True)
        return source + "-" + self.tokens[1], tuple(order), protected_ids


@st.composite
def branch_cases(draw):
    count = draw(st.integers(3, 9))
    tokens = [draw(st.sampled_from(("K", "Lys"))) for _ in range(count)]
    tokens[1] = draw(st.sampled_from(TOKENS))
    children, available = [[] for _ in tokens], [(0, 4)]
    for child in range(2, count):
        host, slot = draw(st.sampled_from(available))
        available.remove((host, slot))
        children[host].append((slot, child))
        available.extend(((child, 1), (child, 4)))
    protected = draw(st.lists(st.booleans(), min_size=count, max_size=count))
    return BranchCase(tuple(tokens), tuple(map(tuple, children)), tuple(protected))
