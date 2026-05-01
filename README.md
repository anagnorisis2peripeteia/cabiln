# pyPept — CABILN Fork

**Chemistry Aware BILN (CABILN) — a single-string notation for any peptide you can make in a lab**

> Development fork of [pyPept](https://github.com/Boehringer-Ingelheim/pyPept) by Boehringer Ingelheim.  
> Original publication: [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2), *J. Cheminformatics*, 2023.  
> Original authors: Rodrigo Ochoa, J.B. Brown, Thomas Fox.  
> Fork extensions: Cameron Beeley, 2026.

---

## Why CABILN?

BILN already handles crosslinks and cyclic peptides. CABILN adds three things BILN has no notation for: **sequential multi-step conjugation** (each reaction step refers to the last fragment, not the original residue), **sidechain branches** (a separate chain hanging off one residue), and **unambiguous crosslink syntax** (bond ID and R-group are visually separated so you can read the string without a decoder ring).

The result is that complex peptides that previously required separate files or hand-crafted CHUCKLES fragments can be written as a single string:

```
# GLP-1 agonist — 39-residue backbone, C20 isopeptide lipidation on Lys
# bracket reads synthesis order: K → gGlu(γ-COOH isopeptide) → AEEA → C20FA
Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am

# Maleimide conjugation then DBCO loading on Cys — two sequential reactions,
# each step's R-group refers to the preceding fragment, not to Cys directly
G-C[.Mal(4,1).DBCO(2,1)]-A-G-AzK-G

# Cyclic peptide with a disulfide bridge — head-to-tail (!1) + sidechain crosslink (!2)
# Second endpoint is just .!2 — R-groups are declared once on the first endpoint
!1-C.!2(4,4)-A-G-K-A-C.!2-!1

# Hydrocarbon-stapled helix (RCM, i,i+7)
ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am

# Full Fmoc-SPPS protection pattern — three residues, three different protecting groups
fmoc-C.trt(4,1)-K.boc(4,1)-R.pbf(4,1)-am
```

Every string above feeds directly into `Sequence(...)` and produces a valid RDKit molecule.

---

## The four CABILN extensions

### 1 — Inline sidechain modification `.Token(host_r, cap_r)`

Attach a cap, protecting group, or small molecule to **any R-group of any residue**, inline, without a separate bond table.

```
Residue.Modifier(host_Rgroup, modifier_Rgroup)
```

| Sequence | What it means |
|----------|---------------|
| `fmoc-C.trt(4,1)-am` | Cys, thiol (R4) capped with trityl R1 |
| `fmoc-K.boc(4,1)-am` | Lys, ε-amine (R4) capped with Boc R1 |
| `fmoc-R.pbf(4,1)-am` | Arg, guanidinium (R4) capped with Pbf R1 |
| `fmoc-C.acm(4,1)-am` | Cys, thiol (R4) with acetamidomethyl |
| `ac-pSer-am` | Pre-formed phosphoserine (no inline notation needed) |

Multiple residues capped in one string:

```python
seq = Sequence('fmoc-C.trt(4,1)-K.boc(4,1)-R.pbf(4,1)-am')
```

### 2 — Terminal cyclisation `!n-...-!n`

Head-to-tail and backbone cyclisation uses marker tokens at both endpoints. The bond ID `n` is just a label — any integer or short string.

```python
# Head-to-tail cyclic tetrapeptide
seq = Sequence('!1-A-G-K-A-!1')

# Lactam staple: Lys ε-amine (R4) → Asp sidechain carboxyl (R3)
seq = Sequence('ac-K.!1(4,3)-A-A-A-D.!1(3,4)-am')

# Disulfide bridge — first endpoint declares (r_self, r_other); second just .!1
seq = Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')

# Hydrocarbon staple (i, i+4): two olefinic residues, RCM closes the ring
seq = Sequence('ac-A-S5.!1(4,4)-A-A-A-S5.!1(4,4)-G-am')

# Hydrocarbon staple (i, i+7): mixed S5/R8 pair for longer helix
seq = Sequence('ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1(4,4)-G-am')
```

The `.!n(r1,r2)` suffix on the **first** endpoint names both R-groups of the bond; the second endpoint carries just `.!n` (the R-group is already determined). This avoids the ambiguity of old BILN `C(1,3)` notation where you couldn't tell which number was the bond ID and which was the R-group.

### 3 — Sequential bioconjugation bracket `Res[.A(r,s).B(t,u)…]`

The bracket notation chains modifications onto one residue in order. **Each step's R-group refers to the fragment immediately before it**, not to the original residue — so the string reads in the same order as the synthesis steps. Each bracket is independent, they can be mixed freely with crosslinks and branch notation, and multiple brackets can appear on the same peptide.

```
Host[.CapA(host_r, capA_r).CapB(capA_r, capB_r).CapC(capB_r, capC_r)]
```

**Two-step** — maleimide thiol-Michael then DBCO loading:

```python
# Step 1: Mal R1 occupies Cys R4 (thiol → thioether), Mal R2 is now exposed
# Step 2: DBCO R1 occupies Mal R2 — DBCO goes onto the maleimide, not onto Cys
seq = Sequence('G-C[.Mal(4,1).DBCO(2,1)]-A')
```

**Three-step** — fatty acid lipidation (semaglutide/retatrutide-type isopeptide linker):

```python
# Step 1: gGlu R4 (γ-COOH) bonds to Lys R4 (ε-amine) — isopeptide
# Step 2: AEEA R2 bonds to gGlu R1 (α-amine) — amide
# Step 3: C20FA R2 bonds to AEEA R1 (amine) — amide to fatty diacid
seq = Sequence('ac-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')
```

Reading the bracket left to right matches synthesis order from Lys outward: K → gGlu (isopeptide) → AEEA → C20FA.

**Brackets compose freely with everything else:**

```python
# Protected Cys with bracket, plus head-to-tail cyclisation
seq = Sequence('!1-C[.Mal(4,1)]-A-G-K-!1')

# Multiple independent brackets on one peptide (auto bond IDs: 100, 101, …)
seq = Sequence('fmoc-C[.trt(4,1)]-G-K[.boc(4,1)]-am')

# Bracket lipidation + disulfide staple on the same peptide
seq = Sequence('ac-C.!1(4,4)-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
```

### 4 — Sidechain branch `mainchain%%branch`

The `%%` separator attaches a pendant chain to the main chain via a crosslink. The branch string reads N→C, bonding residues through their backbone R-groups (R1→R2 at each `-`).

```python
# PEG branch: Lys ε-amine (R4) → PEG chain N-terminus (R1)
seq = Sequence('ac-K.!1(4,1)-G-am%!1-PEG-am')
```

Branch and bracket are interchangeable when the pendant chain attaches via a **backbone R-group** of the junction residue. For sidechain attachment (e.g. gGlu's γ-COOH is R4), bracket is more natural because it reads proximal→distal from K; branch works too but the string order is reversed (distal→proximal) because C20FA — a diacid cap with no R1 — must lead the branch:

```python
# Bracket: reads K → gGlu → AEEA → C20FA (synthesis order)
seq = Sequence('ac-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')

# Branch: same molecule, string reads C20FA → AEEA → gGlu (reversed)
seq = Sequence('ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1')
```

Full Retatrutide (GIP/GLP-1/glucagon triple agonist, 39 residues):

```python
seq = Sequence('Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am')
```

---

## Combining notations

The four extensions compose freely. Any crosslink, bracket, and branch can appear together on the same string:

```python
# Disulfide staple + 3-step lipidation on the same peptide
Sequence('ac-C.!1(4,4)-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
#              ^---disulfide crosslink---^  ^---isopeptide lipid linker---^

# Full Fmoc-SPPS protection: Pbf on Arg, Boc on Lys, disulfide between two Cys
Sequence('fmoc-R.pbf(4,1)-A-C.!1(4,4)-K.boc(4,1)-A-C.!1-am')
#               ^Arg guard  ^-disulfide-^  ^Lys guard  ^2nd Cys

# Hydrocarbon staple (RCM, i,i+4) plus isopeptide lipidation on the same helix
Sequence('ac-S5.!1(4,4)-A-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-S5.!1-am')
#              ^-------RCM staple-------^  ^--------lipid linker--------^

# Cyclic peptide with lipidation inside the ring
Sequence('!1-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1')
#          ^head-to-tail cycle^  ^lipidation bracket^
```

Bioconjugation handles (Mal, DBCO, TCO, Tz) follow the same pattern once the monomers are registered via `pyPept-monomer-add`.

---

## Chemistry types and bond validation

CABILN knows what chemistry each R-group represents. When you write a bond, the library checks whether the reaction is chemically sensible and warns you if it isn't.

### Protecting groups and caps

Standard Fmoc-SPPS protection pattern — these all validate silently:

```python
Sequence('fmoc-C.trt(4,1)-K.boc(4,1)-R.pbf(4,1)-am')
#         ^N-term  ^thiol    ^ε-amine   ^guanidinium ^C-term
```

### Crosslinks and macrolactams

```python
# Disulfide (thiol R4 ↔ thiol R4)
Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')

# Lactam staple (ε-amine R4 ↔ sidechain carboxyl R3)
Sequence('ac-K.!1(4,3)-A-A-A-D.!1(3,4)-am')

# Aspartimide (backbone ↔ sidechain carboxyl — five-membered ring)
Sequence('ac-A-D.!1(4,3)-G.!1-A-am')
# Glutarimide variant (six-membered ring)
Sequence('ac-A-E.!1(4,3)-G.!1-A-am')
```

### Bioconjugation handles

SPAAC (copper-free click — cyclooctyne + azide):

```python
# Azide handle on Lys sidechain, ready for strain-promoted cycloaddition
Sequence('G-AzK-G')

# Two-step: maleimide thiol-Michael then DBCO loading; azide-Lys in the same chain
# DBCO attaches to Mal R2, not to Cys — bracket notation makes the order explicit
Sequence('G-C[.Mal(4,1).DBCO(2,1)]-A-G-AzK-G')
```

IEDDA (tetrazine ligation — tetrazine + TCO):

```python
# Tetrazine on Lys, trans-cyclooctene on adjacent Lys
Sequence('ac-K.Tz(4,1)-A-A-K.TCO(4,1)-am')
```

Oxime ligation (aminooxy + aldehyde):

```python
# Pre-formed oxime linkage in a modified residue
Sequence('ac-A-Aoa-G-Ald-A-am')
```

Depsipeptide (ester backbone):

```python
# O-linked ester in place of one amide
Sequence('ac-A-Hser.!1(3,2)-G.!1-A-am')
```

Thioester:

```python
# S–C(=O) linkage — thiol R4 to carboxyl R2
Sequence('ac-C.!1(4,2)-A-G-!1-am')
```

### Phosphopeptides

Pre-formed phospho residues slot in like any other monomer — no special notation:

```python
Sequence('ac-pSer-am')          # phosphoserine
Sequence('ac-A-pSer-A-am')      # phospho-tripeptide
Sequence('ac-pTyr-am')          # phosphotyrosine
```

### Fatty acid lipidation (GLP-1 agonists)

Semaglutide/retatrutide-type lipidation routes through gGlu's γ-carboxyl (R4). Both bracket and branch produce the same molecule; bracket reads in synthesis order, branch string is reversed because C20FA (no R1) must lead:

```python
# Bracket — reads K → gGlu (isopeptide) → AEEA → C20FA
seq = Sequence('ac-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')

# Branch — same molecule, string reversed: C20FA → AEEA → gGlu ← K
seq = Sequence('ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1')
```

---

## How chemistry-aware validation works

Every bond in CABILN is validated before assembly. The system runs in two stages:

**Stage 1 — R-group labelling.** When a monomer is registered, SMARTS patterns are matched against the CHUCKLES fragment to label each R-group with a chemistry type. The full type set:

| Type | What it matches | Example monomer |
|------|-----------------|-----------------|
| `backbone_n` | α-amine (R1) | all amino acids |
| `backbone_c` | α-carboxyl (R2) | all amino acids |
| `backbone_n_mod` | N-methyl or secondary backbone N (R3) | Pro, Sar |
| `thiol` | –SH | Cys |
| `amine_primary` | –NH₂ sidechain | Lys, Orn |
| `carboxyl` | –COOH sidechain | Asp, Glu, gGlu |
| `guanidinium` | –C(=NH)NH₂ | Arg |
| `hydroxyl` | aliphatic –OH | Ser, Thr |
| `hydroxyl_phenolic` | phenol –OH | Tyr |
| `maleimide_c` | maleimide C= | MalAla, MalLys |
| `alkyne_c` | terminal alkyne –C≡CH | Pra, Hpg |
| `azide_alpha_c` | –CH₂N₃ | AzAla, AzHal, AzK |
| `cyclooctyne_c` | ring-strained alkyne | CyoAla, CyoLys |
| `tetrazine_c` | s-tetrazine | TzAla, TzLys |
| `tco_c` | trans-cyclooctene | TcoAla, TcoLys |
| `terminal_alkene` | –CH=CH₂ | S5, R8 (stapling olefins) |
| `aldehyde` | –CHO | aldehyde handles |
| `aminooxy` | –NHOH | oxime handles |
| `hydrazide` | –CONHNH₂ | hydrazone handles |
| `nhs_ester` | NHS-activated ester | activated linkers |
| `phosphate_p` | –PO(OH)₂ | pSer, pTyr |

**Stage 2 — bond validation.** When a bond is declared (crosslink, bracket, or inline cap), the pair of chemistry types is looked up in a reaction table. Reactions are classified as:

- **Silent** — standard peptide chemistry, no warning
- **Warned** — unusual but valid, `UserWarning` emitted
- **Rejected** — raises `ValueError`

### Supported reaction types

All 15 chemistry reactions in the YAML library, with 5 CABILN examples each:

#### backbone_amide — standard peptide amide bond

Reactant pairs: `backbone_n` + `backbone_c` · `amine_primary` + `backbone_c` · `amine_primary` + `carboxyl`

Forms at every `-` junction and as isopeptides when an amine sidechain bonds to a carboxyl.

<details>
<summary>5 examples</summary>

```python
# 1. Simple linear — backbone_n + backbone_c at every junction
Sequence('ac-A-G-K-L-V-am')

# 2. Head-to-tail cyclic — first backbone_n bonds to last backbone_c
Sequence('!1-A-G-K-A-G-!1')

# 3. Macrolactam — Lys ε-amine (R4) to Asp sidechain carboxyl (R3)
Sequence('ac-K.!1(4,3)-A-A-A-D.!1(3,4)-am')

# 4. Glutamate lactam — six-membered ring
Sequence('ac-K.!1(4,3)-A-A-A-A-E.!1(3,4)-am')

# 5. Lactam inside head-to-tail ring — bicyclic scaffold
Sequence('!1-K.!2(4,3)-A-G-A-D.!2(3,4)-!1')
```

</details>

#### backbone_ester — depsipeptide ester bond

Reactant pair: `backbone_o` + `backbone_c`

An α-hydroxy acid carries its R1 as `backbone_o`; the `-` junction with the next residue forms an ester instead of an amide. Register α-hydroxy acid monomers before use.

<details>
<summary>5 examples — Glc (glycolic acid) is in the library; register others from SMILES</summary>

```python
# Glc (glycolic acid, Gly-sized α-hydroxy acid) is already in the library.
# To add lactic acid (Ala-sized):
# register_monomer(smiles='OC(C)C(=O)O', symbol='Lac', name='Lactic acid')

# 1. Glycolic acid ester junction in a tetrapeptide
Sequence('ac-G-Glc-A-G-am')

# 2. Two glycolic acid junctions — alternating ester/amide backbone
Sequence('ac-A-Glc-A-Glc-A-am')

# 3. Glycolic acid flanked by charged residues
Sequence('ac-K-Glc-D-am')

# 4. Cyclic depsipeptide — head-to-tail ring with ester junction
Sequence('!1-Glc-A-G-A-!1')

# 5. Depsipeptide + disulfide staple
Sequence('ac-C.!1(4,4)-Glc-G-A-C.!1-am')
```

</details>

#### disulfide — disulfide bridge (S–S)

Reactant pair: `thiol` + `thiol`

Cys thiol R4 pairs with another Cys thiol R4. Selenol–selenol (selenocysteine) uses the same reaction.

<details>
<summary>5 examples</summary>

```python
# 1. Classic disulfide — five-residue loop
Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')

# 2. Short loop — adjacent cysteines
Sequence('ac-C.!1(4,4)-G-C.!1-am')

# 3. Two disulfide bridges — cystine ladder
Sequence('ac-C.!1(4,4)-C.!2(4,4)-A-G-A-C.!2-C.!1-am')

# 4. Disulfide inside head-to-tail ring — bicyclic
Sequence('!1-C.!2(4,4)-A-G-A-C.!2-!1')

# 5. Disulfide + isopeptide lipidation on same peptide
Sequence('ac-C.!1(4,4)-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
```

</details>

#### thioester — thioester crosslink (S–C=O)

Reactant pair: `thiol` + `backbone_c` · `thiol` + `carboxyl`

Sulphur displaces the amide nitrogen. Used in native chemical ligation models and cyclic depsipeptide analogues.

<details>
<summary>5 examples</summary>

```python
# 1. Cys thiol to backbone C-terminus — macrolactone mimic
Sequence('ac-C.!1(4,2)-A-G-!1-am')

# 2. Cys thiol to Asp sidechain carboxyl
Sequence('ac-C.!1(4,4)-A-A-D.!1-am')

# 3. Cys thiol to Glu sidechain — six-membered thioester ring
Sequence('ac-C.!1(4,4)-G-G-E.!1-am')

# 4. Longer macrocycle via thioester
Sequence('ac-C.!1(4,2)-A-A-A-G-!1-am')

# 5. Thioester inside head-to-tail cyclic peptide
Sequence('!1-C.!2(4,2)-A-G-!2-!1')
```

</details>

#### thioether_halide — thioether via SN2 on alkyl halide

Reactant pair: `thiol` + `alkyl_halide_c`

Cys thiol displaces a halide leaving group. Register a halide-bearing monomer with `alkyl_halide_c` at R4 before use.

<details>
<summary>5 examples — register halide monomer first</summary>

```python
# Register chloroacetyl-alanine
register_monomer(
    smiles='N[C@@H](CC(=O)Cl)C(=O)O',
    symbol='ClAcAla', name='Chloroacetyl-Ala',
)

# 1. Standard Cys alkylation by chloroacetyl handle
Sequence('ac-C.!1(4,4)-A-A-ClAcAla.!1-am')

# 2. Short loop — adjacent Cys and halide
Sequence('ac-C.!1(4,4)-G-ClAcAla.!1-am')

# 3. Thioether inside head-to-tail ring
Sequence('!1-C.!2(4,4)-A-G-ClAcAla.!2-!1')

# 4. Two thioether crosslinks — dialkylation bicycle
Sequence('ac-C.!1(4,4)-A-ClAcAla.!1-G-C.!2(4,4)-A-ClAcAla.!2-am')

# 5. Thioether + lipidation
Sequence('ac-C.!1(4,4)-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-ClAcAla.!1-am')
```

</details>

#### thiol_maleimide — thiol-maleimide Michael addition

Reactant pair: `thiol` + `maleimide_c`

Cys thiol adds across the maleimide double bond; the succinimide ring is retained. Standard ADC and PEGylation chemistry. `MalAla` and `MalLys` are in the library.

<details>
<summary>5 examples</summary>

```python
# 1. Cys + MalAla — succinimide crosslink
Sequence('ac-C.!1(4,4)-A-G-MalAla.!1-am')

# 2. Cys + MalLys — larger ring
Sequence('ac-C.!1(4,4)-A-A-A-MalLys.!1-am')

# 3. Bracket form — maleimide as first conjugation step
Sequence('ac-A-C[.MalAla(4,1)]-G-am')

# 4. Two thiol-maleimide crosslinks — parallel bridges
Sequence('ac-C.!1(4,4)-A-MalAla.!1-G-C.!2(4,4)-A-MalAla.!2-am')

# 5. Thiol-maleimide inside head-to-tail ring
Sequence('!1-C.!2(4,4)-A-G-MalAla.!2-A-!1')
```

</details>

#### nhs_ester_amide — NHS ester + amine → amide (NHS departs)

Reactant pair: `amine_primary` + `nhs_ester` · `amine_secondary` + `nhs_ester`

NHS-activated ester reacts with primary or secondary amine; NHS ring leaves as byproduct. Register NHS ester monomers before use.

<details>
<summary>5 examples — register NHS ester monomer first</summary>

```python
# Register NHS-glycine
register_monomer(
    smiles='NCC(=O)ON1C(=O)CCC1=O',
    symbol='GlyNHS', name='NHS-glycine',
)

# 1. Lys ε-amine + NHS-glycine — amide ligation
Sequence('ac-K.!1(4,4)-G-A-GlyNHS.!1-am')

# 2. Short loop — Lys adjacent to NHS ester
Sequence('ac-K.!1(4,4)-G-GlyNHS.!1-am')

# 3. Two NHS ligations — dual Lys modification
Sequence('ac-K.!1(4,4)-A-K.!2(4,4)-G-GlyNHS.!1-A-GlyNHS.!2-am')

# 4. NHS amide inside head-to-tail ring
Sequence('!1-K.!2(4,4)-A-G-GlyNHS.!2-!1')

# 5. NHS amide + lipidation on same peptide
Sequence('ac-K.!1(4,4)-A-GlyNHS.!1-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am')
```

</details>

#### oxime_ligation — aminooxy + aldehyde → oxime

Reactant pair: `aminooxy` + `aldehyde`

Aminooxy condenses with aldehyde to form a stable oxime (C=N–O); water is implicit. Bioorthogonal under mildly acidic conditions. Register `Aoa` and `Ald` monomers before use.

<details>
<summary>5 examples — register Aoa and Ald monomers first</summary>

```python
# Register aminooxy-alanine
register_monomer(smiles='N[C@@H](CON)C(=O)O', symbol='Aoa', name='Aminooxy-alanine')
# Register aldehyde-alanine
register_monomer(smiles='N[C@@H](CC=O)C(=O)O', symbol='Ald', name='Aldehyde-alanine')

# 1. Aoa + Ald — intramolecular oxime bridge
Sequence('ac-Aoa.!1(4,4)-A-G-A-Ald.!1-am')

# 2. Shorter oxime loop
Sequence('ac-Aoa.!1(4,4)-G-Ald.!1-am')

# 3. Oxime inside head-to-tail ring
Sequence('!1-Aoa.!2(4,4)-A-G-Ald.!2-A-!1')

# 4. Two oxime bonds — bis-oxime staple
Sequence('ac-Aoa.!1(4,4)-A-Ald.!1-G-Aoa.!2(4,4)-A-Ald.!2-am')

# 5. Oxime + disulfide on same peptide
Sequence('ac-C.!1(4,4)-A-Aoa.!2(4,4)-G-C.!1-Ald.!2-am')
```

</details>

#### hydrazone — hydrazide + aldehyde → hydrazone

Reactant pair: `hydrazide` + `aldehyde`

Hydrazide condenses with aldehyde forming a C=N–N hydrazone. Reversible at low pH — useful for stimuli-responsive linkers.

<details>
<summary>5 examples — register hydrazide monomer first</summary>

```python
# Register hydrazide-glycine
register_monomer(smiles='NCC(=O)NN', symbol='GlyHz', name='Hydrazide-Gly')
# Ald monomer: same SMILES as oxime examples above

# 1. GlyHz + Ald — hydrazone bridge
Sequence('ac-GlyHz.!1(4,4)-A-G-A-Ald.!1-am')

# 2. Shorter loop
Sequence('ac-GlyHz.!1(4,4)-G-Ald.!1-am')

# 3. Hydrazone inside head-to-tail ring
Sequence('!1-GlyHz.!2(4,4)-A-G-Ald.!2-!1')

# 4. Two hydrazone bonds — bifunctional staple
Sequence('ac-GlyHz.!1(4,4)-A-Ald.!1-G-GlyHz.!2(4,4)-A-Ald.!2-am')

# 5. Hydrazone + disulfide combination
Sequence('ac-C.!1(4,4)-A-GlyHz.!2(4,4)-G-C.!1-Ald.!2-am')
```

</details>

#### cuaac_1_4_triazole — CuAAC (terminal alkyne + azide → 1,4-triazole)

Reactant pair: `alkyne_c` + `azide_alpha_c`

Copper-catalysed azide–alkyne cycloaddition. `Pra` (propargylglycine) and `Hpg` (homopropargylglycine) carry the alkyne; `AzAla`, `AzHal`, and `AzK` carry the azide.

<details>
<summary>5 examples</summary>

```python
# 1. Pra × AzK — triazole staple via Lys sidechain azide
Sequence('ac-Pra.!1(4,4)-A-A-A-AzK.!1-am')

# 2. Hpg × AzAla — homologated alkyne, shorter ring
Sequence('ac-Hpg.!1(4,4)-G-AzAla.!1-am')

# 3. Pra × AzAla — flanked by Ala
Sequence('ac-A-Pra.!1(4,4)-G-AzAla.!1-A-am')

# 4. Two CuAAC crosslinks — dual triazole bicycle
Sequence('ac-Pra.!1(4,4)-A-Pra.!2(4,4)-A-AzK.!1-G-AzAla.!2-am')

# 5. CuAAC + isopeptide lipidation
Sequence('ac-Pra.!1(4,4)-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-AzK.!1-am')
```

</details>

#### spaac_triazole — SPAAC (strained cyclooctyne + azide, copper-free)

Reactant pair: `cyclooctyne_c` + `azide_alpha_c`

Strain-promoted [3+2] cycloaddition — no copper required, bioorthogonal in live cells. `CyoAla` and `CyoLys` carry the strained alkyne.

<details>
<summary>5 examples</summary>

```python
# 1. CyoAla × AzAla — copper-free click staple
Sequence('ac-CyoAla.!1(4,4)-G-G-AzAla.!1-am')

# 2. CyoLys × AzK — both Lys sidechains, longer ring
Sequence('ac-CyoLys.!1(4,4)-A-A-AzK.!1-am')

# 3. CyoAla × AzK — mixed pair
Sequence('ac-A-CyoAla.!1(4,4)-G-A-AzK.!1-A-am')

# 4. Two SPAAC crosslinks — parallel bicyclisation
Sequence('ac-CyoAla.!1(4,4)-A-CyoAla.!2(4,4)-A-AzAla.!1-G-AzAla.!2-am')

# 5. SPAAC + disulfide on same peptide
Sequence('ac-C.!1(4,4)-A-CyoAla.!2(4,4)-G-C.!1-AzAla.!2-am')
```

</details>

#### iedda_tetrazine_tco — IEDDA (tetrazine + TCO → dihydropyridazine)

Reactant pair: `tetrazine_c` + `tco_c`

Inverse-electron-demand Diels–Alder; fastest bioorthogonal reaction. Two-step: [4+2] then retro-[4+2] expelling N₂. `TzAla`/`TzLys` carry the tetrazine; `TcoAla`/`TcoLys` carry the TCO.

<details>
<summary>5 examples</summary>

```python
# 1. TzAla × TcoAla — IEDDA staple
Sequence('ac-TzAla.!1(4,4)-A-A-TcoAla.!1-am')

# 2. TzLys × TcoLys — both Lys sidechains, longer ring
Sequence('ac-TzLys.!1(4,4)-A-A-A-TcoLys.!1-am')

# 3. TzAla × TcoLys — mixed pair
Sequence('ac-A-TzAla.!1(4,4)-G-A-TcoLys.!1-A-am')

# 4. Two IEDDA crosslinks — bis-tetrazine scaffold
Sequence('ac-TzAla.!1(4,4)-A-TzLys.!2(4,4)-A-TcoAla.!1-G-TcoLys.!2-am')

# 5. IEDDA + disulfide — orthogonal bioorthogonal handles
Sequence('ac-C.!1(4,4)-A-TzAla.!2(4,4)-G-C.!1-TcoAla.!2-am')
```

</details>

#### phosphorylation — O-phosphorylation

Reactant pair: `hydroxyl` + `phosphate_p`

Phosphorylation of Ser, Thr, or Tyr. Pre-formed `pSer`, `pThr`, `pTyr` monomers cover the common case — no bond notation required.

<details>
<summary>5 examples</summary>

```python
# 1. Phosphoserine — pre-formed monomer, drops straight in
Sequence('ac-A-pSer-A-am')

# 2. Phosphothreonine
Sequence('ac-A-pThr-A-am')

# 3. Phosphotyrosine
Sequence('ac-A-pTyr-A-am')

# 4. Doubly phosphorylated peptide — pSer + pThr
Sequence('ac-pSer-A-G-pThr-A-am')

# 5. Phosphopeptide with disulfide staple
Sequence('ac-C.!1(4,4)-A-pSer-G-pThr-C.!1-am')
```

</details>

#### rcm_alkene — ring-closing olefin metathesis (RCM, all-hydrocarbon staple)

Reactant pair: `terminal_alkene` + `terminal_alkene`

Ruthenium-catalysed RCM closes two terminal alkenes; ethylene byproduct is discarded. `S5` (i,i+4) and `R8` (i,i+7) are in the library.

<details>
<summary>5 examples</summary>

```python
# 1. S5–S5 pair, i,i+4 — standard one-turn helix staple
Sequence('ac-A-S5.!1(4,4)-A-A-A-S5.!1-G-am')

# 2. S5–R8 pair, i,i+7 — two-turn helix staple
Sequence('ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1-G-am')

# 3. S5–S5 at N-terminus
Sequence('ac-S5.!1(4,4)-A-A-A-S5.!1-am')

# 4. Hydrocarbon staple + isopeptide lipidation — stapled GLP-1 analogue motif
Sequence('ac-S5.!1(4,4)-A-A-K[.gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-S5.!1-am')

# 5. Staple + disulfide — orthogonal RCM and thiol crosslinks
Sequence('ac-C.!1(4,4)-A-S5.!2(4,4)-A-A-A-S5.!2-A-C.!1-am')
```

</details>

#### imide_n_acylation — cyclic imide (backbone N acylated by sidechain carboxyl)

Reactant pair: `backbone_n_mod` (R3) + `carboxyl` (R4)

Asp or Glu sidechain carboxyl acylates the backbone amide NH of the adjacent residue, forming a succinimide (Asp → 5-membered) or glutarimide (Glu → 6-membered) ring embedded in the backbone.

<details>
<summary>5 examples</summary>

```python
# 1. Aspartimide — Asp sidechain carboxyl acylates adjacent backbone N
Sequence('ac-A-D.!1(4,3)-G.!1-A-am')

# 2. Glutarimide — Glu gives a 6-membered imide ring
Sequence('ac-A-E.!1(4,3)-G.!1-A-am')

# 3. Aspartimide with Ala as the acylated residue
Sequence('ac-D.!1(4,3)-A.!1-G-am')

# 4. Two aspartimide rings in tandem — bis-imide
Sequence('ac-D.!1(4,3)-G.!1-A-D.!2(4,3)-G.!2-am')

# 5. Aspartimide inside head-to-tail cyclic peptide
Sequence('!1-A-D.!2(4,3)-G.!2-A-!1')
```

</details>

---

## How the R-group system works

Every monomer in the library has numbered R-groups (attachment points). The standard mapping for α-amino acids is:

| R-group | Role |
|---------|------|
| R1 | Backbone N (incoming amide bond) |
| R2 | Backbone C (outgoing amide bond / carboxyl) |
| R3 | N-modification / backbone N cap |
| R4+ | Sidechain (thiol, amine, carboxyl, …) |

When you write `C.trt(4,1)`, you are saying: form a bond between **Cys R4** (the thiol) and **trt R1** (the attachment point of the trityl group). The library looks up both monomers, finds their chemical types, and validates that thiol–C is a sensible bond.

The same `(host_r, cap_r)` pair controls all bond formation — caps, crosslinks, and bracket chains all use identical syntax.

---

## Getting started

### Installation

```bash
pip install git+https://github.com/anagnorisis2peripeteia/pyPept.git
```

Development install:

```bash
git clone https://github.com/anagnorisis2peripeteia/pyPept.git
cd pyPept
pip install -e ".[dev]"
```

Requires Python ≥ 3.9, RDKit, and BioPython.

### Build a molecule

```python
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit import Chem

seq = Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')
mol = Molecule(seq)
rdmol = mol.get_molecule(fmt='ROMol')
print(Chem.MolToSmiles(rdmol))
```

### Validate before building

```python
report = Sequence.validate('ac-K.!1(4,3)-A-A-A-D.!1(3,4)-am')
print(report.ok)        # True
print(report.bonds)     # [(K, R4, D, R3)]
print(report.warnings)  # any chemistry alerts
```

---

## Monomer pipeline — add any structure

If a residue isn't in the library, register it from SMILES in one call. The pipeline auto-detects backbone N, backbone C, and all sidechain attachment points using a SMARTS graph-topology rule — no manual CHUCKLES authoring.

### Pre-activate (inspect R-groups first)

```python
from pyPept.interfaces.monomer_pipeline import pre_activate

result = pre_activate('N[C@@H](CS)C(=O)O')
print(result.chuckles)    # '[1*]N([3*])[C@@H](CS[4*])C([2*])=O'
print(result.leaving)     # {1: '[H]', 2: '[OH]', 3: '[H]', 4: '[H]'}
print(result.chem_types)  # {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod', 4: 'thiol'}
```

R1/R2 are assigned by graph topology (the unique path between the amino N and the carboxyl C), so the rule handles α-, β-, and γ-amino acids, N-methyl, and Aib-type residues without hard-coded stereo assumptions.

### Register to library

```python
from pyPept.interfaces.cli_monomer import register_monomer

register_monomer(
    smiles='N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O',
    symbol='Lys_Boc',
    name='Boc-Lysine',
)
# Appends to monomers.sdf — existing records untouched
```

### CLI registration

```bash
# From SMILES
pyPept-monomer-add --smiles "N[C@@H](CS)C(=O)O" --symbol Cys_check

# From an existing CABILN fragment (assembles first, then registers the product)
pyPept-monomer-add --from-cabiln "K.boc(4,1)" --symbol Lys_Boc --name "Boc-lysine"
```

---

## Legacy BILN crosslink notation

Old BILN used bare integers: `C(1,3)-A-A-A-C(1,3)`. The `1` was the bond ID and the `3` was the R-group with no delimiter — scanning left to right you couldn't tell which was which.

CABILN separates them with `.!`:

| Old BILN | CABILN equivalent |
|----------|-------------------|
| `C(1,3)-A-C(1,3)` | `C.!1(3,3)-A-C.!1` |
| `K(2,4)-A-D(2,3)` | `K.!2(4,3)-A-D.!2` |

Old BILN is accepted, but you must declare it explicitly — `Sequence()` without a format flag raises `ValueError` to avoid silent misparse:

```python
from pyPept.sequence import Sequence, biln_to_cabiln

# Raises ValueError — old notation must be declared
Sequence('C(1,3)-A-A-A-C(1,3)')

# Pass fmt='biln' to auto-convert and continue
seq = Sequence('C(1,3)-A-A-A-C(1,3)', fmt='biln')

# Or convert upfront and use CABILN directly
cabiln = biln_to_cabiln('C(1,3)-A-A-A-C(1,3)')
seq = Sequence(cabiln)   # C.!1(3,3)-A-A-A-C.!1
```

---

## Compatibility

| Input | Status | Notes |
|-------|--------|-------|
| Linear BILN `A-G-K` | Works | Unchanged from upstream |
| Capped BILN `fmoc-A-G-K-am` | Works | Unchanged |
| CABILN crosslinks `.!n(r,s)` | Works | New in this fork |
| CABILN brackets `[.A(r,s).B(t,u)]` | Works | New in this fork |
| CABILN branches `%%` | Works | New in this fork |
| Old BILN crosslinks `C(1,3)` | Requires `fmt='biln'` | Raises `ValueError` otherwise; auto-converts when declared |
| FASTA `ACDEFGHIKLMNPQRSTVWY` | Works | Via `run_pyPept --fasta` |
| HELM linear `PEPTIDE1{...}$$$$V2.0` | Works | Converter produces CABILN |
| HELM crosslinks `$PEPTIDE1,PEPTIDE1,…` | Works | Converter emits `.!n(r,s)` CABILN |

---

## CLI tools

| Command | Description |
|---------|-------------|
| `run_pyPept --biln <seq>` | 2D/3D structure from BILN or CABILN |
| `run_pyPept --fasta <seq>` | Structure from FASTA single-letter codes |
| `run_pyPept --helm <seq>` | Structure from HELM (linear or crosslinked) |
| `pyPept-BILN-validate --biln <seq>` | Validate a CABILN string, list bonds |
| `pyPept-monomer-add --smiles <smi> --symbol <tok>` | Register monomer from SMILES |
| `pyPept-monomer-add --from-cabiln <seq> --symbol <tok>` | Register from CABILN fragment |

---

## Tests

```bash
pytest tests/test_bond_validation_and_assembly.py -v
```

271 tests covering: bond validation, inline caps, bracket notation (1/2/3-step), crosslinks, RCM staples, SPAAC/IEDDA/oxime/hydrazone/depsipeptide/thioester chemistry, fatty acid branching, phosphopeptides, monomer pipeline, CLI, and HELM/BILN/FASTA round-trips.

---

## References

- [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2), *Journal of Cheminformatics*, 2023.
- [BILN — A Human-readable Line Notation for Complex Peptides](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00703), *J. Chem. Inf. Model.*, 2022.
