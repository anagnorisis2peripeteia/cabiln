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
# Retatrutide — 39-residue backbone, K16-K17 with C20 lipid conjugate on K17
# bracket reads synthesis order: K17 → AEEA → gGlu(γ-COOH isopeptide) → C20FA
Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K-K.[AEEA(4,2).gGlu(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am

# Maleimide conjugation then AEEA linker on Cys — two sequential reactions,
# each step's R-group refers to the preceding fragment, not to Cys directly
G-C.[Mal(4,4).AEEA(2,1)]-A-G-AzK-G

# Cyclic peptide with a disulfide bridge — head-to-tail (!1) + sidechain crosslink (!2)
# Second endpoint is just .!2 — R-groups are declared once on the first endpoint
!1-C.!2(4,4)-A-G-K-A-C.!2-!1

# Hydrocarbon-stapled helix (RCM, i,i+7)
ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1-G-am

# Full Fmoc-SPPS protection pattern — three residues, three different protecting groups
fmoc-C.trt(4,2)-K.boc(4,2)-R.pbf(4,2)-am
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
| `fmoc-C.trt(4,2)-am` | Cys, thiol (R4) capped with trityl R2 |
| `fmoc-K.boc(4,2)-am` | Lys, ε-amine (R4) capped with Boc R2 |
| `fmoc-R.pbf(4,2)-am` | Arg, guanidinium (R4) capped with Pbf R2 |
| `fmoc-C.acm(4,2)-am` | Cys, thiol (R4) with acetamidomethyl |
| `ac-pSer-am` | Pre-formed phosphoserine (no inline notation needed) |

Multiple residues capped in one string:

```python
seq = Sequence('fmoc-C.trt(4,2)-K.boc(4,2)-R.pbf(4,2)-am')
```

### 2 — Terminal cyclisation `!n-...-!n`

Head-to-tail and backbone cyclisation uses marker tokens at both endpoints. The bond ID `n` is just a label — any integer or short string.

```python
# Head-to-tail cyclic tetrapeptide
seq = Sequence('!1-A-G-K-A-!1')

# Lactam staple: Lys ε-amine (R4) → Asp sidechain carboxyl (R4)
seq = Sequence('ac-K.!1(4,4)-A-A-A-D.!1-am')

# Disulfide bridge
seq = Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')

# Hydrocarbon staple (i, i+4): two olefinic residues, RCM closes the ring
seq = Sequence('ac-A-S5.!1(4,4)-A-A-A-S5.!1-G-am')

# Hydrocarbon staple (i, i+7): mixed S5/R8 pair for longer helix
seq = Sequence('ac-A-S5.!1(4,4)-A-A-A-A-A-R8.!1-G-am')
```

The `.!n(r1,r2)` suffix on the **first** endpoint names both R-groups of the bond; the second endpoint carries just `.!n` (the R-group is already determined). This avoids the ambiguity of old BILN `C(1,3)` notation where you couldn't tell which number was the bond ID and which was the R-group.

### 3 — Sequential bioconjugation bracket `Res.[A(r,s).B(t,u)…]`

The bracket notation chains modifications onto one residue in order. **Each step's R-group refers to the fragment immediately before it**, not to the original residue — so the string reads in the same order as the synthesis steps. Each bracket is independent, they can be mixed freely with crosslinks and branch notation, and multiple brackets can appear on the same peptide.

```
Host.[CapA(host_r, capA_r).CapB(capA_r, capB_r).CapC(capB_r, capC_r)]
```

**Two-step** — maleimide thiol-Michael then AEEA linker loading:

```python
# Step 1: Mal R4 (maleimide) reacts with Cys R4 (thiol → thioether), Mal R2 is now exposed
# Step 2: AEEA R1 (amine) occupies Mal R2 (carboxyl) — amide bond to linker, not to Cys
seq = Sequence('G-C.[Mal(4,4).AEEA(2,1)]-A')
```

**Three-step** — fatty acid lipidation (semaglutide/retatrutide-type isopeptide linker):

```python
# Step 1: gGlu R4 (γ-COOH) bonds to Lys R4 (ε-amine) — isopeptide
# Step 2: gGlu R1 (α-amine) bonds to AEEA R2 (COOH) — amide
# Step 3: AEEA R1 (amine) bonds to C20FA R2 (COOH) — amide to fatty acid
# (1,2) = bracket grows C→N from gGlu toward C20FA
seq = Sequence('ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')
```

Reading the bracket left to right matches synthesis order from Lys outward: K → gGlu (isopeptide) → AEEA → C20FA.

**Brackets compose freely with everything else:**

```python
# Single-step bracket is valid but redundant — C.[Mal(4,4)] == C.Mal(4,4)
# For Mal, R4 is the maleimide (thiol-reactive); R2 is the linker carboxyl
# Brackets only add value when chaining 2+ sequential steps
seq = Sequence('!1-C.[Mal(4,4)]-A-G-K-!1')

# Multiple protecting groups — inline notation is cleaner for single steps
seq = Sequence('fmoc-C.trt(4,2)-G-K.boc(4,2)-am')

# Bracket lipidation + disulfide staple on the same peptide
seq = Sequence('ac-C.!1(4,4)-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
```

### 4 — Sidechain branch `mainchain%%branch`

The `%%` separator attaches a pendant chain to the main chain via a crosslink. The branch string reads N→C, bonding residues through their backbone R-groups (R1→R2 at each `-`).

```python
# AEEA linker branch: Lys ε-amine (R4) → AEEA C-terminus (R2), N-capped with ac
seq = Sequence('ac-K.!1(4,2)-G-am%ac-AEEA.!1')

# Bracket equivalent (same molecule)
seq = Sequence('ac-K.[AEEA(4,2).ac(1,2)]-G-am')
```

Branch and bracket are interchangeable when the pendant chain attaches via a **backbone R-group** of the junction residue. For sidechain attachment (e.g. gGlu's γ-COOH is R4), bracket reads proximal→distal from K; branch string is reversed (distal→proximal) because C20FA — a cap with only R2 (no R1) — must lead the branch. The bracket uses `(1,2)` continuations because it grows C→N from gGlu toward C20FA (prev.R1 → this.R2):

```python
# Bracket: reads K → gGlu → AEEA → C20FA (synthesis order)
seq = Sequence('ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')

# Branch: same molecule, string reads C20FA → AEEA → gGlu (reversed)
seq = Sequence('ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1')
```

Full Retatrutide (GIP/GLP-1/glucagon triple agonist, 39 residues) — bracket notation (synthesis order):

```python
seq = Sequence('Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K-K.[AEEA(4,2).gGlu(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am')
```

Branch notation (same molecule — fatty acid chain reversed, C20FA leads because it has no R1):

```python
seq = Sequence('Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.!1(4,4)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am%C20FA-AEEA-gGlu.!1')
```

---

## Combining notations

The four extensions compose freely. Any crosslink, bracket, and branch can appear together on the same string:

```python
# Disulfide staple + 3-step lipidation on the same peptide
Sequence('ac-C.!1(4,4)-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
#              ^---disulfide crosslink---^  ^---isopeptide lipid linker---^

# Full Fmoc-SPPS protection: Pbf on Arg, Boc on Lys, disulfide between two Cys
Sequence('fmoc-R.pbf(4,2)-A-C.!1(4,4)-K.boc(4,2)-A-C.!1-am')
#               ^Arg guard  ^-disulfide-^  ^Lys guard  ^2nd Cys

# Hydrocarbon staple (RCM, i,i+4) plus isopeptide lipidation on the same helix
Sequence('ac-S5.!1(4,4)-A-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-S5.!1-am')
#              ^-------RCM staple-------^  ^--------lipid linker--------^

# Cyclic peptide with lipidation inside the ring
Sequence('!1-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1')
#          ^head-to-tail cycle^  ^lipidation bracket^
```

Bioconjugation handles (Mal, TCO, Tz, and all click-chemistry amino acids) are in the library and follow the same pattern.

---

## Chemistry types and bond validation

CABILN knows what chemistry each R-group represents. When you write a bond, the library checks whether the reaction is chemically sensible and warns you if it isn't.

### Protecting groups and caps

Standard Fmoc-SPPS protection pattern — these all validate silently:

```python
Sequence('fmoc-C.trt(4,2)-K.boc(4,2)-R.pbf(4,2)-am')
#         ^N-term  ^thiol    ^ε-amine   ^guanidinium ^C-term
```

### Crosslinks and macrolactams

```python
# Disulfide (thiol R4 ↔ thiol R4)
Sequence('ac-C.!1(4,4)-A-G-A-C.!1-am')

# Lactam staple (ε-amine R4 ↔ sidechain carboxyl R4)
Sequence('ac-K.!1(4,4)-A-A-A-D.!1-am')

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

# Two-step: maleimide thiol-Michael then AEEA linker; azide-Lys in the same chain
# AEEA attaches to Mal R2, not to Cys — bracket notation makes the order explicit
Sequence('G-C.[Mal(4,4).AEEA(2,1)]-A-G-AzK-G')
```

IEDDA (tetrazine ligation — tetrazine + TCO):

```python
# Tetrazine on Lys, trans-cyclooctene on adjacent Lys
Sequence('ac-K.Tz(4,2)-A-A-K.TCO(4,2)-am')
```

Oxime ligation (aminooxy + aldehyde):

```python
# Pre-formed oxime linkage in a modified residue
Sequence('ac-A-Aoa-G-Ald-A-am')
```

Depsipeptide (ester backbone):

```python
# Glycolic acid (Glc) forms an ester junction instead of amide
Sequence('ac-G-Glc-A-G-am')
```

Thioester:

```python
# S–C(=O) linkage — Cys thiol R4 to Asp sidechain carboxyl R4
Sequence('ac-C.!1(4,4)-A-A-D.!1-am')
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
seq = Sequence('ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am')

# Branch — same molecule, string reversed: C20FA → AEEA → gGlu ← K
seq = Sequence('ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1')
```

---

## Bicyclic peptides

Two simultaneous crosslinks on a single chain produce a bicyclic topology. CABILN represents each ring with an independent crosslink label; the notation composes without extra syntax.

**SFTI-1 (Sunflower Trypsin Inhibitor-1)** — 14-residue head-to-tail cyclic peptide with an internal disulfide, prototypical of the naturally occurring bicyclic Bowman-Birk scaffold:

```python
# !1 closes the head-to-tail backbone ring (N→C peptide bond)
# !2 closes a disulfide bridge between the two Cys residues (R4↔R4)
seq = Sequence('!1-G-R-C.!2(4,4)-T-K-S-I-P-P-I-C.!2-F-P-D-!1')
```

**Bicyclic peptides** with two crosslinks use four reactive residues — two pairs, each closed by its own bond ID. This gives two loops sharing a backbone segment:

```python
# Two independent thioether crosslinks — dialkylation bicycle
# loop1 = Cys1–ClAcAla1  |  loop2 = Cys2–ClAcAla2
seq = Sequence('ac-C.!1(4,4)-A-A-ClAcAla.!1-G-C.!2(4,4)-A-A-ClAcAla.!2-am')
```

Each `-A-A-` stretch is a user-defined diversity region.

**TBMB bicycle** — 1,3,5-tris(bromomethyl)benzene bridges three Cys thiols via thioether bonds. TBMB is registered as a crosslinker monomer with three degenerate R-groups (R4/R5/R6, all `alkyl_halide_c`):

```python
# TBMB hub on a pendant chain, three thioether crosslinks to Cys residues
# Each Cys R4 (thiol) bonds to a different TBMB arm (R4, R5, R6)
seq = Sequence('ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3')
```

---

## How chemistry-aware validation works

Every bond in CABILN is validated before assembly. The system runs in two stages:

**Stage 1 — R-group labelling.** When a monomer is registered via `pre_activate()`, three detection passes assign a chemistry type to every attachment point:

1. **Backbone detection** (graph topology) — shortest path between amino N and carboxyl C → R1 (`backbone_n`), R2 (`backbone_c`), R3 (`backbone_n_mod` if backbone N has ≥ 2 H)
2. **Sidechain scan** (24 SMARTS patterns in priority order) — first match wins per atom; multi-atom patterns protect non-attachment atoms from re-matching
3. **Second-H pass** — sidechain primary amines (`amine_primary`) get an extra slot for the second N–H → `amine_secondary`

The complete type set, in priority order:

| # | Type | SMARTS | LG | label_only | Example |
|---|------|--------|----|------------|---------|
| — | `backbone_n` | *(graph topology)* | `[H]` | — | all amino acids R1 |
| — | `backbone_c` | *(graph topology)* | `[OH]` | — | all amino acids R2 |
| — | `backbone_n_mod` | *(≥ 2 H on backbone N)* | `[H]` | — | all except Pro R3 |
| 1 | `thiol` | `[SX2H1:1]` | `[H]` | no | Cys |
| 2 | `selenol` | `[SeX2H1:1]` | `[H]` | no | Sec |
| 3 | `alkyl_halide_c` | `[CX4;!H0:1][Cl,Br,I]` | `[H]` | no | ClAc, BrAc |
| 4 | `aminooxy` | `[NH2:1][OX2H0]` | `[H]` | no | Aoa |
| 5 | `hydrazide` | `[NX3H1:1][NX3H2]` | `[H]` | no | HzAla |
| 6 | `amine_primary` | `[NX3;H2:1]` | `[H]` | no | Lys, Orn |
| 7 | `guanidinium` | `[NX3;H1:1][CX3](=N)` | `[H]` | **yes** | Arg |
| 8 | `guanidinium_imine` | `[NX2H1:1]=[CX3]([NX3])[NX3]` | `[H]` | no | Arg =NH |
| 9 | `carboxyl` | `[CX3:1](=O)[OX2H1]` | `[OH]` | no | Asp, Glu |
| 10 | `hydroxyl_phenolic` | `[OX2H1:1][c]` | `[H]` | **yes** | Tyr |
| 11 | `hydroxyl` | `[OX2H1:1][CX4]` | `[H]` | no | Ser, Thr |
| 12 | `aromatic_nh` | `[nH:1]` | `[H]` | **yes** | His, Trp |
| 13 | `amide_nh` | `[NX3;H1:1][CX3]=O` | `[H]` | **yes** | Lys_Boc, Cit |
| 14 | `phosphate_p` | `[P:1](=O)([OH])[OH]` | `[H]` | no | pSer, pTyr |
| 15 | `cyclooctyne_c` | `[CX4;!H0:1][C;r]#[C;r]` | `[H]` | no | CyoAla |
| 16 | `alkyne_c` | `[CX4;!H0:1]C#[CH]` | `[H]` | no | Pra |
| 17 | `azide_alpha_c` | `[CX4;!H0:1][N]=[N+]=[N-]` | `[H]` | no | AzAla |
| 18 | `terminal_alkene` | `[CH1:1]=[CH2]` | `[H]` | no | S5, R8 |
| 19 | `tetrazine_c` | `[CX4;!H0:1][c]1nnc(nn1)` | `[H]` | no | TzAla |
| 20 | `tco_c` | `[CX4;!H0;!r:1][C;r]=[C;r]` | `[H]` | no | TcoAla |
| 21 | `aldehyde` | `[CX3H1:1](=O)[!#7;!#1]` | `[H]` | no | Ald |
| 22 | `formamide_c` | `[CX3H1:1](=O)[#7X3]` | `[H]` | no | Lys_For |
| 23 | `nhs_ester` | `[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O` | `[H]` | no | NHSAla |
| 24 | `maleimide_c` | `[CH1:1]1=[C][C](=[O])[N][C]1=[O]` | `[H]` | no | MalAla |
| — | `amine_secondary` | *(Pass 3: second H)* | `[H]` | — | Lys R5 |

**`label_only`** types get a visible R-group slot but no entry in the bond reaction table. This prevents incorrect bonding of detectable-but-unreactive atoms (e.g. Arg's guanidinium NH, Tyr's phenolic OH, His's imidazole NH). The slot is still useful for cap attachment and for tools that inspect monomer connectivity.

**`label_only` protection rule**: multi-atom matches for `label_only` types do NOT add non-attachment atoms to the `protected` set. This allows reactive types to detect their attachment points even when a label_only pattern matches a neighboring atom. For example, Arg's guanidinium (#7, label_only) reserves the ε-NH slot without blocking guanidinium_imine (#8) from detecting the =NH. Similarly, amide_nh (#13, label_only) reserves the amide N–H without blocking formamide_c (#22) from detecting the formyl C.

Each type is documented below with 5 walkthrough examples showing step-by-step detection logic.

---

#### Backbone detection — graph topology

R1/R2 are found by the shortest path between the amino N and carboxyl C. R3 is assigned if the backbone N has ≥ 2 hydrogens. No SMARTS pattern — pure graph search.

<details>
<summary>5 examples</summary>

```python
# 1. Alanine — standard α-amino acid
pre_activate('N[C@@H](C)C(=O)O')
# Backbone path: N → Cα → C(=O) → O. R1=backbone_n, R2=backbone_c.
# Backbone N has 2 H → R3=backbone_n_mod.
# No sidechain reactive group. Result: {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'}

# 2. Proline — ring nitrogen, no R3
pre_activate('O=C(O)[C@@H]1CCCN1')
# Backbone path found through pyrrolidine ring.
# Ring N has only 1 H → R3 threshold (≥ 2 H) NOT met. No R3 assigned.
# Result: {1: 'backbone_n', 2: 'backbone_c'}

# 3. β-alanine (bAla) — β-amino acid
pre_activate('NCCC(=O)O')
# Backbone path: N → C → C → C(=O). Longer path than α-amino acid.
# Graph topology handles β-, γ-amino acids automatically.
# Result: {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod'}

# 4. Acetyl cap (ac) — no nitrogen
pre_activate('CC(=O)O')
# No amino N found → cap detection. COOH present → R2=backbone_c.
# No R1 or R3. Single-ended cap for N-terminal capping.
# Result: {2: 'backbone_c'}

# 5. Glycolic acid — depsipeptide building block
pre_activate('OCC(=O)O')
# No amino N → cap detection. Hydroxyl path to COOH.
# R2=backbone_c, R4=hydroxyl (the sidechain -OH).
# Result: {2: 'backbone_c', 4: 'hydroxyl'}
```

</details>

#### thiol — `[SX2H1:1]` (priority 1)

Sulfhydryl –SH. Highest sidechain priority. LG=`[H]` — the S–H hydrogen departs during disulfide or thioether bond formation.

<details>
<summary>5 examples</summary>

```python
# 1. Cysteine (C) — canonical thiol
pre_activate('N[C@@H](CS)C(=O)O')
# Sidechain scan: [SX2H1:1] matches the -SH sulfur (degree 2, 1 H). → R4=thiol
# CHUCKLES: [1*]N([3*])[C@@H](C[S][4*])C([2*])=O

# 2. Penicillamine (Pen) — β,β-dimethyl cysteine
pre_activate('NC(C)(CS)C(=O)O')
# Same [SX2H1:1] match on -SH. gem-dimethyl backbone doesn't affect sidechain.
# CHUCKLES: [1*]N([3*])C(C)(C[S][4*])C([2*])=O

# 3. Homocysteine (Hcy) — one extra CH₂ in sidechain
pre_activate('N[C@@H](CCS)C(=O)O')
# [SX2H1:1] matches regardless of chain length to sulfur.
# CHUCKLES: [1*]N([3*])[C@@H](CC[S][4*])C([2*])=O

# 4. D-Cysteine (DCys) — opposite chirality
pre_activate('N[C@H](CS)C(=O)O')
# SMARTS is chirality-agnostic. Same R4=thiol assignment.

# 5. Deamino-cysteine — cap with thiol sidechain
pre_activate('OC(=O)CS')
# No amino N → cap backbone. R2=backbone_c. Sidechain: R4=thiol.
# Result: {2: 'backbone_c', 4: 'thiol'}
```

</details>

#### selenol — `[SeX2H1:1]` (priority 2)

Selenohydryl –SeH. Diselenide chemistry analogue of thiol. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Selenocysteine (Sec) — the 21st amino acid
pre_activate('N[C@@H](C[SeH])C(=O)O')
# [SeX2H1:1] matches -SeH (degree 2, 1 H). → R4=selenol

# 2. D-Selenocysteine (seC) — D-form
pre_activate('N[C@H](C[SeH])C(=O)O')
# Same match, opposite chirality.

# 3. Homo-selenocysteine (Se_Hcy) — longer sidechain
pre_activate('N[C@@H](CC[SeH])C(=O)O')
# Chain length doesn't affect [SeX2H1:1] match. → R4=selenol

# 4. D-Homo-selenocysteine (D_Se_Hcy)
pre_activate('N[C@H](CC[SeH])C(=O)O')
# D-form of Se_Hcy. Same detection.

# 5. Selenopenicillamine (Se_Pen) — gem-dimethyl
pre_activate('NC(C)(C[SeH])C(=O)O')
# [SeX2H1:1] on Se, unaffected by β,β-disubstitution. → R4=selenol
```

</details>

#### alkyl_halide_c — `[CX4;!H0:1][Cl,Br,I]` (priority 3)

sp3 carbon bearing a halide. The dummy replaces a C–H; the halide stays in the monomer and departs as part of the SMIRKS reaction (SN₂). LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Chloroacetyl (ClAc) — chloroacetic acid cap
pre_activate('OC(=O)CCl')
# [CX4;!H0:1] matches the CH₂ (sp3, has H); [Cl] matches chlorine. → R4=alkyl_halide_c
# No amino N → cap backbone. Result: {2: 'backbone_c', 4: 'alkyl_halide_c'}

# 2. Bromoacetyl (BrAc) — bromine variant
pre_activate('OC(=O)CBr')
# Same SMARTS, [Br] matches. → R4=alkyl_halide_c

# 3. β-Bromoalanine (Ala_3Br) — amino acid with halide sidechain
pre_activate('N[C@@H](CBr)C(=O)O')
# Full amino acid backbone + sidechain halide. R4=alkyl_halide_c

# 4. β-Chloroalanine (Ala_3Cl) — chlorine on Ala
pre_activate('N[C@@H](CCl)C(=O)O')
# Same pattern, Cl variant. R4=alkyl_halide_c

# 5. D-β-Bromoalanine (D_Ala_3Br) — D-form
pre_activate('N[C@H](CBr)C(=O)O')
# Chirality-agnostic. R4=alkyl_halide_c
```

</details>

#### aminooxy — `[NH2:1][OX2H0]` (priority 4)

Aminooxy –O–NH₂. Checked *before* amine_primary (#6) so `-ONH₂` is not misclassified as a generic amine. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Aminooxy-serine (Aoa) — canonical aminooxy
pre_activate('N[C@@H](CON)C(=O)O')
# [NH2:1] matches the terminal -NH₂; [OX2H0] matches the O (2 bonds, no H).
# ⚠ PRIORITY: checked before amine_primary (#6), so -ONH₂ → aminooxy, not amine.
# CHUCKLES: [1*]N([3*])[C@@H](CON[4*])C([2*])=O → R4=aminooxy

# 2. Aminooxy-lysine (AoaLys) — Lys-length linker
pre_activate('N[C@@H](CCCCON)C(=O)O')
# Longer chain, same -ONH₂ terminus. → R4=aminooxy

# 3. Aminooxy-ornithine (AoaOrn) — Orn-length
pre_activate('N[C@@H](CCCON)C(=O)O')
# R4=aminooxy. Chain length irrelevant.

# 4. Aminooxy-Dab (AoaDab) — short linker
pre_activate('N[C@@H](CCON)C(=O)O')
# R4=aminooxy.

# 5. D-Aminooxy-serine (D_Aoa) — D-form
pre_activate('N[C@H](CON)C(=O)O')
# Same detection, opposite chirality. R4=aminooxy
```

</details>

#### hydrazide — `[NX3H1:1][NX3H2]` (priority 5)

Hydrazide –CO–NH–NH₂. The N–N bond distinguishes it from plain amine. The `:1` atom is the NH (not the terminal NH₂). Multi-atom match protects the terminal NH₂ from being re-classified as amine_primary. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Hydrazide-alanine (HzAla) — canonical hydrazide
pre_activate('N[C@@H](CC(=O)NN)C(=O)O')
# [NX3H1:1] matches the -NH- (1 H); [NX3H2] matches the -NH₂.
# ⚠ PROTECTION: The terminal -NH₂ (match[1]) is added to `protected` set,
#   so amine_primary (#6) cannot re-match it. → R4=hydrazide only.
# CHUCKLES: [1*]N([3*])[C@@H](CC(=O)N([4*])N)C([2*])=O

# 2. Hydrazide-lysine (HzLys) — long linker
pre_activate('N[C@@H](CCCCC(=O)NN)C(=O)O')
# Same -CONHNH₂ terminus. → R4=hydrazide

# 3. Hydrazide-ornithine (HzOrn)
pre_activate('N[C@@H](CCCC(=O)NN)C(=O)O')
# R4=hydrazide.

# 4. Hydrazide-Dab (HzDab) — short linker
pre_activate('N[C@@H](CCC(=O)NN)C(=O)O')
# R4=hydrazide.

# 5. D-Hydrazide-alanine (D_HzAla) — D-form
pre_activate('N[C@H](CC(=O)NN)C(=O)O')
# R4=hydrazide. Chirality-agnostic.
```

</details>

#### amine_primary — `[NX3;H2:1]` (priority 6)

Generic sidechain –NH₂. Catches all primary amines not already claimed by aminooxy (#4) or hydrazide (#5). In Pass 3, a second slot is added for the second N–H → `amine_secondary`. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Lysine (K) — ε-amine gets two slots
pre_activate('NCCCC[C@@H](N)C(=O)O')
# [NX3;H2:1] matches the ε-NH₂. → R4=amine_primary
# Pass 3: same atom gets R5=amine_secondary (second H).
# CHUCKLES: [1*]N([3*])[C@@H](CCCCN([4*])[5*])C([2*])=O

# 2. Ornithine (Orn) — one CH₂ shorter than Lys
pre_activate('N[C@@H](CCCN)C(=O)O')
# ε-NH₂ → R4=amine_primary, R5=amine_secondary.

# 3. 2,3-Diaminopropionic acid (Dap) — shortest amine sidechain
pre_activate('N[C@@H](CN)C(=O)O')
# β-NH₂ → R4=amine_primary, R5=amine_secondary.

# 4. D-Lysine (DLys) — D-form
pre_activate('NCCCC[C@H](N)C(=O)O')
# Same detection. R4=amine_primary, R5=amine_secondary.

# 5. Aminobutyric acid (Dab) — γ-amine
pre_activate('N[C@@H](CCN)C(=O)O')
# γ-NH₂ → R4=amine_primary, R5=amine_secondary.
```

</details>

#### guanidinium — `[NX3;H1:1][CX3](=N)` (priority 7, label_only)

Guanidinium NH on arginine. **label_only**: slot is visible but has no bond-table reaction — prevents the guanidinium N–H from participating in unintended bond formation. The multi-atom match protects the guanidinium C from later pattern re-matching.

<details>
<summary>5 examples</summary>

```python
# 1. Arginine (R) — label_only protection allows guanidinium_imine
pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# Sidechain scan (in order):
#   amine_primary (#6): terminal -NH₂ matches → R4=amine_primary
#   guanidinium (#7): the ε-NH matches [NX3;H1:1][CX3](=N) → R5=guanidinium (label_only)
#     ⚠ label_only → match tuple NOT added to `protected`.
#   guanidinium_imine (#8): =NH matches [NX2H1:1]=[CX3]... → R6=guanidinium_imine
# Pass 3: R4 amine_primary gets second H → R7=amine_secondary
# Result: {4: 'amine_primary', 5: 'guanidinium', 6: 'guanidinium_imine', 7: 'amine_secondary'}

# 2. D-Arginine (DArg)
pre_activate('N[C@H](CCCNC(=N)N)C(=O)O')
# Same detection. {4: amine_primary, 5: guanidinium, 6: guanidinium_imine, 7: amine_secondary}

# 3. Homoarginine (hArg) — one extra CH₂
pre_activate('N[C@@H](CCCCNC(=N)N)C(=O)O')
# Longer chain, same detection.

# 4. N-Methyl-arginine (meR) — slot numbering shift
pre_activate('CN[C@@H](CCCNC(=N)N)C(=O)O')
# Backbone N has 1 H → no R3. Sidechain starts at R3 (not R4).
# Result: {3: 'amine_primary', 4: 'guanidinium', 5: 'guanidinium_imine', 6: 'amine_secondary'}

# 5. Arginine aldehyde (Arg_al) — C-terminal aldehyde
pre_activate('N[C@@H](CCCNC(=N)N)C=O')
# Backbone C is aldehyde (not COOH). R2 still assigned.
# Sidechain guanidinium detection unaffected.
```

</details>

#### guanidinium_imine — `[NX2H1:1]=[CX3]([NX3])[NX3]` (priority 8)

Guanidinium =NH imine nitrogen. Matches the imine N in Arg-like guanidinium groups. Guanidinium (#7, label_only) matches the ε-NH first but does NOT protect the imine N (label_only types skip protection). So guanidinium_imine detects the =NH as a reactive slot. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Arginine (R) — all three guanidinium nitrogens detected
pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# Terminal -NH₂ → R4=amine_primary
# ε-NH → R5=guanidinium (label_only)
# =NH → R6=guanidinium_imine (reactive)
# Pass 3: R7=amine_secondary
# Result: {4: 'amine_primary', 5: 'guanidinium', 6: 'guanidinium_imine', 7: 'amine_secondary'}

# 2. D-Arginine (DArg)
pre_activate('N[C@H](CCCNC(=N)N)C(=O)O')
# Same detection, opposite chirality.

# 3. Homoarginine (hArg) — one extra CH₂
pre_activate('N[C@@H](CCCCNC(=N)N)C(=O)O')
# Same guanidinium detection regardless of chain length.

# 4. N-Methyl-arginine (meR) — slot numbering shift
pre_activate('CN[C@@H](CCCNC(=N)N)C(=O)O')
# No R3 (backbone N has 1H). Sidechain starts at R3.
# Result: {3: 'amine_primary', 4: 'guanidinium', 5: 'guanidinium_imine', 6: 'amine_secondary'}

# 5. ADMA (asymmetric dimethylarginine) — dimethylated terminal N
pre_activate('CN(C)C(=N)NCCC[C@@H](N)C(=O)O')
# Terminal N is N(CH₃)₂ — no H₂, so amine_primary doesn't match.
# ε-NH → guanidinium (label_only). =NH → guanidinium_imine.
```

</details>

#### carboxyl — `[CX3:1](=O)[OX2H1]` (priority 9)

Sidechain –COOH. LG=`[OH]` (not `[H]`) — this distinguishes carboxyl from aldehyde at assembly time. Backbone COOH is excluded by graph topology claim before sidechain scan runs.

<details>
<summary>5 examples</summary>

```python
# 1. Aspartate (D) — two COOH groups, only sidechain detected
pre_activate('N[C@@H](CC(=O)O)C(=O)O')
# ⚠ EDGE CASE: Both backbone and sidechain have COOH.
# Backbone COOH claimed first by graph topology → excluded from scan.
# Sidechain COOH matches [CX3:1](=O)[OX2H1] → R4=carboxyl, LG=[OH]

# 2. Glutamate (E) — longer sidechain COOH
pre_activate('N[C@@H](CCC(=O)O)C(=O)O')
# Same logic. Backbone COOH excluded; sidechain → R4=carboxyl.

# 3. γ-Glutamate (gGlu) — γ-linked backbone
pre_activate('N[C@@H](CCC(=O)O)CC(=O)O')
# Graph topology picks the γ-carboxyl as backbone C.
# The α-carboxyl becomes sidechain → R4=carboxyl.

# 4. D-Aspartate (DAsp)
pre_activate('N[C@H](CC(=O)O)C(=O)O')
# Same as Asp, D-chirality. R4=carboxyl.

# 5. Glutaric acid derivative (Gla) — dicarboxylic sidechain
pre_activate('N[C@@H](CC(=O)CC(=O)O)C(=O)O')
# Backbone COOH excluded. Sidechain COOH matches → R4=carboxyl.
```

</details>

#### hydroxyl_phenolic — `[OX2H1:1][c]` (priority 10, label_only)

Phenolic –OH on aromatic ring. **label_only**: phenolic OH is unreactive under standard peptide conditions. Checked before aliphatic hydroxyl (#11) so Tyr doesn't get the wrong type.

<details>
<summary>5 examples</summary>

```python
# 1. Tyrosine (Y) — canonical phenol
pre_activate('N[C@@H](Cc1ccc(O)cc1)C(=O)O')
# [OX2H1:1][c] matches the phenolic -OH on the aromatic ring. → R4=hydroxyl_phenolic
# ⚠ PRIORITY: matched before hydroxyl (#11), so [OX2H1][c] wins over [OX2H1][CX4].

# 2. D-Tyrosine (DTyr)
pre_activate('N[C@H](Cc1ccc(O)cc1)C(=O)O')
# Same detection, D-form.

# 3. 3-Iodotyrosine (Tyr_3I) — halogenated ring
pre_activate('N[C@@H](Cc1cc(I)c(O)cc1)C(=O)O')
# Phenolic -OH still on aromatic ring. Iodo substituent doesn't block.

# 4. 3-Nitrotyrosine (Tyr_3NO2) — nitrated
pre_activate('N[C@@H](Cc1cc([N+](=O)[O-])c(O)cc1)C(=O)O')
# Still phenolic OH on aromatic carbon. → R4=hydroxyl_phenolic

# 5. 2,6-Dimethyltyrosine (Tyr_26diMe) — sterically hindered
pre_activate('N[C@@H](Cc1c(C)cc(O)cc1C)C(=O)O')
# Methyl groups don't change the O-H on aromatic ring match.
```

</details>

#### hydroxyl — `[OX2H1:1][CX4]` (priority 11)

Aliphatic –OH on sp3 carbon. LG=`[H]`. Matches Ser, Thr, hydroxyproline. Does not match phenolic OH (caught by hydroxyl_phenolic at #10).

<details>
<summary>5 examples</summary>

```python
# 1. Serine (S) — canonical aliphatic hydroxyl
pre_activate('N[C@@H](CO)C(=O)O')
# [OX2H1:1][CX4] matches the -OH on sp3 CH₂. → R4=hydroxyl

# 2. Threonine (T) — β-methyl serine
pre_activate('N[C@@H]([C@@H](C)O)C(=O)O')
# Same match on -OH. Methyl branch irrelevant.

# 3. Homoserine (Hser) — one extra CH₂
pre_activate('N[C@@H](CCO)C(=O)O')
# Chain length doesn't affect [OX2H1][CX4]. → R4=hydroxyl

# 4. Hydroxyproline (Hyp) — ring hydroxyl
pre_activate('O=C(O)[C@@H]1C[C@@H](O)CN1')
# Ring -OH on sp3 C. [CX4] is the ring carbon. → R4=hydroxyl

# 5. D-Serine (DSer) — D-form
pre_activate('N[C@H](CO)C(=O)O')
# Same hydroxyl detection, opposite chirality. R4=hydroxyl.
# Note: tert-butyl ether (COC(C)(C)C) does NOT match — ether O has no H.
```

</details>

#### aromatic_nh — `[nH:1]` (priority 12, label_only)

Aromatic ring NH (imidazole, indole). **label_only**: aromatic NH is unreactive under standard conditions but the slot is useful for monitoring and cap attachment.

<details>
<summary>5 examples</summary>

```python
# 1. Histidine (H) — imidazole NH
pre_activate('N[C@@H](Cc1c[nH]cn1)C(=O)O')
# [nH:1] matches the aromatic ring nitrogen bearing H. → R4=aromatic_nh

# 2. Tryptophan (W) — indole NH
pre_activate('N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O')
# Same [nH:1] match on the indole ring nitrogen. → R4=aromatic_nh

# 3. D-Histidine (DHis)
pre_activate('N[C@H](Cc1c[nH]cn1)C(=O)O')
# D-form. Same detection.

# 4. D-Tryptophan (DTrp)
pre_activate('N[C@H](Cc1c[nH]c2ccccc12)C(=O)O')
# D-form. R4=aromatic_nh.

# 5. 5-Hydroxytryptophan (Trp_5OH) — two sidechain slots
pre_activate('N[C@@H](Cc1c[nH]c2ccc(O)cc12)C(=O)O')
# R4=aromatic_nh (indole NH), R5=hydroxyl_phenolic (5-OH on aromatic ring).
# Both detected independently — different atoms, no priority conflict.
```

</details>

#### amide_nh — `[NX3;H1:1][CX3]=O` (priority 13, label_only)

Amide N–H (sidechains of Asn, Gln; Boc protecting group; formamide). **label_only**: prevents amide NH from participating in unintended bond formation. Multi-atom match: `[CX3]=O` is added to `protected`, blocking the adjacent carbonyl C from matching later patterns (notably formamide_c #22).

<details>
<summary>5 examples</summary>

```python
# 1. Boc-lysine (Lys_Boc) — carbamate NH
pre_activate('N[C@@H](CCCCNC(=O)OC(C)(C)C)C(=O)O')
# The carbamate -NHC(=O)O- has NH (1H). [NX3;H1:1][CX3]=O matches. → R4=amide_nh
# ⚠ WHY NOT amine_primary? amine_primary (#6) requires H2, but this N has only 1H.

# 2. Citrulline (Cit) — ureido sidechain, two N–H types
pre_activate('N[C@@H](CCCNC(=O)N)C(=O)O')
# -NHC(=O)NH₂ has BOTH N-H types:
#   amine_primary (#6) matches terminal -NH₂ (H2) → R4=amine_primary
#   amide_nh (#13) matches the ε-NH (H1) bonded to C=O → R5=amide_nh (label_only)
# Pass 3: R4 gets second H → R6=amine_secondary.
# Result: {4: 'amine_primary', 5: 'amide_nh', 6: 'amine_secondary'}

# 3. Formyl-lysine (Lys_For) — formamide sidechain, two slots
pre_activate('N[C@@H](CCCCNC=O)C(=O)O')
# amide_nh matches the N-H → R4=amide_nh (label_only).
# ⚠ label_only → formamide C is NOT protected.
# formamide_c (#22) matches the C-H → R5=formamide_c (reactive).
# Result: {4: 'amide_nh', 5: 'formamide_c'}

# 4. Acetyl-lysine (Lys_Ac) — acetamide sidechain
pre_activate('N[C@@H](CCCCNC(=O)C)C(=O)O')
# -NHC(=O)CH₃: the NH matches amide_nh. → R4=amide_nh (label_only)

# 5. Asparagine (N) — sidechain amide, NOT amide_nh
pre_activate('N[C@@H](CC(=O)N)C(=O)O')
# ⚠ SURPRISE: Asn's -C(=O)NH₂ has NH₂ (2H), not NH (1H).
# amine_primary (#6) matches first: [NX3;H2:1] → R4=amine_primary
# amide_nh (#13) requires H1 — does NOT match NH₂.
# Pass 3: R5=amine_secondary. Result: {4: 'amine_primary', 5: 'amine_secondary'}
```

</details>

#### phosphate_p — `[P:1](=O)([OH])[OH]` (priority 14)

Phosphate group. The dummy is placed on phosphorus. LG=`[H]`. Matched on raw SMILES; after CHUCKLES conversion the –OH groups may not be present, so `infer_chem_type` falls back to element-based heuristic (`element_15`) at assembly time.

<details>
<summary>5 examples</summary>

```python
# 1. Phosphoserine (pSer)
pre_activate('N[C@@H](COP(=O)(O)O)C(=O)O')
# [P:1](=O)([OH])[OH] matches the phosphate P. → R4=phosphate_p

# 2. Phosphothreonine (pThr)
pre_activate('N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)O')
# Same phosphate match. → R4=phosphate_p

# 3. Phosphotyrosine (pTyr)
pre_activate('N[C@@H](Cc1ccc(OP(=O)(O)O)cc1)C(=O)O')
# Phosphate on aromatic ring via ether. → R4=phosphate_p

# 4. Ser_PO3H2 — alternative naming
pre_activate('N[C@@H](COP(=O)(O)O)C(=O)O')
# Same as pSer. Naming convention varies.

# 5. Thr_PO3H2
pre_activate('N[C@@H]([C@@H](C)OP(=O)(O)O)C(=O)O')
# Same as pThr.
```

</details>

#### cyclooctyne_c — `[CX4;!H0:1][C;r]#[C;r]` (priority 15)

sp3 C–H adjacent to a ring-strained alkyne (cyclooctyne). Two ring carbons flank the triple bond → two slots per cyclooctyne. LG=`[H]`. Used in strain-promoted azide–alkyne cycloaddition (SPAAC).

<details>
<summary>5 examples</summary>

```python
# 1. Cyclooctyne-alanine (CyoAla) — two slots from one ring
pre_activate('N[C@@H](CC1CCCCC#CC1)C(=O)O')
# The 8-membered ring has C#C. Two ring carbons adjacent to the triple bond
# each match [CX4;!H0:1][C;r]#[C;r]. → R4=cyclooctyne_c, R5=cyclooctyne_c

# 2. Cyclooctyne-lysine (CyoLys) — long linker
pre_activate('N[C@@H](CCCCC1CCCCC#CC1)C(=O)O')
# Same ring, longer sidechain. Two cyclooctyne_c slots.

# 3. Cyclooctyne-ornithine (CyoOrn) — Orn-length
pre_activate('N[C@@H](CCCC1CCCCC#CC1)C(=O)O')
# R4=cyclooctyne_c, R5=cyclooctyne_c.

# 4. Cyclooctyne-Dab (CyoDab) — short linker
pre_activate('N[C@@H](CCC1CCCCC#CC1)C(=O)O')
# Two slots regardless of linker length.

# 5. Cyclooctyne cap (Cyo) — no amino group
pre_activate('OC(=O)CC1CCCCC#CC1')
# Cap backbone: R2=backbone_c. Ring → R4=cyclooctyne_c, R5=cyclooctyne_c.
```

</details>

#### alkyne_c — `[CX4;!H0:1]C#[CH]` (priority 16)

sp3 C–H adjacent to a terminal alkyne –C≡CH. LG=`[H]`. The dummy goes on the sp3 neighbor, not the alkyne C. Used in CuAAC click chemistry.

<details>
<summary>5 examples</summary>

```python
# 1. Propargylglycine (Pra) — canonical alkyne handle
pre_activate('N[C@@H](CC#C)C(=O)O')
# [CX4;!H0:1] matches the sp3 -CH₂-; C#[CH] matches the terminal alkyne.
# → R4=alkyne_c. Dummy on sp3 C, not on the alkyne carbon.

# 2. Homopropargylglycine (Hpg) — one extra CH₂
pre_activate('N[C@@H](CCC#C)C(=O)O')
# The sp3 C adjacent to C≡CH is the last CH₂ before the triple bond.
# → R4=alkyne_c

# 3. Phe-4-O-propargyl (Phe_4OPrg) — propargyl ether on Phe
pre_activate('N[C@@H](Cc1ccc(OCC#C)cc1)C(=O)O')
# The sp3 -OCH₂- matches [CX4;!H0:1]. → R4=alkyne_c

# 4. D-Homopropargylglycine (D_HPG) — D-form
pre_activate('N[C@H](CCC#C)C(=O)O')
# R4=alkyne_c.

# 5. Cysteine-propargyl (Cys_Prg) — thioether, not thiol
pre_activate('N[C@@H](CSCC#C)C(=O)O')
# ⚠ The S is a thioether (S bonded to two C's), NOT a thiol (no S-H).
# [SX2H1] does NOT match (H0). Only alkyne_c matches. → R4=alkyne_c
```

</details>

#### azide_alpha_c — `[CX4;!H0:1][N]=[N+]=[N-]` (priority 17)

sp3 C–H bearing an azide –N₃. LG=`[H]`. Dummy on the sp3 carbon, not on the azide nitrogen. Used in CuAAC and SPAAC.

<details>
<summary>5 examples</summary>

```python
# 1. Azidoalanine (AzAla) — canonical azide handle
pre_activate('N[C@@H](CN=[N+]=[N-])C(=O)O')
# [CX4;!H0:1] on the sp3 CH₂; azide [N]=[N+]=[N-] matched. → R4=azide_alpha_c

# 2. Azidolysine (AzK) — Lys scaffold
pre_activate('N[C@@H](CCCCN=[N+]=[N-])C(=O)O')
# Terminal azide on ε-carbon. → R4=azide_alpha_c

# 3. Azido-ornithine (AzOrn) — Orn scaffold
pre_activate('N[C@@H](CCCN=[N+]=[N-])C(=O)O')
# → R4=azide_alpha_c

# 4. D-Azidoalanine (D_AzAla) — D-form
pre_activate('N[C@H](CN=[N+]=[N-])C(=O)O')
# R4=azide_alpha_c.

# 5. Azidohomoalanine (AzHal) — methionine surrogate
pre_activate('N[C@@H](CCN=[N+]=[N-])C(=O)O')
# Extra CH₂. → R4=azide_alpha_c
```

</details>

#### terminal_alkene — `[CH1:1]=[CH2]` (priority 18)

Vinyl –CH=CH₂ for ring-closing metathesis (RCM). LG=`[H]`. The atom map `:1` is on the =CH–, not on =CH₂.

<details>
<summary>5 examples</summary>

```python
# 1. (S)-pentenylglycine (S5) — i,i+4 staple
pre_activate('N[C@@H](CCC=C)C(=O)O')
# [CH1:1]=[CH2] matches the terminal vinyl. → R4=terminal_alkene

# 2. (R)-octenylglycine (R8) — i,i+7 staple
pre_activate('N[C@@H](CCCCCC=C)C(=O)O')
# Same vinyl pattern, longer chain. → R4=terminal_alkene

# 3. Allylglycine (AllGly) — shortest vinyl sidechain
pre_activate('N[C@@H](CC=C)C(=O)O')
# CH₂-CH=CH₂. → R4=terminal_alkene

# 4. D-Allylglycine (D_AllGly) — D-form
pre_activate('N[C@H](CC=C)C(=O)O')
# R4=terminal_alkene.

# 5. Gly_allyl — glycine-backbone allyl cap
pre_activate('N[C@@H](C=C)C(=O)O')
# Vinyl directly on Cα. [CH1:1]=[CH2] still matches. → R4=terminal_alkene
```

</details>

#### tetrazine_c — `[CX4;!H0:1][c]1[n][n][c][n][n]1` (priority 19)

sp3 C–H adjacent to an s-tetrazine ring. LG=`[H]`. Used in inverse electron-demand Diels–Alder (IEDDA) with TCO.

<details>
<summary>5 examples</summary>

```python
# 1. Tetrazine-alanine (TzAla) — short linker
pre_activate('N[C@H](Cc1nncnn1)C(=O)O')
# [CX4;!H0:1] matches the CH₂ adjacent to the tetrazine ring.
# [c]1nnc(nn1) matches the aromatic tetrazine. → R4=tetrazine_c

# 2. Tetrazine-lysine (TzLys) — Lys-length linker
pre_activate('N[C@@H](CCCCCc1nncnn1)C(=O)O')
# Last sp3 CH₂ before the ring. → R4=tetrazine_c

# 3. Tetrazine cap (Tz) — no amino group
pre_activate('OC(=O)Cc1nncnn1')
# Cap backbone: R2=backbone_c. → R4=tetrazine_c

# 4. Tetrazine-ornithine (TzOrn) — Orn-length
pre_activate('N[C@@H](CCCCc1nncnn1)C(=O)O')
# → R4=tetrazine_c

# 5. Tetrazine-Dab (TzDab) — short linker
pre_activate('N[C@@H](CCCc1nncnn1)C(=O)O')
# → R4=tetrazine_c
```

</details>

#### tco_c — `[CX4;!H0;!r:1][C;r]=[C;r]` (priority 20)

Exocyclic sp3 C–H adjacent to a trans-cyclooctene ring C=C. LG=`[H]`. The `!r` constraint means the matched C must NOT be in the ring. Used in IEDDA with tetrazine.

<details>
<summary>5 examples</summary>

```python
# 1. TCO-alanine (TcoAla) — short linker
pre_activate('N[C@H](CC1=CCCCCCC1)C(=O)O')
# [CX4;!H0;!r:1] matches the exocyclic sp3 CH adjacent to ring C=C.
# [C;r]=[C;r] matches the ring double bond. → R4=tco_c

# 2. TCO-lysine (TcoLys) — Lys-length
pre_activate('N[C@@H](CCCCC1=CCCCCCC1)C(=O)O')
# Last exocyclic sp3 CH₂ before ring. → R4=tco_c

# 3. TCO cap (TCO) — no amino group
pre_activate('OC(=O)CC1=CCCCCCC1')
# Cap backbone. → R4=tco_c

# 4. TCO-ornithine (TcoOrn) — Orn-length
pre_activate('N[C@@H](CCCC1=CCCCCCC1)C(=O)O')
# → R4=tco_c

# 5. TCO-Dab (TcoDab) — short linker
pre_activate('N[C@@H](CCC1=CCCCCCC1)C(=O)O')
# → R4=tco_c
```

</details>

#### aldehyde — `[CX3H1:1](=O)[!#7;!#1]` (priority 21)

Aldehyde –CHO. LG=`[H]`. The `[!#7;!#1]` exclusion prevents matching formamide C–H (bonded to nitrogen) and ensures H atoms don't falsely satisfy the constraint when explicit. At assembly time, carboxyl and aldehyde are distinguished by their leaving groups: `[OH]` → carboxyl, `[H]` → aldehyde.

<details>
<summary>5 examples</summary>

```python
# 1. 4-Formyl-phenylalanine (Ald) — aromatic aldehyde handle
pre_activate('N[C@@H](Cc1ccc(C=O)cc1)C(=O)O')
# [CX3H1:1](=O) matches the aldehyde C-H; [!#7;!#1] confirms neighbor is not N or H.
# → R4=aldehyde

# 2. Alanine-aldehyde (Ala_al) — C-terminal aldehyde mimic
pre_activate('N[C@@H](C)C=O')
# Backbone C is now an aldehyde (not COOH). The graph-topology still finds it.
# The aldehyde C matches [CX3H1](=O)[!#7]. → R2=backbone_c (aldehyde-type)

# 3. Aspartate-aldehyde (Asp_al) — sidechain aldehyde
pre_activate('N[C@@H](CC=O)C(=O)O')
# Sidechain -CH=O matches. → R4=aldehyde

# 4. Glycine-aldehyde (Gly_al) — simplest aldehyde
pre_activate('NCC=O')
# Backbone C as aldehyde.

# 5. Arginine-aldehyde (Arg_al) — aldehyde + guanidinium
pre_activate('N[C@@H](CCCNC(=N)N)C=O')
# ⚠ TWO SYSTEMS: sidechain has guanidinium types; backbone C is aldehyde.
# Both detected independently — sidechain scan is separate from backbone.
```

</details>

#### formamide_c — `[CX3H1:1](=O)[#7X3]` (priority 22)

Formamide C–H. Matches the carbonyl carbon in –N–CHO groups. The `[#7X3]` matches any nitrogen (aromatic or aliphatic) with degree 3. LG=`[H]`. Amide_nh (#13, label_only) matches the adjacent N first but does NOT protect the carbonyl C (label_only types skip protection), so formamide_c correctly detects the reactive electrophilic center.

<details>
<summary>5 examples</summary>

```python
# 1. Formyl-lysine (Lys_For) — both amide_nh and formamide_c detected
pre_activate('N[C@@H](CCCCNC=O)C(=O)O')
# Amide_nh (#13, label_only) matches the N-H → R4=amide_nh
#   ⚠ label_only → does NOT protect the adjacent C
# Formamide_c (#22) matches the C-H → R5=formamide_c (reactive)
# Result: {4: 'amide_nh', 5: 'formamide_c'}

# 2. Formyl-ornithine (Orn_For)
pre_activate('N[C@@H](CCCNC=O)C(=O)O')
# Same dual detection. {4: 'amide_nh', 5: 'formamide_c'}

# 3. Formyl-Dab (Dab_For) — short linker
pre_activate('N[C@@H](CCNC=O)C(=O)O')
# {4: 'amide_nh', 5: 'formamide_c'}

# 4. D-Trp-formyl (D_Trp_For) — formyl on aromatic indole N
pre_activate('N[C@H](Cc1cn(C=O)c2ccccc12)C(=O)O')
# The indole N is aromatic — [#7X3] matches it (unlike [NX3] which is aliphatic-only).
# No amide_nh here (indole N has H0). → R4=formamide_c only.

# 5. Formyl-Dap (Dap_For) — shortest formamide sidechain
pre_activate('N[C@@H](CNC=O)C(=O)O')
# {4: 'amide_nh', 5: 'formamide_c'}
```

</details>

#### nhs_ester — `[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O` (priority 23)

NHS-activated ester. The entire succinimide ring is encoded in the SMARTS. Dummy on the sp3 C adjacent to the ester carbonyl. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. NHS-alanine (NHSAla) — canonical NHS handle
pre_activate('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O')
# The entire SMARTS matches as one unit. [CX4;!H0:1] on the sp3 CH₂.
# → R4=nhs_ester

# 2. NHS-lysine (NHSLys) — Lys-length linker
pre_activate('N[C@@H](CCCCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# Longer chain, same NHS terminus. → R4=nhs_ester

# 3. NHS-ornithine (NHSOrn) — Orn-length
pre_activate('N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester

# 4. NHS-Dab (NHSDab) — short linker
pre_activate('N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester

# 5. D-NHS-alanine (D_NHSAla) — D-form
pre_activate('N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O')
# → R4=nhs_ester
```

</details>

#### maleimide_c — `[CH1:1]1=[C][C](=[O])[N][C]1=[O]` (priority 24)

Maleimide ring C=C hydrogen. LG=`[H]`. Used in thiol–maleimide conjugation.

<details>
<summary>5 examples</summary>

```python
# 1. Maleimide-alanine (MalAla) — short linker
pre_activate('N[C@@H](CN1C(=O)C=CC1=O)C(=O)O')
# The maleimide ring pattern matches the =CH-. → R4=maleimide_c

# 2. Maleimide-lysine (MalLys) — Lys-length
pre_activate('N[C@@H](CCCCN1C(=O)C=CC1=O)C(=O)O')
# Same ring pattern. → R4=maleimide_c

# 3. Maleimide cap (Mal) — no amino group
pre_activate('OC(=O)CCN1C(=O)C=CC1=O')
# Cap backbone. → R4=maleimide_c

# 4. Maleimide-ornithine (MalOrn) — Orn-length
pre_activate('N[C@@H](CCCN1C(=O)C=CC1=O)C(=O)O')
# → R4=maleimide_c

# 5. Maleimide-Dab (MalDab) — short linker
pre_activate('N[C@@H](CCN1C(=O)C=CC1=O)C(=O)O')
# → R4=maleimide_c
```

</details>

#### amine_secondary — Pass 3 (after registry scan)

Not in the SMARTS registry. For every sidechain atom already assigned `amine_primary`, if the nitrogen had 2 hydrogens, a second slot is created for the second N–H. LG=`[H]`.

<details>
<summary>5 examples</summary>

```python
# 1. Lysine (K) — R4 + R5 on same nitrogen
pre_activate('NCCCC[C@@H](N)C(=O)O')
# Pass 2: ε-NH₂ → R4=amine_primary
# Pass 3: same ε-N had 2H, already amine_primary → R5=amine_secondary

# 2. Ornithine (Orn) — δ-amine
pre_activate('N[C@@H](CCCN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 3. Diaminopropionic acid (Dap) — β-amine
pre_activate('N[C@@H](CN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 4. Diaminobutyric acid (Dab) — γ-amine
pre_activate('N[C@@H](CCN)C(=O)O')
# R4=amine_primary, R5=amine_secondary

# 5. Arginine (R) — second H on the guanidinium -NH₂
pre_activate('N[C@@H](CCCNC(=N)N)C(=O)O')
# The terminal -NH₂ got R4=amine_primary in Pass 2.
# Pass 3: same atom gets R7=amine_secondary (R5=guanidinium, R6=guanidinium_imine).
```

</details>

**Stage 2 — bond validation.** When a bond is declared (crosslink, bracket, or inline cap), the pair of chemistry types is looked up in a reaction table. Reactions are classified as:

- **Silent** — standard peptide chemistry, no warning
- **Warned** — unusual but valid, `UserWarning` emitted
- **Rejected** — raises `ValueError`

### Supported reaction types

All 19 bond-forming reactions in the YAML library (plus 6 terminal-restoration and 100 cap reactions applied automatically), with 5 CABILN examples each:

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

# 3. Macrolactam — Lys ε-amine (R4) to Asp sidechain carboxyl (R4)
Sequence('ac-K.!1(4,4)-A-A-A-D.!1-am')

# 4. Glutamate lactam — six-membered ring
Sequence('ac-K.!1(4,4)-A-A-A-A-E.!1-am')

# 5. Lactam inside head-to-tail ring — bicyclic scaffold
Sequence('!1-K.!2(4,4)-A-G-A-D.!2-!1')
```

</details>

#### sidechain_amide — amide from primary amine + sidechain carboxyl

Reactant pair: `amine_primary` + `carboxyl`

Lys ε-amine (R4) attacks Asp/Glu sidechain carboxyl (R4). Same amide bond as backbone, but both partners are sidechains — gives isopeptide (lactam) bridges.

<details>
<summary>5 examples</summary>

```python
# 1. Lys ε-amine to Asp sidechain — isopeptide bridge
Sequence('ac-K.!1(4,4)-A-A-D.!1-am')

# 2. Lys to Glu — six-membered lactam ring
Sequence('ac-K.!1(4,4)-A-A-A-E.!1-am')

# 3. Short isopeptide loop — Lys adjacent to Asp
Sequence('ac-K.!1(4,4)-G-D.!1-am')

# 4. Two isopeptide bridges — bis-lactam bicycle
Sequence('ac-K.!1(4,4)-A-D.!1-G-K.!2(4,4)-A-E.!2-am')

# 5. Isopeptide inside head-to-tail ring
Sequence('!1-K.!2(4,4)-A-G-D.!2-!1')
```

</details>

#### backbone_ester — depsipeptide ester bond

Reactant pair: `backbone_o` + `backbone_c`

An α-hydroxy acid carries its R2 as `backbone_o`; the `-` junction with the preceding residue forms an ester instead of an amide. Glc (glycolic acid) is in the library.

<details>
<summary>5 examples</summary>

```python
# 1. Glycolic acid ester junction in a tetrapeptide
Sequence('ac-G-Glc-A-G-am')

# 2. Two glycolic acid junctions — alternating ester/amide backbone
Sequence('ac-A-Glc-A-Glc-A-am')

# 3. Glycolic acid flanked by charged residues
Sequence('ac-K-Glc-D-am')

# 4. Depsipeptide + disulfide staple
Sequence('ac-C.!1(4,4)-Glc-G-A-C.!1-am')

# 5. Homoserine in a peptide chain (sidechain hydroxyl at R4)
Sequence('ac-A-Hser-G-am')
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
Sequence('ac-C.!1(4,4)-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-C.!1-am')
```

</details>

#### diselenide — diselenide bridge (Se–Se)

Reactant pair: `selenol` + `selenol`

Selenocysteine (Sec) pairs form Se–Se bridges — redox-active analogues of disulfides with lower reduction potential.

<details>
<summary>5 examples</summary>

```python
# 1. Classic diselenide — Sec pair
Sequence('ac-Sec.!1(4,4)-A-G-Sec.!1-am')

# 2. Short loop — adjacent selenocysteines
Sequence('ac-Sec.!1(4,4)-G-Sec.!1-am')

# 3. Diselenide inside head-to-tail ring
Sequence('!1-Sec.!2(4,4)-A-G-A-Sec.!2-!1')

# 4. Two diselenide bridges
Sequence('ac-Sec.!1(4,4)-A-Sec.!2(4,4)-G-Sec.!2-A-Sec.!1-am')

# 5. Mixed disulfide + diselenide on same peptide
Sequence('ac-C.!1(4,4)-A-Sec.!2(4,4)-G-C.!1-Sec.!2-am')
```

</details>

#### thioester — thioester crosslink (S–C=O)

Reactant pair: `thiol` + `backbone_c` · `thiol` + `carboxyl`

Sulphur displaces the amide nitrogen. Used in native chemical ligation models and cyclic depsipeptide analogues.

<details>
<summary>5 examples</summary>

```python
# 1. Cys thiol to Asp sidechain carboxyl — short ring
Sequence('ac-C.!1(4,4)-G-D.!1-am')

# 2. Cys thiol to Asp sidechain carboxyl — medium ring
Sequence('ac-C.!1(4,4)-A-A-D.!1-am')

# 3. Cys thiol to Glu sidechain — six-membered thioester ring
Sequence('ac-C.!1(4,4)-G-G-E.!1-am')

# 4. Longer thioester macrocycle via Asp
Sequence('ac-C.!1(4,4)-A-A-A-D.!1-am')

# 5. Thioester inside head-to-tail cyclic peptide
Sequence('!1-C.!2(4,4)-A-G-D.!2-!1')
```

</details>

#### thioether_halide — thioether via SN2 on alkyl halide

Reactant pair: `thiol` + `alkyl_halide_c`

Cys thiol displaces a halide leaving group. `ClAcAla` (chloroacetyl-alanine) is in the library with `alkyl_halide_c` at R4.

<details>
<summary>5 examples</summary>

```python
# 1. Standard Cys alkylation by chloroacetyl handle (ClAcAla is in the library)
Sequence('ac-C.!1(4,4)-A-A-ClAcAla.!1-am')

# 2. Short loop — adjacent Cys and halide
Sequence('ac-C.!1(4,4)-G-ClAcAla.!1-am')

# 3. Thioether inside head-to-tail ring
Sequence('!1-C.!2(4,4)-A-G-ClAcAla.!2-!1')

# 4. Two thioether crosslinks — dialkylation bicycle
Sequence('ac-C.!1(4,4)-A-ClAcAla.!1-G-C.!2(4,4)-A-ClAcAla.!2-am')

# 5. Thioether + lipidation
Sequence('ac-C.!1(4,4)-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-ClAcAla.!1-am')
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

# 3. Inline notation — equivalent, cleaner for a single step
Sequence('ac-A-C.MalAla(4,4)-G-am')

# 4. Two thiol-maleimide crosslinks — parallel bridges
Sequence('ac-C.!1(4,4)-A-MalAla.!1-G-C.!2(4,4)-A-MalAla.!2-am')

# 5. Thiol-maleimide inside head-to-tail ring
Sequence('!1-C.!2(4,4)-A-G-MalAla.!2-A-!1')
```

</details>

#### nhs_ester_amide — NHS ester + amine → amide (NHS departs)

Reactant pair: `amine_primary` + `nhs_ester` · `amine_secondary` + `nhs_ester`

NHS-activated ester reacts with primary or secondary amine; NHS ring leaves as byproduct. `NHSAla` (NHS-alanine) is in the library with `nhs_ester` at R4.

<details>
<summary>5 examples</summary>

```python
# 1. Lys ε-amine + NHSAla — amide ligation
Sequence('ac-K.!1(4,4)-G-A-NHSAla.!1-am')

# 2. Short loop — Lys adjacent to NHS ester
Sequence('ac-K.!1(4,4)-G-NHSAla.!1-am')

# 3. Two NHS ligations — dual Lys modification
Sequence('ac-K.!1(4,4)-A-K.!2(4,4)-G-NHSAla.!1-A-NHSAla.!2-am')

# 4. NHS amide inside head-to-tail ring
Sequence('!1-K.!2(4,4)-A-G-NHSAla.!2-!1')

# 5. NHS amide + lipidation on same peptide
Sequence('ac-K.!1(4,4)-A-NHSAla.!1-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am')
```

</details>

#### oxime_ligation — aminooxy + aldehyde → oxime

Reactant pair: `aminooxy` + `aldehyde`

Aminooxy condenses with aldehyde to form a stable oxime (C=N–O); water is implicit. Bioorthogonal under mildly acidic conditions. `Aoa` (aminooxy-alanine) and `Ald` (aldehyde-alanine) are in the library.

<details>
<summary>5 examples</summary>

```python
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

Hydrazide condenses with aldehyde forming a C=N–N hydrazone. Reversible at low pH — useful for stimuli-responsive linkers. `HzAla` (hydrazide-alanine) and `Ald` (aldehyde-alanine) are in the library.

<details>
<summary>5 examples</summary>

```python
# 1. HzAla + Ald — hydrazone bridge
Sequence('ac-HzAla.!1(4,4)-A-G-A-Ald.!1-am')

# 2. Shorter loop
Sequence('ac-HzAla.!1(4,4)-G-Ald.!1-am')

# 3. Hydrazone inside head-to-tail ring
Sequence('!1-HzAla.!2(4,4)-A-G-Ald.!2-!1')

# 4. Two hydrazone bonds — bifunctional staple
Sequence('ac-HzAla.!1(4,4)-A-Ald.!1-G-HzAla.!2(4,4)-A-Ald.!2-am')

# 5. Hydrazone + disulfide combination
Sequence('ac-C.!1(4,4)-A-HzAla.!2(4,4)-G-C.!1-Ald.!2-am')
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
Sequence('ac-Pra.!1(4,4)-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-AzK.!1-am')
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
Sequence('ac-S5.!1(4,4)-A-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-S5.!1-am')

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

#### backbone_n_alkylation — N-alkylation of backbone amide nitrogen

Reactant pair: `backbone_n_mod` (R3) + `carbon` (R1)

N-alkylation of the backbone amide NH via the R3 slot. `meA` (N-methyl cap) attaches to any non-Pro residue's R3, converting the backbone N to an N-methyl amide.

<details>
<summary>5 examples</summary>

```python
# 1. N-methyl Ala — meA cap on backbone N
Sequence('ac-A.meA(3,2)-G-am')

# 2. N-methyl Gly (sarcosine equivalent)
Sequence('ac-G.meA(3,2)-A-am')

# 3. Two N-methylations — peptoid-like backbone
Sequence('ac-A.meA(3,2)-A.meA(3,2)-G-am')

# 4. N-methylation inside head-to-tail ring
Sequence('!1-A.meA(3,2)-G-A-!1')

# 5. N-methylation + disulfide on same peptide
Sequence('ac-C.!1(4,4)-A.meA(3,2)-G-C.!1-am')
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

When you write `C.trt(4,2)`, you are saying: form a bond between **Cys R4** (the thiol) and **trt R2** (the attachment point of the trityl group). The library looks up both monomers, finds their chemical types, and validates that thiol–C is a sensible bond.

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
report = Sequence.validate('ac-K.!1(4,4)-A-A-A-D.!1-am')
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

### SMARTS matching reference

The three-pass detection pipeline (backbone → sidechain SMARTS → second-H) is documented with step-by-step walkthrough examples in [**Stage 1 — R-group labelling**](#how-chemistry-aware-validation-works) above. That section covers all 24 SMARTS patterns, priority ordering, multi-atom protection, and label_only types.

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
pyPept-monomer-add --from-cabiln "K.boc(4,2)" --symbol Lys_Boc --name "Boc-lysine"
```

---

## Legacy BILN crosslink notation

Old BILN used bare integers: `C(1,3)-A-A-A-C(1,3)`. The `1` was the bond ID and the `3` was the R-group with no delimiter — scanning left to right you couldn't tell which was which.

CABILN separates them with `.!`:

| Old BILN | CABILN equivalent |
|----------|-------------------|
| `C(1,3)-A-C(1,3)` | `C.!1(4,4)-A-C.!1` |
| `K(2,3)-A-D(2,3)` | `K.!2(4,4)-A-D.!2` |

Old BILN is accepted, but you must declare it explicitly — `Sequence()` without a format flag raises `ValueError` to avoid silent misparse:

```python
from pyPept.sequence import Sequence, biln_to_cabiln

# Raises ValueError — old notation must be declared
Sequence('C(1,3)-A-A-A-C(1,3)')

# Pass fmt='biln' to auto-convert and continue
seq = Sequence('C(1,3)-A-A-A-C(1,3)', fmt='biln')

# Or convert upfront and use CABILN directly
cabiln = biln_to_cabiln('C(1,3)-A-A-A-C(1,3)')
seq = Sequence(cabiln)   # C.!1(4,4)-A-A-A-C.!1
```

---

## Compatibility

| Input | Status | Notes |
|-------|--------|-------|
| Linear BILN `A-G-K` | Works | Unchanged from upstream |
| Capped BILN `fmoc-A-G-K-am` | Works | Unchanged |
| CABILN crosslinks `.!n(r,s)` | Works | New in this fork |
| CABILN brackets `.[A(r,s).B(t,u)]` | Works | New in this fork |
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

572 tests (571 pass, 1 xfail) covering: bond validation, inline caps, bracket notation (1/2/3-step), crosslinks, RCM staples, SPAAC/IEDDA/oxime/hydrazone/depsipeptide/thioester chemistry, fatty acid branching, phosphopeptides, monomer pipeline, full library round-trip (1029 monomers), CLI, and HELM/BILN/FASTA round-trips.

---

## References

- [pyPept: a python library to generate atomistic 2D and 3D representations of peptides](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-023-00748-2), *Journal of Cheminformatics*, 2023.
- [BILN — A Human-readable Line Notation for Complex Peptides](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00703), *J. Chem. Inf. Model.*, 2022.
