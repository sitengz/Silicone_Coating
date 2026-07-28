# V22 formulation generator with explicit silicone oil

`v22_generator.cpp` uses the shared formulation/oil implementation and the V22
defaults:

| Component | Default |
|---|---|
| Network strands | `N1=128`, `M1=900` |
| Cross-linkers | `N2=32`, functionality 8, stoichiometric `M2=225` |
| Silicone oil | None: `N3=0`, `M3=0` |
| Star moderators | `N4=5`, `M4=6` |

## Compile

From the repository root:

```bash
cd V22
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    v22_generator.cpp -o v22_generator
```

## Default no-oil reference

```bash
./v22_generator
```

```text
data.V22_no_oil
in.V22_no_oil
submit.V22_no_oil.sh
V22_no_oil.info
```

## Explicit oil examples

```bash
# PDMS N32 at 10 wt%
./v22_generator \
    --oil pdms \
    --oil-length 32 \
    --oil-wt 10
```

```bash
# PMPS N16 at 10 wt%
./v22_generator \
    --oil pmps \
    --oil-length 16 \
    --oil-wt 10
```

```bash
# Random copolymer, 32 repeats, 50% MPS positions, 10 wt% total oil
./v22_generator \
    --oil copolymer \
    --oil-length 32 \
    --oil-wt 10 \
    --mps-percent 50 \
    --sequence random
```

For copolymers, `--mps-wt X` may replace `--mps-percent X`, and the sequence
may be `random`, `alternating`, or `block`.

Omitting `--oil` always produces no component 3. Legacy `--n3`, `--m3`,
`--filler-length`, and `--filler-wt` options are rejected.

## Film example

```bash
./v22_generator \
    --thickness 100 \
    --oil pmps \
    --oil-length 16 \
    --oil-wt 10
```

Film mode preserves fixed `Lz`, uses periodic x/y boundaries, applies repulsive
z walls with separate `zlo_wall` and `zhi_wall` fixes, and controls pressure
only in x/y after compression.

## Force-field namespace

| Atom type | Meaning |
|---:|---|
| 1 | Neutral DMS |
| 2 | Reactive strand ends and moderator arms |
| 3 | Reactive cross-linker sites |
| 4 | MPS backbone |
| 5 | MPS pendant |

Bond type 2 remains reserved for formulation crosslinks and moderator bonds.
Oil bonds use type 1 for DMS/mixed backbone, type 3 for MPS–MPS backbone, and
type 4 for MPS pendant bonds.

All systems use the complete Oil-generator pair matrices and:

```lammps
special_bonds lj 0 0 0.5
```

This includes the no-oil reference so comparisons use the same global 1–4 LJ
convention. This is an intentional change from older formulation inputs.

## Command reference

| Option | Meaning |
|---|---|
| `--oil pdms|pmps|copolymer` | Explicit component-3 oil model |
| `--oil-length N` | Repeat units per chain |
| `--oil-wt X` | Oil percentage of total formulation mass |
| `--mps-percent X` | MPS percentage of copolymer repeat positions |
| `--mps-wt X` | MPS repeat-unit weight percentage |
| `--sequence MODE` | `random`, `alternating`, or `block` |
| `--oil-seed N` | Oil sequence/conformation/placement seed |
| `--oil-min-separation X` | Initial oil-to-other minimum distance |
| `--n1 N`, `--m1 M` | Network-strand length and count |
| `--n2 N` | Cross-linker length |
| `--n4 N`, `--m4 M` | Moderator size and count |
| `--functionality F` | Cross-linker functionality, 3–16 |
| `--crosslink-distribution MODE` | `random` or `regular` |
| `--crosslink-seed N` | Cross-linker site seed |
| `--mass X` | Common DMS bead mass |
| `--density X` | Initial mass density |
| `--target-density X` | Compression target density |
| `--bond-length X` | Initial DMS bond length |
| `--spacing X` | Formulation placement spacing |
| `--thickness X` | Fixed film thickness |
| `--seed N` | Moderator seed |
| `--output FILE` | Override generated case name |
| `--help` | Print built-in help |

`M2` is always calculated from `M2=2*M1/functionality`. Oil chain count is
calculated from the requested weight percentage using the actual PDMS/PMPS
masses. The realized composition and topology are recorded in the `.info`
manifest.
