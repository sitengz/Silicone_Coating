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
V22_no_oil/
├── data.V22_no_oil
├── in.V22_no_oil
├── submit.V22_no_oil.sh
└── V22_no_oil.info
```

Every run creates a case-named folder in the current directory. To submit the
generated case:

```bash
cd V22_no_oil
sbatch submit.V22_no_oil.sh
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

Bond type 2 is reserved exclusively for formulation crosslinks. Moderator
center–arm bonds use the ordinary DMS bond type 1. Oil bonds use type 1 for
DMS/mixed backbone, type 3 for MPS–MPS backbone, and type 4 for MPS pendant
bonds.

All systems use the complete Oil-generator pair matrices and:

```lammps
special_bonds lj 0 0 0.5
```

Every DMS–MPS cross interaction (types 1–3 with types 4–5) uses:

```text
epsilon_ij = 0.579966 * sqrt(epsilon_ii * epsilon_jj)
sigma_ij   = (sigma_ii + sigma_jj) / 2
```

The epsilon prefactor is applied at both 800 K and 300 K. Pure DMS and pure
PMPS terms are unchanged.

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
| `--oil-min-separation X` | Initial oil/moderator-to-other minimum distance |
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
| `--output FILE` | Override the data filename and generated case-folder name |
| `--help` | Print built-in help |

`M2` is always calculated from `M2=2*M1/functionality`. Oil chain count is
calculated from the requested weight percentage using the actual PDMS/PMPS
masses. The realized composition and topology are recorded in the `.info`
manifest.

## Initial component placement

- Component 1 is placed from the bottom of the box upward.
- Component 2 is placed from the top of the box downward.
- Component 3 oil is confined to the central 40% of `Lz`, from `-0.20 Lz` to
  `+0.20 Lz`.
- Component 4 moderators are placed in an upper-middle band, from `+0.28 Lz`
  to `+0.38 Lz`.

Oil placement rejects beads closer than `--oil-min-separation` to components 1
and 2. Moderator placement applies the same rejection distance against all
previously placed beads, including the oil and earlier moderators.

## MSD production trajectory

After the seven-million-step equilibration workflow, the generated input runs
an additional one million steps at 300 K under NVT. It writes 1,001 frames to
`dump.msd.lammpstrj`: the timestep-zero frame followed by one frame every
1,000 steps through timestep 1,000,000. Each frame contains
`id mol type x y z ix iy iz`; the wrapped coordinates and image flags can be
combined during MSD analysis.

For example, `--output results/data.test_case` creates
`results/test_case/` and puts all four generated files inside it.

## Generated file descriptions

- `data.<case>` is the initial LAMMPS data file containing the simulation box,
  atoms, bonds, angles, and dihedrals.
- `in.<case>` is the complete LAMMPS workflow for 800 K relaxation and
  crosslinking, cooling, 300 K equilibration, and the final MSD trajectory.
  V22 uses a bond-creation probability of `0.1`.
- `submit.<case>.sh` is the one-node, 96-task Slurm submission script.
- `<case>.info` is a JSON manifest containing composition, geometry,
  force-field, topology, random-seed, and simulation settings.
