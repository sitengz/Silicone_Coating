# V35 formulation generator with explicit silicone oil

`v35_generator.cpp` generates the V35 network formulation as a LAMMPS data
file, matching equilibration input, Slurm submission script, and JSON `.info`
manifest.

The default is now a **no-oil V35 reference**. Component 3 is added only when
the oil model, chain length, and weight percentage are explicitly supplied.

## Compile

From the repository root:

```bash
cd V35
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    v35_generator.cpp -o v35_generator
```

## V35 defaults

| Component | Default |
|---|---|
| Network strands | `N1=384`, `M1=306`, reactive type-2 ends |
| Cross-linkers | `N2=32`, functionality 4, stoichiometric `M2=153` |
| Silicone oil | None: `N3=0`, `M3=0` |
| Star moderators | `N4=5`, `M4=6` |
| Geometry | Cubic bulk |
| Initial density | `0.1 g/cm3` |
| Compression target | `0.8 g/cm3` |
| 1–4 LJ scaling | `special_bonds lj 0 0 0.5` |

Generate the default:

```bash
./v35_generator
```

It writes:

```text
data.V35_no_oil
in.V35_no_oil
submit.V35_no_oil.sh
V35_no_oil.info
```

## Explicit oil commands

PDMS:

```bash
./v35_generator \
    --oil pdms \
    --oil-length 32 \
    --oil-wt 10
```

PMPS:

```bash
./v35_generator \
    --oil pmps \
    --oil-length 16 \
    --oil-wt 10
```

Random 50/50 monomer copolymer:

```bash
./v35_generator \
    --oil copolymer \
    --oil-length 32 \
    --oil-wt 10 \
    --mps-percent 50 \
    --sequence random
```

Copolymer composition can instead be requested by MPS repeat-unit weight:

```bash
./v35_generator \
    --oil copolymer \
    --oil-length 32 \
    --oil-wt 10 \
    --mps-wt 50 \
    --sequence block
```

The program converts the requested composition to an integer number of MPS
repeats per chain and converts oil loading to an integer number of chains. It
prints and records the realized values.

## Oil rules

- `--oil`, `--oil-length`, and `--oil-wt` are all required to add component 3.
- Omitting `--oil` produces no oil.
- `pdms` uses only type-1 DMS beads.
- `pmps` uses type-4 backbone and type-5 pendant beads.
- `copolymer` requires exactly one of `--mps-percent` or `--mps-wt`.
- Copolymer sequences may be `random`, `alternating`, or `block`.
- Every oil chain is placed wholly inside the primary box.
- Oil chains initially occupy the region from the z midplane toward positive z.
- Oil chain count is calculated from real molecular masses, not bead count.

Legacy `--n3`, `--m3`, `--filler-length`, and `--filler-wt` inputs are rejected
to prevent an unspecified filler from being generated accidentally.

## Atom and bonded types

| Atom type | Meaning | Mass (g/mol) |
|---:|---|---:|
| 1 | Neutral DMS | 74.0 by default |
| 2 | Reactive DMS strand ends and moderator arms | 74.0 |
| 3 | Reactive DMS cross-linker sites | 74.0 |
| 4 | MPS backbone | 59.1204 |
| 5 | MPS pendant | 77.106 |

| Bond type | Meaning |
|---:|---|
| 1 | DMS or mixed DMS–MPS backbone |
| 2 | Formulation crosslink and moderator bond |
| 3 | MPS–MPS backbone |
| 4 | MPS backbone–pendant |

Angles use type 1 DMS harmonic, type 2 MPS-backbone quartic, and type 3
pendant-containing quartic interactions. Dihedral types 1–4 retain the rough
oil rules used by the standalone Oil generator.

The complete 800 K and 300 K pair matrices are written explicitly. Types 1–3
share the DMS nonbonded model. DMS–MPS interactions use:

```text
epsilon_ij = sqrt(epsilon_ii * epsilon_jj)
sigma_ij   = (sigma_ii + sigma_jj) / 2
```

## Important 1–4 convention

Every generated system, including the no-oil reference, uses:

```lammps
special_bonds lj 0 0 0.5
```

This matches the tested standalone Oil model and ensures oil/no-oil comparisons
use one global convention. It differs from older V22/V35 inputs, which relied
on the LAMMPS default of zero 1–4 LJ weight.

## Film generation

```bash
./v35_generator \
    --thickness 100 \
    --oil pmps \
    --oil-length 16 \
    --oil-wt 10
```

Film mode uses `boundary p p f`, fixed `Lz`, lateral compression, lateral NPT,
and separate `zlo_wall` and `zhi_wall` repulsive fixes. The requested thickness
should come from the corresponding equilibrated bulk result.

## All command-line inputs

| Option | Meaning |
|---|---|
| `--oil pdms|pmps|copolymer` | Explicit component-3 oil model |
| `--oil-length N` | Repeat units per oil chain |
| `--oil-wt X` | Oil percentage of total formulation mass |
| `--mps-percent X` | MPS percentage of copolymer repeat positions |
| `--mps-wt X` | MPS repeat-unit weight percentage within copolymer |
| `--sequence random|alternating|block` | Copolymer sequence |
| `--oil-seed N` | Oil sequence, conformation, and placement seed |
| `--oil-min-separation X` | Initial oil-to-other minimum distance |
| `--n1 N`, `--m1 M` | Network-strand length and count |
| `--n2 N` | Cross-linker length |
| `--n4 N`, `--m4 M` | Moderator size and count |
| `--functionality F` | Cross-linker functionality, 3–16 |
| `--crosslink-distribution MODE` | `random` or `regular` reactive sites |
| `--crosslink-seed N` | Cross-linker site seed |
| `--mass X` | Common DMS bead mass |
| `--density X` | Initial mass density |
| `--target-density X` | Compression target density |
| `--bond-length X` | Initial DMS bond length |
| `--spacing X` | Formulation placement spacing |
| `--thickness X` | Enable film geometry with fixed `Lz` |
| `--seed N` | Star-moderator seed |
| `--output FILE` | Override data filename and companion case name |
| `--help` | Print built-in command help |

`M2` is always calculated from `M2=2*M1/functionality`.
