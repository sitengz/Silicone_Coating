# Standalone Silicone Oil Generator

This directory contains the first standalone generator for neutral PDMS oil,
neutral PMPS oil, and rough PDMS-PMPS copolymer oils. It is intentionally
separate from the V22 and V35 formulation generators. The oil model can be
tested independently before it is merged into those systems.

The generator writes four companion files for every case:

```text
data.<case>
in.<case>
submit.<case>.sh
<case>.info
```

The `.info` file is valid JSON and records the composition, realized MPS
monomer and weight percentages, type populations, topology counts, box size,
random seeds, mixing-rule choice, and Green-Kubo production settings. During
the LAMMPS run, a fifth file named `gk_stress.<case>.dat` is created for
viscosity analysis.

## Compile

From the repository root:

```bash
cd Oil
g++ -std=c++14 -O2 -Wall -Wextra -Wpedantic \
    oil_generator.cpp -o oil_generator
```

The code is compatible with C++14 and newer compilers. C++17 can also be used.

## Default model

Running without options:

```bash
./oil_generator
```

generates 625 PMPS oil chains with 16 MPS repeat units per chain:

```text
data.Oil_PMPS_N16_M625
in.Oil_PMPS_N16_M625
submit.Oil_PMPS_N16_M625.sh
Oil_PMPS_N16_M625.info
```

The default system contains 20,000 beads and reproduces the counts of the
supplied PMPS N16 model:

| Quantity | Count |
|---|---:|
| Atoms | 20,000 |
| Bonds | 19,375 |
| Angles | 27,500 |
| Dihedrals | 35,000 |

Unlike the old PMPS generator, the initial chains are non-collinear, molecules
are randomly rotated, every complete chain is placed wholly inside the primary
box with zero image flags, and intermolecular placements below 4.5 Å are
rejected. The default initial density is deliberately low at 0.1 g/cm³.

## Atom types

The combined type namespace is reserved now so that the oil can later be added
to V22 or V35 without renumbering:

| Atom type | Meaning | Mass (g/mol) | Used by oil-only generator |
|---:|---|---:|---|
| 1 | Neutral DMS bead | 74.0 | Yes |
| 2 | Reactive DMS type from formulation model | 74.0 | Reserved |
| 3 | Reactive DMS type from formulation model | 74.0 | Reserved |
| 4 | MPS backbone bead | 59.1204 | Yes |
| 5 | MPS phenyl pendant bead | 77.106 | Yes |

One DMS repeat contains one type-1 bead and has a mass of 74.0 g/mol.

One MPS repeat contains one type-4 backbone bead plus one type-5 pendant bead:

```text
M(MPS repeat) = 59.1204 + 77.106 = 136.2264 g/mol
```

Types 2 and 3 have zero atoms in an oil-only system. They are retained in the
data-file header and pair matrix for compatibility with the formulation models.

## Temperature-dependent nonbonded interactions

All pair coefficients are calculated by the generator and written explicitly.
LAMMPS does not apply an automatic mixing rule.

The DMS-MPS combining rule is:

```text
epsilon(DMS,MPS) = 0.579966
                   * sqrt[epsilon(DMS,DMS) * epsilon(MPS,MPS)]
sigma(DMS,MPS)   = [sigma(DMS,DMS)   + sigma(MPS,MPS)]   / 2
```

The epsilon prefactor is applied to DMS interactions with both the MPS
backbone and MPS pendant beads, at both 300 K and 800 K.

The pure PMPS backbone-pendant term remains the supplied PMPS cross
interaction; it is not recalculated by the DMS-MPS combining rule.

### 300 K attractive matrix

At 300 K, the input uses:

```lammps
pair_style lj/gromacs 12 15
special_bonds lj 0 0 0.5
```

| Pair | Epsilon (kcal/mol) | Sigma (Å) | Source |
|---|---:|---:|---|
| DMS-DMS | 1.012878450 | 6.445843660 | V22/V35 DMS model |
| DMS-MPS backbone | 0.778707979 | 5.960508330 | 0.579966 × geometric ε; arithmetic σ |
| DMS-MPS pendant | 0.537772752 | 6.206674705 | 0.579966 × geometric ε; arithmetic σ |
| MPS backbone-backbone | 1.779864125 | 5.475173000 | PMPS model |
| MPS backbone-pendant | 1.331357250 | 5.727290650 | Supplied PMPS cross term |
| MPS pendant-pendant | 0.848858275 | 5.967505750 | PMPS model |

Types 1, 2, and 3 receive the same DMS nonbonded parameters. Therefore all
three types have the same combined interactions with types 4 and 5.

### 800 K repulsive matrix

The DMS parameters at 800 K are calculated with the same temperature equations
used by the V22/V35 generators. Because only 300 K PMPS parameters were
supplied, this standalone model applies the corresponding DMS temperature
ratios to each PMPS epsilon and sigma:

```text
epsilon_PMPS(800) = epsilon_PMPS(300)
                    * epsilon_DMS(800)/epsilon_DMS(300)

sigma_PMPS(800)   = sigma_PMPS(300)
                    * sigma_DMS(800)/sigma_DMS(300)
```

After scaling the pure PMPS terms, the mixed DMS-MPS values are recalculated
using `0.579966 ×` geometric epsilon and arithmetic sigma.

| Pair | Epsilon at 800 K | Sigma at 800 K (Å) | Repulsive cutoff (Å) |
|---|---:|---:|---:|
| DMS-DMS | 0.531850637 | 6.640519401 | 7.453731009 |
| DMS-MPS backbone | 0.408890459 | 6.140526096 | 6.892507499 |
| DMS-MPS pendant | 0.282378187 | 6.394127126 | 7.177165031 |
| MPS backbone-backbone | 0.934585852 | 5.640532791 | 6.331283990 |
| MPS backbone-pendant | 0.699080133 | 5.900264835 | 6.622823352 |
| MPS pendant-pendant | 0.445725560 | 6.147734850 | 6.900599052 |

The high-temperature stage uses `pair_style lj/cut`. Every pair receives its
own cutoff at `2^(1/6)*sigma`, so only the repulsive-force side of that pair's
Lennard-Jones interaction is retained. The global cutoff is the largest
pair-specific value, 7.453731009 Å.

Both complete 800 K and 300 K matrices are also recorded in each `.info` file.

## Rough bonded model

The copolymer rule is implemented by inspecting every bonded interaction.

### Bonds

| Bond type | Interaction | Potential |
|---:|---|---|
| 1 | DMS-DMS or mixed DMS-MPS backbone bond | DMS harmonic bond |
| 2 | MPS-MPS backbone bond | PMPS backbone harmonic bond |
| 3 | MPS backbone-pendant bond | PMPS pendant harmonic bond |

### Angles

| Angle type | Interaction | Potential |
|---:|---|---|
| 1 | Backbone-only angle containing at least one DMS repeat | DMS harmonic angle |
| 2 | Three MPS backbone beads | PMPS backbone quartic angle |
| 3 | Any angle containing an MPS pendant | PMPS pendant quartic angle |

### Dihedrals

| Dihedral type | Interaction | Potential |
|---:|---|---|
| 1 | Backbone-only dihedral containing at least one DMS repeat | DMS `nharmonic` |
| 2 | Four MPS backbone beads | PMPS backbone `nharmonic` |
| 3 | Dihedral containing one MPS pendant | PMPS one-pendant `nharmonic` |
| 4 | Dihedral containing two MPS pendants | PMPS two-pendant `nharmonic` |

Thus a mixed bonded interaction without a pendant follows the DMS model. Once
a pendant participates, it follows the corresponding PMPS rule.

The input uses:

```lammps
bond_style harmonic
angle_style hybrid harmonic quartic
dihedral_style nharmonic
```

## Composition controls

Every chain has the same number of DMS and MPS repeat units. For `random`
sequences, the positions vary between chains but the composition does not.

### Monomer percentage

Use `--mps-percent` to specify the percentage of repeat-unit positions that are
MPS:

```bash
./oil_generator \
    --length 32 \
    --chains 500 \
    --mps-percent 25 \
    --sequence random
```

This requests 8 MPS and 24 DMS repeat units per chain.

### Weight percentage

Use `--mps-wt` to choose the closest integer MPS count per chain based on the
different DMS and MPS repeat masses:

```bash
./oil_generator \
    --length 32 \
    --chains 500 \
    --mps-wt 40 \
    --sequence random
```

For a chain with `N` repeat positions and `k` MPS repeats:

```text
M_chain = (N-k)*74.0 + k*136.2264

MPS wt% = 100 * k*136.2264 / M_chain
```

Because `k` must be an integer, the realized weight percentage can differ from
the request. The generator reports both values and records them in `.info`.

`--mps-percent` and `--mps-wt` cannot be supplied together.

## Sequence modes

| Mode | Behavior |
|---|---|
| `random` | Selects the exact requested number of MPS sites randomly and reproducibly on each chain |
| `alternating` | Distributes the MPS sites as evenly as possible along each chain |
| `block` | Places one contiguous MPS block in the middle of each chain |

Pure PDMS and pure PMPS chains are unaffected by the sequence choice.

Examples:

```bash
# Pure PDMS oil
./oil_generator --length 32 --chains 500 --mps-percent 0

# Pure PMPS oil
./oil_generator --length 16 --chains 625 --mps-percent 100

# Approximately alternating 50:50 copolymer
./oil_generator \
    --length 32 \
    --chains 500 \
    --mps-percent 50 \
    --sequence alternating

# Central PMPS block at a target MPS weight percentage
./oil_generator \
    --length 64 \
    --chains 250 \
    --mps-wt 30 \
    --sequence block
```

## All command-line inputs

| Option | Type | Default | Meaning |
|---|---|---:|---|
| `--length N` | positive integer | 16 | Repeat-unit positions per chain |
| `--chains M` | positive integer | 625 | Number of oil chains |
| `--n N` | positive integer | 16 | Alias for `--length` |
| `--m M` | positive integer | 625 | Alias for `--chains` |
| `--mps-percent X` | 0–100 | 100 | MPS monomer percentage |
| `--mps-wt X` | 0–100 | unset | Target MPS repeat-unit weight percentage |
| `--sequence MODE` | text | `random` | `random`, `alternating`, or `block` |
| `--density X` | positive number | 0.1 | Initial mass density in g/cm³ |
| `--target-density X` | positive number | 0.8 | Density after scripted 800 K compression |
| `--min-separation X` | positive number below 15 | 4.5 | Minimum intermolecular bead distance in Å |
| `--seed N` | positive integer | 20260727 | Sequence, conformation, rotation, and packing seed |
| `--velocity-seed N` | positive integer | 492845 | LAMMPS initial-velocity seed |
| `--output FILE` | path | automatic | Override the data filename; companion files use the same directory |
| `--help` | — | — | Print command help |

## Generated LAMMPS workflow

The companion input first follows the V22/V35 seven-million-step
equilibration structure without the crosslinking operation, then adds a
100 ns Green-Kubo production stage:

```text
1M  initial relaxation at 800 K
1M  isotropic compression at 800 K
1M  relaxation at the compressed dimensions
2M  extended equilibration at 800 K
1M  cooling from 800 K to 300 K under isotropic NPT
1M  final 300 K equilibration under isotropic NPT
20M  100 ns viscosity production at 300 K under NVT
```

Before dynamics, the input minimizes the structure and creates 800 K
velocities. With the default initial and target densities, each box dimension
is scaled by:

```text
(0.1/0.8)^(1/3) = 0.5
```

The high-temperature repulsive matrix is active for the first five million
steps. The input then switches to the explicit 300 K `lj/gromacs` matrix before
cooling. The timestep is 5 fs, matching the V22/V35 workflow. The complete
workflow is 27 million steps: 7 million equilibration steps followed by
20 million production steps.

### Green-Kubo stress output

After the final 300 K NPT equilibration, the barostat is removed and the final
equilibrated volume is held fixed for a 100 ns NVT trajectory at 300 K. The
production timestep counter and accumulated time are reset to zero.

The atom-coordinate dump is stopped before this long production stage. Instead,
LAMMPS writes the instantaneous off-diagonal pressure components every 10
timesteps, corresponding to a 50 fs sampling interval:

```text
gk_stress.<case>.dat
```

The file has one header followed by exactly four columns:

```text
# time_fs pxy_atm pxz_atm pyz_atm
```

Thus a completed 100 ns run contains 2,000,000 stress samples. The pressure
components are in atmospheres because the input uses LAMMPS `real` units. No
time averaging is applied before writing, preserving the instantaneous stress
series needed for Green-Kubo autocorrelation analysis.

The generated Slurm script still requests 48 hours. Confirm from the shorter
equilibration tests that 27 million steps fit within that wall time; otherwise,
increase the Slurm time request within the limits of the selected Nova
partition.

The Slurm script uses the same Nova module configuration and quoted output
redirection as the tested V22/V35 scripts:

```bash
sbatch submit.<case>.sh
```

## Scope of this version

This version generates periodic bulk oil only. It does not yet:

- insert the oil into V22 or V35;
- generate film geometry or walls;
- create reactive oil end groups;
- infer any mixing rule other than the explicitly documented
  `0.579966 ×` geometric-epsilon, arithmetic-sigma DMS-MPS rule.

Those integrations should be done only after the standalone oil structures and
force-field assignment have been tested.

## Generated file descriptions

- `data.<case>` is the initial LAMMPS data file containing the periodic box,
  oil atoms, and bonded topology.
- `in.<case>` is the complete LAMMPS workflow for 800 K relaxation and
  compression, cooling, 300 K equilibration, and the 100 ns Green-Kubo NVT
  production run.
- `submit.<case>.sh` is the one-node, 96-task Slurm submission script.
- `<case>.info` is a JSON manifest containing composition, sequence,
  force-field, topology, random-seed, and production settings.
- `gk_stress.<case>.dat` is created during the LAMMPS run and contains
  `time_fs`, `pxy`, `pxz`, and `pyz` for viscosity analysis.
