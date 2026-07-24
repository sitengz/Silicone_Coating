# V35 LAMMPS data generator

`v35_generator.cpp` creates a complete four-file model package for the coarse-grained V35 formulation: a LAMMPS `atom_style full` data file, matching LAMMPS input, matching Slurm submission script, and a JSON model manifest with the `.info` extension.

1. `N1`, `M1`: bifunctional network strands. The two terminal beads are reactive atom type 2.
2. `N2`, `M2`: multifunctional cross-linkers. Reactive sites are atom type 3.
3. `N3`, `M3`: neutral linear PDMS oil filler. All beads are atom type 1.
4. `N4`, `M4`: five-bead star-like moderators. The center is neutral type 1 and the four arms are reactive type 2.

Atom type 1 is neutral. Types 2 and 3 may react during later simulations. The source preserves the legacy declaration of four atom types, although atom type 4 is currently unused.

## Build

```bash
g++ -std=c++17 -O2 -Wall -Wextra -pedantic v35_generator.cpp -o v35_generator
```

Display the built-in command reference with:

```bash
./v35_generator --help
```

## Generated model package

Every successful run writes four matching files in the same directory. For the
default formulation:

```text
data.V35_PDMS_N32_10wt
in.V35_PDMS_N32_10wt
submit.V35_PDMS_N32_10wt.sh
V35_PDMS_N32_10wt.info
```

The LAMMPS input reads the exact generated data filename. The Slurm script
changes to the directory recorded by `SLURM_SUBMIT_DIR` before launching
LAMMPS and reads the matching input filename, so the package can be copied to
the cluster and submitted from that directory as one unit.

The `.info` file is valid JSON and records composition, realized filler loading,
box dimensions, topology counts, seeds, cross-linker sites, generated filenames,
compression settings, and the selected bulk or film workflow. Its
`format_version` field is intended to support future analysis scripts.

## Default setting

Run the default formulation with:

```bash
./v35_generator
```

| Component | Meaning | Chain size `N` | Molecules `M` | Bead types |
|---|---|---:|---:|---|
| 1 | Bifunctional network strand | `N1=384` | `M1=306` | Type 1 interior; type 2 ends |
| 2 | Tetrafunctional cross-linker | `N2=32` | `M2=153` | Type 1 neutral; type 3 reactive sites |
| 3 | Linear PDMS filler | `N3=32` | `M3=425` | Type 1 |
| 4 | Star-like moderator | `N4=5` | `M4=6` | Type 1 center; type 2 arms |

The default molecule counts follow directly from stoichiometry and the V22-sized bead target:

```text
M2 = 2*M1/functionality = 2*306/4 = 153
base beads = 384*306 + 32*153 + 5*6 = 122430
10 wt% filler: N3=32, M3=425
total beads = 122430 + 32*425 = 136030
```

Other defaults:

| Setting | Default |
|---|---:|
| Cross-linker functionality | `4` |
| Cross-linker site distribution | `random` |
| Cross-linker random seed | `20260722` |
| PDMS loading | approximately `10 wt%` |
| Bead mass | `74 g/mol` |
| Initial packing density | `0.1 g/cm^3` |
| Density after scripted compression | `0.8 g/cm^3` |
| Initial bond length | `2.801 angstrom` |
| Molecule placement spacing | `7.0 angstrom` |
| Geometry | Cubic bulk box |
| Star-moderator coordinate seed | `5489` |
| Data filename | `data.V35_PDMS_N32_10wt` |

`M2=153` is derived from stoichiometry, not independently specified. The V35 defaults were selected so the fixed base and the complete 10 wt% system contain exactly the same bead counts as V22: 122430 base beads and 136030 total beads.

## All command-line inputs

| Option | Value | Default | Meaning and restrictions |
|---|---|---:|---|
| `--n1` | integer | `384` | Beads per network strand; nonnegative |
| `--m1` | integer | `306` | Number of network strands; determines stoichiometric `M2` |
| `--n2` | integer | `32` | Beads per cross-linker; must be at least the functionality |
| `--n3` | integer | `32` | Manual PDMS chain length |
| `--m3` | integer | `425` | Manual PDMS chain count; overridden by `--filler-wt` |
| `--n4` | integer | `5` | Moderator size; must remain 5 when `M4>0` |
| `--m4` | integer | `6` | Number of star-like moderators |
| `--filler-length` | integer | `32` | Recommended PDMS chain-length input; alias for `--n3` |
| `--filler-wt` | number | `10` conceptually | PDMS percent of total weight; range `0 <= X < 100`; calculates `M3` |
| `--functionality` | integer | `4` | Reactive sites per cross-linker; `3 <= F <= min(16,N2)` |
| `--crosslink-distribution` | text | `random` | Either `random` or `regular` |
| `--crosslink-seed` | integer | `20260722` | Reproducible random type-3 site selection |
| `--mass` | number | `74` | Common bead mass in g/mol; must be positive |
| `--density` | number | `0.1` | Initial packing density in g/cm3; must be positive |
| `--target-density` | number | `0.8` | Density reached by scripted compression; must be positive |
| `--bond-length` | number | `2.801` | Initial bond length in angstrom; must be positive |
| `--spacing` | number | `7.0` | Initial molecule-placement spacing; must be positive |
| `--thickness` | number | omitted | Enables film mode with fixed `Lz` in angstrom |
| `--seed` | integer | `5489` | Reproducible star-moderator coordinate seed |
| `--output` | path | automatic | Overrides the data filename; companion files use the same directory and derived case name |
| `--help` | none | — | Prints built-in help and exits |

`--m2` is deliberately rejected because the program calculates `M2` from stoichiometry.

For normal V35 work, leave components 1, 2, and 4 fixed and vary the PDMS using `--filler-length` and `--filler-wt`. The lower-level `--n3` and `--m3` inputs are retained for manual studies.

## PDMS filler controls

Change the PDMS chain length and weight percentage with:

```bash
./v35_generator --filler-length 64 --filler-wt 15
```

The V35 base consists of components 1, 2, and 4. With equal bead masses, the program calculates:

```text
base_beads   = N1*M1 + N2*M2 + N4*M4
filler_beads = base_beads * (filler_wt/100) / (1 - filler_wt/100)
M3           = round(filler_beads / filler_length)
```

The example writes `data.V35_PDMS_N64_15wt`. Since `M3` must be an integer, the program prints the realized weight percentage, which may differ slightly from the request.

Either input may be used alone:

```bash
./v35_generator --filler-length 64  # retains 10 wt%
./v35_generator --filler-wt 15      # retains N32
./v35_generator --filler-wt 0       # no PDMS filler
```

Use `--output FILE` to replace the automatically generated filename.

## Film geometry

Omitting `--thickness` preserves the original cubic bulk system. Supplying a
positive thickness enables film mode:

```bash
./v35_generator --thickness 100
```

The requested value becomes the fixed `Lz`. The program preserves the target
density by calculating equal lateral dimensions:

```text
volume = total_bead_mass / density
Lx = Ly = sqrt(volume / Lz)
```

The example automatically writes
`data.V35_PDMS_N32_10wt_film_Lz100`. Film geometry can be combined with the
other composition controls:

```bash
./v35_generator \
  --thickness 80 \
  --filler-length 64 \
  --filler-wt 15
```

The film thickness should be taken from the corresponding equilibrated bulk
simulation and supplied through `--thickness`. The manifest records this
bulk-to-film provenance expectation.

The generated film input uses `boundary p p f` and repulsive
`wall/lj126` surfaces at `zlo EDGE` and `zhi EDGE`. It never deforms or
barostats `z`. Compression and final NPT control operate only in `x` and
`y`, with `couple xy`.

The compression scale is calculated from the initial and target densities:

```text
bulk scale per dimension = (initial_density / target_density)^(1/3)
film x/y scale           = (initial_density / target_density)^(1/2)
```

With the defaults, bulk uses `0.5` in x/y/z and film uses approximately
`0.35355339` in x/y while keeping Lz fixed.

## Cross-linker stoichiometry

A network requires cross-linker functionality of at least 3. The maximum supported functionality is 16, and functionality cannot exceed `N2`.

Every network strand has two reactive type-2 ends, so the program enforces:

```text
M2 * functionality = 2 * M1
M2 = 2*M1/functionality
```

Exact stoichiometry requires integer `M2`. If `2*M1` is not divisible by the requested functionality, generation stops with an error rather than silently creating an imbalance.

For the default `M1=306`, the exactly compatible functionalities between 3 and 16 are:

```text
3, 4, 6, 9, 12
```

Examples:

```bash
./v35_generator --functionality 4   # M2 = 153
./v35_generator --functionality 12  # M2 = 51
./v35_generator --functionality 8   # error: 612/8 is not an integer
```

## Reactive-site distribution

Regular distribution places type-3 sites at beads `1,3,5,7,...`:

```bash
./v35_generator \
  --functionality 8 \
  --crosslink-distribution regular
```

Regular placement additionally requires `functionality <= ceil(N2/2)`.

Random distribution selects distinct positions from `1...N2`:

```bash
./v35_generator \
  --functionality 12 \
  --crosslink-distribution random \
  --crosslink-seed 314159
```

The same cross-linker seed reproduces the same positions. `--crosslink-seed` affects only reactive-site selection; `--seed` independently controls star-moderator coordinates. Selected cross-linker sites are printed during generation.

## Complete example

This example retains the fixed V35 network strands and moderators, uses functionality-12 cross-linkers with random sites, and adds N64 PDMS at 15 wt%:

```bash
./v35_generator \
  --functionality 12 \
  --crosslink-distribution random \
  --crosslink-seed 314159 \
  --filler-length 64 \
  --filler-wt 15 \
  --seed 5489
```

The program calculates `M2=51`, calculates the required integer `M3`,
reports the realized composition, and writes the data, LAMMPS input, Slurm
script, and `.info` manifest for `V35_PDMS_N64_15wt`.

## Generated simulation workflows

The generated input preserves the supplied V35 workflows while avoiding
temperature-variable redefinition, which is not portable to the 2023 LAMMPS
module used by the Slurm template. Temperature-dependent LJ values are
calculated by the generator and written numerically.

The bulk workflow contains 7,000,000 total steps:

```text
1M  relaxation at 800 K
1M  isotropic compression with crosslinking at 800 K
1M  relaxation with continued crosslinking
2M  continued crosslinking at fixed dimensions
1M  cooling from 800 K to 300 K under isotropic NPT
1M  equilibration at 300 K under isotropic NPT
```

The film workflow contains 7,000,000 total steps:

```text
1M  relaxation at 800 K
1M  lateral compression with crosslinking at fixed Lz
1M  relaxation with continued crosslinking
2M  continued crosslinking at fixed dimensions
1M  cooling from 800 K to 300 K under lateral NPT
1M  equilibration at 300 K under lateral NPT
```

In both geometries, `fix bond/create` remains active for 4,000,000 steps, beginning before compression so reactive sites can bond as compression brings them together.

The generated Slurm script preserves the provided cluster defaults: one node,
48 MPI tasks, 200 GB memory, the `nova` partition, and the
`lammps/20230802.2-py310-openmpi4-ezoqd7f` module.

## `.info` manifest

Although its extension is `.info`, the manifest uses JSON. For example:

```json
{
  "format": "V35-model-info",
  "format_version": 1,
  "geometry": "film",
  "film_thickness_angstrom": 257.07,
  "composition": {
    "requested_filler_weight_percent": 5.0,
    "realized_filler_weight_percent": 4.99138613
  }
}
```

Future analysis tools should check `format` and `format_version` before
reading fields. The manifest intentionally records both requested and realized
filler weight percentages because integer molecule counts introduce rounding.

## Modeling notes

- Topology counts in the LAMMPS header are derived from the topology actually generated.
- Linear molecules use bond type 1. Star-moderator arms use bond type 2.
- N384 network strands use an in-box folded serpentine initial conformation because their 1072.8-angstrom straight contour is longer than the initial box; every initial bond remains exactly 2.801 angstrom.
- The star moderator preserves the original planar five-bead geometry and its six arm-arm angles.
- Star-moderator centers preserve the original narrow z slab, from `0.02L` to `0.12L`.
- Coordinates are initial placements and can overlap under periodic wrapping. Energy minimization or another packing/overlap-removal procedure is needed before production dynamics.
- Neutral PDMS filler chains start at the z=0 midplane and add layers toward positive z after filling the lateral grid.
- Film wall positions are tied to the generated z box edges rather than hardcoded coordinates.
- The filler calculation currently assumes the same coarse-grained bead mass for all components.

