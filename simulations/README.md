# V22/V35 coating series: seven oils × two loadings × two elastomers, plus no-oil controls

The 15 case folders under `V22/` contain V22 configurations, and the
15 folders under `V35/` contain V35 configurations. Fourteen per formulation
use 5 or 10 wt% of each of the seven oils in the separate Silicone_Oil
repository's numbered
`simulations/01`–`07` series. Only the oil **chemistry, chain length, and MPS
repeat fraction** are transferred. The oil-only chain counts (which target
100,000 oil repeats in a pure-oil system) are not used: the V22 generator
calculates its own oil chain count to reach each coating weight percentage.
For copolymers, `mps_distribution = balanced` allocates the rounded **total**
MPS repeat count across chains, as the oil generator does. Individual chains
may differ by one pendant bead; `.info` records the base MPS count, the
number of chains with one extra, and the realized mol%.

| Oil sample | Chemistry | Oil length | MPS mol% | Case suffixes (both bases) |
|---|---|---:|---:|---|
| 01 | PDMS | 30 | 0 | `PDMS_N30_5wt`, `PDMS_N30_10wt` |
| 02 | PMPS | 12 | 100 | `PMPS_N12_5wt`, `PMPS_N12_10wt` |
| 03 | random copolymer | 179 | 5 | `Copol_N179_MPS5mol_5wt`, `Copol_N179_MPS5mol_10wt` |
| 04 | random copolymer | 42 | 10 | `Copol_N42_MPS10mol_5wt`, `Copol_N42_MPS10mol_10wt` |
| 05 | random copolymer | 44 | 10 | `Copol_N44_MPS10mol_5wt`, `Copol_N44_MPS10mol_10wt` |
| 06 | random copolymer | 65 | 10 | `Copol_N65_MPS10mol_5wt`, `Copol_N65_MPS10mol_10wt` |
| 07 | random copolymer | 15 | 50 | `Copol_N15_MPS50mol_5wt`, `Copol_N15_MPS50mol_10wt` |

`MPS5mol` denotes the MPS repeat-unit mole percentage; the final `5wt` or
`10wt` denotes the oil loading in the complete coating formulation.

`NoOil_0wt` is the oil-free elastomer control under each formulation. These
configs retain the normal network strands, crosslinkers, and star moderators;
only component 3 is absent (`N3=0`, `M3=0`). They do not specify `oil`,
`oil_length`, or `oil_wt`, because omission is how the generators select no oil.
The generated cases are `V22_NoOil_0wt` and `V35_NoOil_0wt`.

Each `V22/*/model.conf` fixes the V22 base (`N1=128`, `M1=900`, functionality 8).
The oil-bearing configs additionally set oil identity, chain length, MPS
fraction, and requested oil wt%. Values can be overridden on the generator
command line. The generated `.info` records the actual integer oil chain count
and realized weight percentage. The `V35/*/model.conf` files instead fix
`N1=384`, `M1=306`, and functionality 4; the oil definitions and loadings
are otherwise identical.

From the repository root, compile and generate one bulk case:

```bash
make generators
./build/v22_generator --config simulations/V22/Copol_N179_MPS5mol_5wt/model.conf
./build/v35_generator --config simulations/V35/Copol_N179_MPS5mol_5wt/model.conf
```

For the two oil-free controls, use:

```bash
./build/v22_generator --config simulations/V22/NoOil_0wt/model.conf
./build/v35_generator --config simulations/V35/NoOil_0wt/model.conf
```

To build every C++ generator and analyzer, run `make` at the repository root.
To generate all 30 bulk packages from the committed configs, run
`make -C simulations`. Use `make -C simulations bulk-v22` or `bulk-v35`
for only one elastomer series. Make skips packages whose `.info` is newer
than their config and generator; `make -B -C simulations bulk` forces a
regeneration. Regenerating a case overwrites its generated initial data,
input, submit script, and `.info`, so avoid doing that in a case with
hand-edited inputs or an active simulation.

The generated data, LAMMPS input, Slurm script, and `.info` are in the
sample's child folder, e.g. `simulations/V22/Copol_N179_MPS5mol_5wt/V22_Copol_N179_MPS5mol_5wt/` or
`simulations/V35/Copol_N179_MPS5mol_5wt/V35_Copol_N179_MPS5mol_5wt/`.
Earlier generated packages in the former top-level folders, if present
locally, are left untouched; these configs generate packages under `V22/`
and `V35/`.
Only the configs and generator tooling are intended for GitHub; generated
data and job files can be recreated with the commands above.
Submit only after reviewing the inputs and cluster resources.

After every bulk run has written its final `data.<case>.npt_eq`, generate the
corresponding independently built film with the bulk's measured `Lz`:

```bash
make -C simulations films-dry-run
make -C simulations films
```

These targets inspect both elastomer series. `films-v22` and `films-v35`
are also available separately. The film script skips cases with no completed
bulk NPT data. It generates files but does **not** submit jobs. Each film
package contains the wall-bounded curing input and a separate
`in.<case>.surface` / `submit.<case>.surface.sh` measurement pair. The surface
input expands the box by 50 Å per face by default and places separate 300 K
LJ guard walls at the new z edges. It also
contains `submit.<case>.chain.sh`, a login-node launcher that queues the
surface job with a Slurm `afterok` dependency on the film job. For example:

```bash
bash simulations/V22/PDMS_N30_5wt/V22_PDMS_N30_5wt_film/submit.V22_PDMS_N30_5wt_film.chain.sh --dry-run
bash simulations/V22/PDMS_N30_5wt/V22_PDMS_N30_5wt_film/submit.V22_PDMS_N30_5wt_film.chain.sh
```

The first command checks the files without submitting anything. The second
submits both jobs and prints their IDs; it returns immediately rather than
waiting for either simulation. Run it **once per case** from a Slurm login
environment, after reviewing the generated inputs and resource requests.
Running it again would submit a duplicate pair. If the film job fails, the
surface job does not run. The surface input reads the film's final
`data.<case>.npt_eq` file; do not use the periodic bulk network as a slab.
Regenerating film packages with `make ... films` overwrites generated inputs
and initial data, so do not repeat it for active or hand-edited cases.

For a completed surface job, analyze its production stress file with:

```bash
python3 Analysis/surface_stress.py \
    simulations/V22/Copol_N179_MPS5mol_5wt/V22_Copol_N179_MPS5mol_5wt_film/V22_Copol_N179_MPS5mol_5wt_film.info
```

The analyzer creates a sample-named output folder with raw time series,
block means, and a JSON summary. The measured quantity is *apparent
mechanical surface stress*, not automatically thermodynamic surface tension.
Inspect block drift, film integrity, oil enrichment, and a matched bulk
stress reference before interpreting it.
