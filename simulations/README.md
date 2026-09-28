# V22/V35 coating series: seven oils × two loadings × two elastomers

The 14 top-level numbered case folders contain V22 configurations, and the
14 numbered folders under `V35/` contain V35 configurations, for 5 and 10 wt% of
each of the seven oils in the separate Silicone_Oil repository's numbered
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
| 01 | PDMS | 30 | 0 | `01_5wt`, `01_10wt` |
| 02 | PMPS | 12 | 100 | `02_5wt`, `02_10wt` |
| 03 | random copolymer | 179 | 5 | `03_5wt`, `03_10wt` |
| 04 | random copolymer | 42 | 10 | `04_5wt`, `04_10wt` |
| 05 | random copolymer | 44 | 10 | `05_5wt`, `05_10wt` |
| 06 | random copolymer | 65 | 10 | `06_5wt`, `06_10wt` |
| 07 | random copolymer | 15 | 50 | `07_5wt`, `07_10wt` |

Each top-level `model.conf` fixes the V22 base (`N1=128`, `M1=900`, functionality 8),
the oil identity, chain length, MPS fraction, and requested oil wt%. These
values can be overridden on the generator command line. The generated `.info`
file records the actual integer chain count and realized weight percentage.
The `V35/*/model.conf` files instead fix `N1=384`, `M1=306`, and
functionality 4; all seven oil definitions and loadings are identical.

From the repository root, compile and generate one bulk case:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic V22/v22_generator.cpp -o V22/v22_generator
./V22/v22_generator --config simulations/03_5wt/model.conf
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic V35/v35_generator.cpp -o V35/v35_generator
./V35/v35_generator --config simulations/V35/03_5wt/model.conf
```

The generated data, LAMMPS input, Slurm script, and `.info` are in the
sample's child folder, e.g. `simulations/03_5wt/V22_Oil03_5wt/` or
`simulations/V35/03_5wt/V35_Oil03_5wt/`. All 28 bulk
packages have been generated locally as initial configurations and **have not been run**.
Only the configs and generator tooling are intended for GitHub; generated
data and job files can be recreated with the commands above.
Submit only after reviewing the inputs and cluster resources.

After every bulk run has written its final `data.<case>.npt_eq`, generate the
corresponding independently built film with the bulk's measured `Lz`:

```bash
python3 simulations/generate_films.py --generator V22/v22_generator --dry-run
python3 simulations/generate_films.py --generator V22/v22_generator
python3 simulations/generate_films.py --formulation V35 --generator V35/v35_generator --dry-run
python3 simulations/generate_films.py --formulation V35 --generator V35/v35_generator
```

The script selects V22 by default; `--formulation V35` selects the V35
series and checks that the supplied generator matches. It skips cases with
no completed bulk NPT data. It does not submit
jobs. The film case gets a regular wall-bounded curing input and a separate
`in.<case>.surface` / `submit.<case>.surface.sh` measurement pair. Run the
film curing job first; then run the surface job, which requires its final
`data.<case>.npt_eq` file. Do not use the periodic bulk network as a slab.

For a completed surface job, analyze its production stress file with:

```bash
python3 Analysis/surface_stress.py \
    simulations/03_5wt/V22_Oil03_5wt_film/V22_Oil03_5wt_film.info
```

The analyzer creates a sample-named output folder with raw time series,
block means, and a JSON summary. The measured quantity is *apparent
mechanical surface stress*, not automatically thermodynamic surface tension.
Inspect block drift, film integrity, oil enrichment, and a matched bulk
stress reference before interpreting it.
