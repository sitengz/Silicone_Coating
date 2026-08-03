# V22/V35 analysis tools

This directory contains four independent analyzers. `z_profile.cpp` measures
the through-thickness composition, `final_snapshot_analyzer.cpp` measures
realized crosslinking and strand structure, and `msd_analyzer.cpp` measures
translational motion and diffusion from `dump.msd.lammpstrj`.
`phase_separation_analyzer.cpp` measures DMS/MPS composition fluctuations,
domain length scales, and PMPS-rich chain clusters. Every analyzer uses the
matching generator `.info` file to identify the model correctly.

`z_profile.cpp` calculates the through-thickness composition and density
profile of a four-component V22 or V35 system from one final LAMMPS data file.
It emphasizes the neutral-oil filler density and enrichment needed to compare
periodic bulk systems with wall-bounded films. It is designed for final
equilibrated bulk or film configurations, so a trajectory containing many
frames is not required.

The analyzer reads:

1. A LAMMPS `atom_style full` data file, such as
   `data.V22_PDMS_N32_10wt.npt_eq`.
2. The matching JSON `.info` file generated with the model.

The `.info` file supplies the molecule counts for the network strands,
cross-linkers, neutral oil filler, and star-like moderators. The data file
supplies the final box dimensions, bead coordinates, atom types, and bead
masses.

## Why molecule IDs are used

Atom type 1 cannot identify the neutral oil by itself because neutral beads
also occur in the network strands, cross-linkers, and moderators. The
generators create molecules component by component, so the analyzer uses these
molecule-ID ranges:

```text
network:     1 ... M1
crosslinker: M1+1 ... M1+M2
filler:      M1+M2+1 ... M1+M2+M3
moderator:   M1+M2+M3+1 ... M1+M2+M3+M4
```

The analyzer checks the bead count of every range against the `.info` file.
This also detects a mismatched data and info file before a misleading profile
is written.

## Compile

From the repository root:

```bash
cd Analysis
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    z_profile.cpp -o z_profile
```

Display the command reference with:

```bash
./z_profile --help
```

## Run

The required inputs are the final data file and its matching `.info` file:

```bash
./z_profile DATA_FILE INFO_FILE
```

For a V22 test system, replace `DATA_DIRECTORY` with the directory containing
the matching data and `.info` files:

```bash
./z_profile \
    DATA_DIRECTORY/data.V22_PDMS_N32_10wt.npt_eq \
    DATA_DIRECTORY/V22_PDMS_N32_10wt.info
```

Each run creates an output folder named from the sample's `case_name` in the
`.info` file. The two numeric results and a text copy of the run report are
written together:

```text
V22_PDMS_N32_10wt/
├── z_profile.V22_PDMS_N32_10wt.dat
├── z_regions.V22_PDMS_N32_10wt.dat
└── z_report.V22_PDMS_N32_10wt.txt
```

This keeps independently analyzed samples separated when several analyses are
run from the same directory. The report contains the same summary printed to
the terminal, including the realized spacing, window width, densities,
component molecule ranges, and output paths.

Choose explicit output paths with:

```bash
./z_profile DATA_FILE INFO_FILE \
    --output my_profile.dat \
    --region-output my_regions.dat \
    --report-output my_report.txt
```

An explicitly supplied path is used exactly as written. Any output whose path
is not supplied still uses the sample-named folder. No source or input files
are modified.

## Sample spacing and overlapping windows

The number of profile points and the amount of averaging are controlled
independently:

```text
sample spacing = distance between neighboring profile points
window width   = z-range used to count beads at each profile point
```

For example:

```bash
./z_profile DATA_FILE INFO_FILE \
    --sample-spacing 1.0 \
    --window-width 5.0
```

This retains approximately one profile point per angstrom. At each point, the
program counts beads in a window approximately 5 Å wide. Neighboring windows
overlap by approximately 4 Å, producing a smoother profile without reducing
the number of plotted points.

Both values default to `1.0 Å`. Changing only `--window-width` does not change
the sample points:

```bash
./z_profile DATA_FILE INFO_FILE --window-width 1.0
./z_profile DATA_FILE INFO_FILE --window-width 2.0
./z_profile DATA_FILE INFO_FILE --window-width 5.0
```

`--bin-width` is retained as an alias for `--window-width`:

```bash
./z_profile DATA_FILE INFO_FILE --bin-width 5.0
```

The program divides `Lz` into the nearest whole number of sample spacings.
The window width is then rounded to the nearest whole number of those actual
spacings:

```text
number_of_points = max(1, round(Lz / requested_sample_spacing))
actual_sample_spacing = Lz / number_of_points
window_points = max(1, round(requested_window_width /
                             actual_sample_spacing))
actual_window_width = window_points * actual_sample_spacing
```

Requested and actual values are printed after each run.

For a bulk system, a window that crosses `zlo` or `zhi` wraps around the
periodic z boundary. For a film, the window is clipped at the lower or upper
wall. Film density normalization uses that point's actual clipped window
width, preventing the wall-adjacent density from being artificially reduced
by nonexistent volume outside the simulation box.

Overlapping points are correlated samples. This is appropriate for a smooth
profile, but they should not be treated as independent measurements when
estimating statistical uncertainty.

## Numeric output

The profile output is a headerless, whitespace-separated numeric table so it
can be loaded directly with MATLAB `load`. Each row is one sample point and
contains 21 columns. Columns 1–14 retain the original layout:

| Column | Quantity |
|---:|---|
| 1 | One-based sample index |
| 2 | Lower edge of the counting window (Å) |
| 3 | Sample-center z coordinate (Å) |
| 4 | Upper edge of the counting window (Å) |
| 5 | Total bead count in the window |
| 6 | Local total mass density divided by the mean system mass density |
| 7 | Network bead count |
| 8 | Cross-linker bead count |
| 9 | Filler bead count |
| 10 | Moderator bead count |
| 11 | Local network mass fraction |
| 12 | Local cross-linker mass fraction |
| 13 | Local filler mass fraction |
| 14 | Local moderator mass fraction |
| 15 | Local total mass density (g/cm³) |
| 16 | Local filler mass density (g/cm³) |
| 17 | Local non-filler matrix mass density (g/cm³) |
| 18 | Local filler density divided by the whole-system mean filler density |
| 19 | Local filler mass fraction divided by the global filler mass fraction |
| 20 | Distance from the nearest wall (Å); `NaN` for periodic bulk |
| 21 | Region code: `-1` lower surface, `0` interior/bulk, `1` upper surface |

The raw component counts are retained in columns 7–10, allowing alternative
normalization in MATLAB. When windows overlap, a bead can appear in several
rows. Therefore, summing a component-count column over all rows does not give
the number of unique beads in the system.

For a periodic bulk profile, columns 2 or 4 may extend beyond the printed box
bounds when a window wraps around. The corresponding count correctly includes
beads from the opposite side of the periodic box. Film window edges are
clipped to the physical box bounds.

## Normalization

For a sampling window centered at `z`, the component mass density is:

```text
rho_i(z) = mass_i(z) / (Lx * Ly * effective_window_width)
```

The relative total density in column 6 is:

```text
relative_density(z) = rho_total(z) / rho_system
```

where `rho_system` is the mean mass density calculated from all atoms and the
final data-file volume. For film points close to a wall,
`effective_window_width` is the clipped width reported by columns 2 and 4.

The final four columns are normalized by the local total mass density:

```text
fraction_i(z) = rho_i(z) / rho_total(z)
```

For every occupied window:

```text
network_fraction
  + crosslinker_fraction
  + filler_fraction
  + moderator_fraction
  = 1
```

Mass fractions are calculated from the `Masses` section rather than assuming
equal bead masses. This matters for the copolymer fillers, whose MPS bead types
do not have the same masses as DMS beads. Raw bead counts remain available in
columns 7–10.

An empty window has zero counts and zero relative density. Its four component
fractions are written as `NaN` because a `0/0` local composition is undefined.

The two filler enrichment measures answer different questions:

```text
column 18 = local filler mass density / mean filler mass density
column 19 = local filler mass fraction / global filler mass fraction
```

Column 18 includes local total-density changes. Column 19 isolates composition:
a value above one means the local material is filler-enriched relative to the
overall formulation.

## Film surface and interior regions

For film geometry, the lower and upper layers are assigned as operational
surface regions:

```text
lower surface: zlo <= z < zlo + surface_width
interior:      zlo + surface_width <= z < zhi - surface_width
upper surface: zhi - surface_width <= z <= zhi
```

Set the layer thickness explicitly when comparing systems:

```bash
./z_profile FILM_DATA FILM_INFO \
    --surface-width 20 \
    --sample-spacing 1 \
    --window-width 5
```

If omitted, `surface_width` is the smaller of 20 Å and 20% of the film
thickness. This is an operational definition, not a claim that every
interfacial effect ends at exactly that distance. Plot columns 18 or 19 against
column 20 to check where the profile reaches its interior plateau, then rerun
with a revised common width if needed.

The region output is also a headerless numeric table. It contains one row for
periodic bulk. A film contains lower-surface, interior, upper-surface, and
combined-surface rows:

| Column | Quantity |
|---:|---|
| 1 | Region code: `-1`, `0`, `1`, or `2` for both surfaces combined |
| 2 | Region lower z bound (Å); `NaN` for the combined row |
| 3 | Region upper z bound (Å); `NaN` for the combined row |
| 4 | Total region width (Å) |
| 5 | Total unique bead count |
| 6 | Unique filler bead count |
| 7 | Region total mass density (g/cm³) |
| 8 | Region filler mass density (g/cm³) |
| 9 | Region matrix mass density (g/cm³) |
| 10 | Region filler mass fraction |
| 11 | Region filler density divided by whole-system mean filler density |
| 12 | Region filler mass fraction divided by global filler mass fraction |

Unlike overlapping profile windows, each atom belongs to exactly one primary
region. The combined-surface row is the sum of the lower and upper regions.

## MATLAB example

```matlab
profile = load(fullfile('V22_PDMS_N32_10wt', ...
    'z_profile.V22_PDMS_N32_10wt.dat'));

z = profile(:,3);
relativeDensity = profile(:,6);
networkFraction = profile(:,11);
crosslinkerFraction = profile(:,12);
fillerFraction = profile(:,13);
moderatorFraction = profile(:,14);
fillerDensity = profile(:,16);
fillerEnrichment = profile(:,19);
wallDistance = profile(:,20);
regionCode = profile(:,21);

figure;
plot(z, fillerFraction, 'LineWidth', 1.5);
xlabel('z (\AA)');
ylabel('Local filler mass fraction');
```

For a film, plot both surfaces on a common distance-from-wall axis:

```matlab
filmRows = regionCode ~= 0;

figure;
scatter(wallDistance(filmRows), fillerEnrichment(filmRows), 20, ...
        regionCode(filmRows), 'filled');
xlabel('Distance from nearest wall (\AA)');
ylabel('Filler enrichment');
```

Load the non-overlapping region summary with:

```matlab
regions = load(fullfile('V22_PDMS_N32_10wt', ...
    'z_regions.V22_PDMS_N32_10wt.dat'));
```

To verify the local normalization:

```matlab
fractionSum = sum(profile(:,11:14), 2);
max(abs(fractionSum - 1), [], 'omitnan')
```

For a profile that is still noisy after overlapping-window averaging, MATLAB
can apply additional post-processing without losing the original counts:

```matlab
smoothFillerFraction = movmean(fillerFraction, 5, 'omitnan');
```

## Intended workflow

1. Generate and equilibrate a bulk system.
2. Use the final bulk result to define the film thickness.
3. Generate and equilibrate the corresponding film.
4. Run this analyzer on the final `*.npt_eq` data file and matching `.info`.
5. Load the numeric profile into MATLAB for plotting, smoothing, averaging
   between independent samples, or further interfacial analysis.

## Crosslinking, connectivity, and end-to-end analyzer

`final_snapshot_analyzer.cpp` is intended for the `*.npt_eq` snapshot produced
by the current V22 or V35 generator and the matching model-info version-2
`.info` file. It works for bulk and film geometries and for no-oil, PDMS,
PMPS, and copolymer formulations.

Do not pair a data file and `.info` file from different generated cases. The
analyzer uses molecule-ID ranges, initial topology counts, reactive-site
counts, crosslink bond mapping, crosslinker functionality, geometry, and run
length from the `.info` file. Several exact checks deliberately stop the run
instead of reporting misleading conversion data.

Version-1 `.info` files are not accepted. Those files predate the rule that
reserves bond type 2 exclusively for bonds created by `fix bond/create`; their
moderator bonds can therefore make a final crosslink count ambiguous.

### Compile

```bash
cd Analysis
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    final_snapshot_analyzer.cpp -o final_snapshot_analyzer
```

### Run

```bash
./final_snapshot_analyzer FINAL_DATA_FILE MODEL_INFO_FILE
```

For example:

```bash
./final_snapshot_analyzer \
    DATA_DIRECTORY/data.V22_PDMS_N32_10wt.npt_eq \
    DATA_DIRECTORY/V22_PDMS_N32_10wt.info
```

By default, both outputs are placed in a folder named from `case_name`:

```text
V22_PDMS_N32_10wt/
├── final_snapshot_report.V22_PDMS_N32_10wt.txt
└── strand_end_to_end.V22_PDMS_N32_10wt.dat
```

Paths can be overridden independently:

```bash
./final_snapshot_analyzer DATA_FILE INFO_FILE \
    --report-output my_report.txt \
    --strand-output my_strands.dat
```

`--output` remains an alias for `--report-output`.

### Crosslink calculations

The current generators define:

```text
atom type 2 = unreacted network terminal or moderator arm
atom type 3 = unreacted crosslinker site
bond type 2 = newly created formulation crosslink only
```

`fix bond/create` changes each reacted type-2 and type-3 atom to type 1.
Consequently, every valid crosslink consumes exactly one site of each type.
The analyzer verifies all of the following before reporting conversion:

```text
net increase in bond count = number of bond-type-2 bonds
number of bond-type-2 bonds = intermolecular formulation crosslinks
final type-2 atoms + crosslinks = initial type-2 sites
final type-3 atoms + crosslinks = initial type-3 sites
```

It also checks the site counts by component. A type-2 bond must join a
crosslinker to either a network strand or a moderator. The report contains:

- overall bond conversion using the smaller of the type-2 and type-3 site
  totals as the stoichiometric maximum;
- crosslinker conversion using all initial type-3 sites;
- total type-2 site conversion;
- network-terminal conversion and moderator-site conversion separately;
- network--crosslinker and moderator--crosslinker bond counts; and
- molecule-degree distributions for strands, crosslinkers, and moderators.

For the current stoichiometry, `M2 * functionality = 2 * M1`. Moderator arms
add extra type-2 sites, so the crosslinker/type-3 side is normally limiting.
The crosslinker conversion is therefore the value most directly comparable to
the intended 95% conversion target.

### Reactive-network connectivity

The connectivity graph includes all network strands, crosslinkers, and
star-like moderators. Silicone-oil molecules are excluded because they do not
participate in `fix bond/create`. Including moderators is important because a
multifunctional moderator can connect more than one crosslinker and therefore
change network connectivity.

The report gives the number of linked and isolated reactive components, the
network/crosslinker/moderator composition of the largest component, and its
bead and mass fractions relative to all reactive material.

### Network-strand end-to-end distances

The end-to-end vector is defined from the first to the last bead of every
bifunctional network strand. It is not restricted to reacted strands. The
program reconstructs the vector by summing bead-to-bead displacements along
the generator-defined ordered covalent path. Minimum-image displacement is
used in periodic directions; z is not wrapped for a film. This remains valid
when endpoint image flags are absent or a long strand spans a periodic box
more than once.

The report provides statistics for:

- all network strands;
- strands whose two terminal sites both reacted; and
- fully reacted strands belonging to the largest reactive component.

The strand table is headerless and numeric for direct MATLAB loading. Its
columns are:

| Column | Quantity |
|---:|---|
| 1 | Network-strand molecule ID |
| 2 | First terminal atom ID |
| 3 | Last terminal atom ID |
| 4 | Intermolecular crosslink degree |
| 5 | Both terminal sites reacted: `1` yes, `0` no |
| 6 | Strand belongs to largest reactive component: `1` yes, `0` no |
| 7 | End-to-end `dx` (Å) |
| 8 | End-to-end `dy` (Å) |
| 9 | End-to-end `dz` (Å) |
| 10 | End-to-end distance (Å) |
| 11 | Squared end-to-end distance (Å²) |

MATLAB example:

```matlab
ree = load('strand_end_to_end.V22_PDMS_N32_10wt.dat');
R = ree(:,10);
R2 = ree(:,11);
fullyReacted = ree(:,5) == 1;
largestNetwork = ree(:,6) == 1;
```

### Scope

One final snapshot supports realized composition, topology, crosslink
conversion, static connectivity, and static end-to-end statistics. Reaction
kinetics require the LAMMPS log or time-resolved trajectory. MSD, diffusion,
and relaxation calculations require the generated production trajectory.

## MSD and diffusion analyzer

`msd_analyzer.cpp` reads the trajectory generated by the final 1M-step, 300 K
NVT stage. The standard formulation trajectory contains 1,001 frames: one at
timestep zero and one every 1,000 simulation steps through timestep 1,000,000.
With the 5 fs timestep, neighboring frames are 5 ps apart and the complete
trajectory spans 5 ns.

The analyzer derives the expected frame count from the production steps and
dump interval, including the timestep-zero snapshot. Consequently, older
`.info` files that declare 1,000 frames are interpreted as 1,001 when their
timing fields specify 1,000,000 steps and a 1,000-step dump interval. The
report preserves the declared value and records a legacy-metadata warning.

The trajectory must contain:

```text
id mol type x y z ix iy iz
```

Wrapped `x y z` coordinates are reconstructed with the image flags. The
analyzer checks that the NVT box is constant and that its periodic directions
match the `bulk` or `film` geometry in the `.info` file.

### Compile

```bash
cd Analysis
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    msd_analyzer.cpp -o msd_analyzer
```

For an exact time-averaged calculation, OpenMP can distribute lag times among
available CPU threads:

```bash
g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic -fopenmp \
    msd_analyzer.cpp -o msd_analyzer
```

### Recommended oil-diffusion calculation

The default calculation is the recommended observable for silicone-oil
translation:

```bash
./msd_analyzer dump.msd.lammpstrj CASE.info
```

It uses:

```text
selection      = filler molecules only
particle       = mass-weighted molecular center of mass
averaging      = all available time origins
drift removal  = whole-system center-of-mass displacement
```

Whole-system drift is removed from every selected particle. The filler COM is
not removed, because doing so would erase the translational motion being
measured. For PMPS and copolymer oil, molecular COMs are mass weighted using
the type-1, type-4, and type-5 masses in the `.info` file.

### Other selections and definitions

Raw, single-origin whole-system bead MSD:

```bash
./msd_analyzer dump.msd.lammpstrj CASE.info \
    --selection all --particle beads --averaging raw
```

Raw filler-bead MSD:

```bash
./msd_analyzer dump.msd.lammpstrj CASE.info \
    --selection filler --particle beads --averaging raw
```

Time-averaged filler-bead MSD:

```bash
./msd_analyzer dump.msd.lammpstrj CASE.info \
    --selection filler --particle beads --averaging time-averaged
```

The raw MSD uses the first available frame as its single origin:

```text
MSD_raw(t) = (1/N) sum_i |r_i(t) - r_i(t_first)|^2
```

The time-averaged MSD uses every available origin by default:

```text
MSD_TA(lag) = [1 / (N * N_origins)]
              sum_origin sum_i |r_i(origin + lag) - r_i(origin)|^2
```

`--origin-stride N` samples every Nth origin. `--max-lag-frames N` limits the
largest lag. These controls are useful for expensive all-bead calculations.
The exact default filler-chain COM calculation is much smaller: a typical
10 wt% system stores only hundreds of COM trajectories rather than more than
100,000 bead trajectories.

Bead MSD and molecular-COM MSD answer different questions. Bead MSD includes
bond vibration and internal chain relaxation. Molecular-COM MSD isolates oil
translation and should normally be used to estimate the silicone-oil
self-diffusion coefficient.

### Diffusion coefficient

The analyzer performs an ordinary linear regression of MSD against lag time.
For time-averaged MSD, the automatic preliminary fit uses 50-90% of the
available lag-time range. The final 10% is excluded because those long lags
have the fewest independent time origins and therefore the greatest sampling
noise. Raw single-origin MSD retains a 50-100% automatic range. A deliberate
range should still be selected after plotting the complete MSD:

```bash
./msd_analyzer dump.msd.lammpstrj CASE.info \
    --fit-start-ns 1.0 --fit-end-ns 4.0
```

The fitted slopes give:

```text
D_x  = slope(MSD_x)  / 2
D_y  = slope(MSD_y)  / 2
D_z  = slope(MSD_z)  / 2
D_xy = slope(MSD_xy) / 4
D_3D = slope(MSD_3D) / 6
```

Diffusion is reported in Å²/ps and cm²/s, using
`1 Å²/ps = 10^-4 cm²/s`. The table also includes `MSD/(2*d*t)` at each lag,
but the fitted long-time slope is the preferred estimate. A negative slope or
poor `R²` indicates that the selected interval is not a valid diffusive
regime.

For films, inspect `MSD_xy` and `D_xy` separately. Long-time z motion is
bounded by the walls, so a 3D diffusion coefficient can become physically
misleading even though it is still included for comparison.

### Output files

Default outputs are grouped by `case_name`, for example:

```text
V22_Copolymer_N15_MPS8of15_random_10wt/
├── msd_filler_molecule_com_time_averaged.V22_Copolymer_N15_MPS8of15_random_10wt.dat
└── msd_report_filler_molecule_com_time_averaged.V22_Copolymer_N15_MPS8of15_random_10wt.txt
```

Explicit paths can be supplied with `--output` and `--report-output`.

The numeric MSD table is headerless for MATLAB. The zero-lag point is not
written because its MSD is always zero and the
instantaneous diffusion ratios divide by zero. The first exported row is
therefore lag index 1 with a positive lag time. Its columns are:

| Column | Quantity |
|---:|---|
| 1 | Zero-based lag index |
| 2 | Lag in LAMMPS steps |
| 3 | Lag time (ps) |
| 4 | Lag time (ns) |
| 5 | Number of time origins |
| 6 | Number of particles/entities |
| 7 | `MSD_x` (Å²) |
| 8 | `MSD_y` (Å²) |
| 9 | `MSD_z` (Å²) |
| 10 | `MSD_xy = MSD_x + MSD_y` (Å²) |
| 11 | `MSD_3D` (Å²) |
| 12 | Instantaneous `MSD_3D/(6t)` (Å²/ps) |
| 13 | Instantaneous `MSD_3D/(6t)` (cm²/s) |
| 14 | Instantaneous `MSD_xy/(4t)` (Å²/ps) |
| 15 | Instantaneous `MSD_z/(2t)` (Å²/ps) |

MATLAB example:

```matlab
msd = load('msd_filler_molecule_com_time_averaged.CASE.dat');
time_ns = msd(:,4);
msd3D = msd(:,11);
msdXY = msd(:,10);

figure;
plot(time_ns, msd3D, 'LineWidth', 1.5);
xlabel('Lag time (ns)');
ylabel('MSD (\AA^2)');
```

A short trajectory excerpt is suitable for testing the parser and formulas,
but not for reporting diffusion. Production diffusion analysis should use the
complete trajectory and a visually confirmed linear fit interval.

## DMS/MPS phase-separation analyzer

`phase_separation_analyzer.cpp` quantifies PMPS-rich segregation in V22 and
V35 formulation trajectories. It combines three measurements that answer
different questions:

1. a local-composition field and excess segregation index;
2. the normalized DMS/MPS concentration structure factor; and
3. intermolecular contact clusters of PMPS-containing oil chains.

The analyzer automatically accepts either the final `data.*.npt_eq` snapshot
or the existing `dump.msd.lammpstrj`, together with the matching model-info
version-2 `.info` file. A final snapshot is sufficient for static morphology
and is the recommended input when a PMPS-rich system is nearly frozen. The
trajectory mode remains available when time-dependent coarsening is useful.
It supports PMPS and DMS/MPS copolymer oils in both bulk and film systems.
PDMS-only and no-oil cases contain no MPS markers and are intentionally
rejected.

### Repeat-marker definition

Raw bead counts would bias the composition because one DMS repeat is one bead
whereas one MPS repeat contains a backbone and a pendant bead. The analyzer
therefore uses exactly one marker per repeat:

```text
DMS marker = atom type 1, 2, or 3
MPS marker = atom type 4 (MPS backbone)
type 5     = validated MPS pendant, not counted as a second repeat
```

The type-4 and type-5 counts must both equal
`oil_chain_count * mps_repeats_per_chain` from the `.info` file. All type-4
and type-5 atoms must lie inside the component-3 molecule-ID range. These
checks prevent a mismatched trajectory and `.info` file from producing a
plausible but incorrect phase metric.

### Compile

```bash
cd Analysis
g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic \
    phase_separation_analyzer.cpp -o phase_separation_analyzer
```

### Run

Recommended static analysis:

```bash
./phase_separation_analyzer data.CASE.npt_eq CASE.info
```

Optional trajectory analysis:

```bash
./phase_separation_analyzer dump.msd.lammpstrj CASE.info
```

The default settings are:

```text
frame stride        = 10
target grid spacing = 8 A
MPS contact cutoff  = 8 A
q maximum           = 80% of the grid Nyquist limit
```

In trajectory mode, the final frame is always analyzed even if it is not on
the requested stride. `--frame-stride` has no effect on a static snapshot.
Controls and explicit output paths are available:

```bash
./phase_separation_analyzer DATA_OR_TRAJECTORY CASE.info \
    --frame-stride 10 \
    --grid-spacing 8.0 \
    --contact-cutoff 8.0 \
    --q-max 0.30 \
    --output phase_metrics.dat \
    --structure-output phase_structure_factor.dat \
    --field-output phase_field_final.dat \
    --report-output phase_report.txt
```

The grid dimensions are the nearest powers of two for the internal FFT.
Consequently, the realized grid spacing can be slightly larger or smaller
than the requested target and is recorded in the report. Use the same
requested grid spacing when comparing formulations.

### Local segregation index

For each voxel, the local MPS repeat fraction is

```text
phi_i = N_MPS,i / (N_MPS,i + N_DMS,i)
```

The analyzer calculates the marker-count-weighted variance of `phi_i`. A
finite system has nonzero variance even when species labels are randomly
mixed, so the exact fixed-composition random-label variance is subtracted.
The reported index is

```text
I = (variance_observed - variance_random)
    / (phi_global*(1-phi_global) - variance_random)
```

`I` near zero is consistent with random mixing at the selected voxel scale;
positive `I` indicates excess local segregation; and `I` approaching one
indicates nearly pure local cells. Statistical fluctuations can make `I`
slightly negative. The analyzer reports both a full 3D index and an XY index
after integrating through z.

For films, the XY index is the preferred phase-separation measure because it
does not mistake ordinary wall-normal layering for lateral PMPS-rich domains.
The existing `z_profile` analyzer should be used separately for surface
segregation. For bulk systems, the 3D index is primary.

### Concentration structure factor

The concentration-fluctuation construction is based on the Bhatia–Thornton
separation of number-density and concentration modes in a binary mixture
[[1]](https://doi.org/10.1103/PhysRevB.2.3004). Partial structure factors and
their interpretation in simulated binary polymer blends are discussed by Cui
et al. [[2]](https://doi.org/10.1021/ma961020h). The analyzer uses the
following finite-system normalization so that formulations with different
numbers of DMS and MPS markers can be compared directly:

The normalized concentration mode is

```text
Scc_normalized(q) = |x_DMS*rho_MPS(q) - x_MPS*rho_DMS(q)|^2
                    / (N*x_MPS*x_DMS)
```

An ideal random-label mixture has a baseline near one. Increasing low-q
intensity indicates large-scale composition fluctuations. The strongest
resolved peak gives the preliminary length scale

```text
domain spacing = 2*pi/q_peak
```

This conversion from the primary structure-factor peak, `q_peak`, to a
real-space spacing follows the usual scattering relation `d = 2*pi/q*`; for
an example combining that interpretation with coarse-grained molecular
simulation, see Srivastava et al. [[3]](https://doi.org/10.1038/ncomms14131).

Bulk systems use a radially averaged 3D `Scc(q)`. Films use in-plane
`Scc(q_xy)` after integrating through z, which isolates lateral domains from
the average z profile. Species are deposited on the composition grid before
the FFT, so the high-q result is grid-smoothed; this analyzer is designed for
the low-q domain regime.

If the peak occurs in the first q shell, the apparent domain size is limited
by the simulation box. In that case `2*pi/q_peak` is a lower bound rather than
a fully resolved domain diameter.

### PMPS-chain contact clusters

This contact graph is an operational definition specific to this analyzer,
not a reproduction of a cluster algorithm from the cited papers. Srivastava
et al. provide a related example in which coarse-grained molecular models and
scattering were used together to distinguish discrete aggregates from
interconnected networks [[3]](https://doi.org/10.1038/ncomms14131).

Two component-3 oil chains are connected when type-4 MPS markers from the two
different molecules lie within `--contact-cutoff`, using minimum-image
distances in periodic directions. Intramolecular contacts are excluded. The
output includes:

- the number of PMPS-containing chain clusters;
- the number and fraction of chains in the largest cluster; and
- the number of unique contacting chain pairs.

The cutoff should eventually be selected from the first minimum of the
intermolecular MPS--MPS radial distribution function. Until that reference is
available, sensitivity should be checked with more than one reasonable
cutoff rather than interpreting the default as a force-field constant.

### Output files

Default outputs are grouped by `case_name`:

```text
CASE/
├── phase_metrics.CASE.dat
├── phase_structure_factor.CASE.dat
├── phase_field_final.CASE.dat
└── phase_report.CASE.txt
```

All `.dat` files are headerless, whitespace-separated numeric tables for
MATLAB.

`phase_metrics` columns:

| Column | Quantity |
|---:|---|
| 1 | Zero-based trajectory frame index |
| 2 | LAMMPS timestep |
| 3 | Time from first trajectory frame (ns) |
| 4 | MPS repeat-marker count |
| 5 | DMS repeat-marker count |
| 6 | Global MPS repeat fraction |
| 7 | Excess 3D segregation index |
| 8 | Observed weighted 3D composition variance |
| 9 | Random-label weighted 3D variance |
| 10 | Excess XY segregation index |
| 11 | Mean normalized `Scc` over the three lowest q shells |
| 12 | q at the strongest resolved `Scc` peak (1/Å) |
| 13 | `2*pi/q_peak` (Å) |
| 14 | Normalized `Scc` at the peak |
| 15 | Peak is in the lowest q shell: `1` yes, `0` no |
| 16 | Number of PMPS-containing chain clusters |
| 17 | Chains in the largest cluster |
| 18 | Fraction of oil chains in the largest cluster |
| 19 | Unique intermolecular PMPS chain contacts |

For a final data snapshot, columns 1–3 are zero because a LAMMPS data file
does not preserve a trajectory frame index or production timestep. The
remaining columns contain the static morphology measurements normally used
for a frozen system.

`phase_structure_factor` columns:

| Column | Quantity |
|---:|---|
| 1 | Zero-based trajectory frame index |
| 2 | LAMMPS timestep |
| 3 | Time from first trajectory frame (ns) |
| 4 | Mean q of the shell (1/Å) |
| 5 | Number of reciprocal-space modes in the shell |
| 6 | Mean normalized `Scc` in the shell |

`phase_field_final` columns:

| Column | Quantity |
|---:|---|
| 1–3 | Integer voxel indices `ix iy iz` |
| 4–6 | Voxel-center coordinates `x y z` (Å) |
| 7 | MPS repeat-marker count |
| 8 | DMS repeat-marker count |
| 9 | Total repeat-marker count |
| 10 | Local MPS repeat fraction |

Empty voxels have zero counts and a reported local fraction of zero. Filter
column 9 greater than zero before using column 10 in a composition histogram.

MATLAB examples:

```matlab
metrics = load('phase_metrics.CASE.dat');
time_ns = metrics(:,3);
seg3D = metrics(:,7);
segXY = metrics(:,10);
lowQ = metrics(:,11);
largestCluster = metrics(:,18);

figure;
tiledlayout(3,1);
nexttile; plot(time_ns, seg3D, time_ns, segXY, 'LineWidth', 1.3);
legend('3D', 'XY'); ylabel('Segregation index');
nexttile; plot(time_ns, lowQ, 'LineWidth', 1.3); ylabel('Low-q S_{cc}');
nexttile; plot(time_ns, largestCluster, 'LineWidth', 1.3);
ylabel('Largest cluster fraction'); xlabel('Time (ns)');
```

Plot the final structure factor with:

```matlab
s = load('phase_structure_factor.CASE.dat');
lastFrame = max(s(:,1));
last = s(:,1) == lastFrame;
plot(s(last,4), s(last,6), 'o-');
xlabel('q (1/\AA)'); ylabel('Normalized S_{cc}(q)');
```

Phase separation should not be assigned from a single number. For a final
snapshot, stronger evidence is the consistent combination of excess local
segregation, enhanced low-q `Scc`, and a large PMPS-rich chain cluster. These
measure morphology but cannot by themselves prove equilibrium. When a useful
trajectory exists, stable time series support a converged morphology; rising
low-q intensity or cluster fraction indicates continued coarsening.

### References for phase-separation analysis

1. A. B. Bhatia and D. E. Thornton, "Structural Aspects of the Electrical
   Resistivity of Binary Alloys," *Physical Review B* **2**, 3004–3012
   (1970). [doi:10.1103/PhysRevB.2.3004](https://doi.org/10.1103/PhysRevB.2.3004)
2. S. T. Cui, H. D. Cochran, P. T. Cummings, and S. K. Kumar, "Computer
   Simulations of the Static Scattering from Model Polymer Blends,"
   *Macromolecules* **30** (11), 3375–3382 (1997).
   [doi:10.1021/ma961020h](https://doi.org/10.1021/ma961020h)
3. S. Srivastava, M. Andreev, A. E. Levi, D. J. Goldfeld, J. Mao, W. T.
   Heller, V. M. Prabhu, J. J. de Pablo, and M. V. Tirrell, "Gel phase
   formation in dilute triblock copolyelectrolyte complexes," *Nature
   Communications* **8**, 14131 (2017).
   [doi:10.1038/ncomms14131](https://doi.org/10.1038/ncomms14131)

References 1–3 support the concentration-mode, polymer-blend scattering, and
peak-spacing interpretations above. The finite-particle-corrected voxel
segregation index and the PMPS chain-contact graph are reproducible analysis
choices introduced for this tool; their definitions and adjustable scales
are documented explicitly rather than attributed to those publications.
