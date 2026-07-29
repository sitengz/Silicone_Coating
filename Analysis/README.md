# Final-snapshot z-profile analyzer

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

Two results are written in the current directory by default:

```text
z_profile.V22_PDMS_N32_10wt.dat
z_regions.V22_PDMS_N32_10wt.dat
```

Choose explicit output paths with:

```bash
./z_profile DATA_FILE INFO_FILE \
    --output my_profile.dat \
    --region-output my_regions.dat
```

No source or input files are modified.

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
profile = load('z_profile.V22_PDMS_N32_10wt.dat');

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
regions = load('z_regions.V22_PDMS_N32_10wt.dat');
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
