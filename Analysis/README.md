# Final-snapshot z-profile analyzer

`z_profile.cpp` calculates the through-thickness composition and density
profile of a four-component V22 or V35 system from one final LAMMPS data file.
It is designed for final equilibrated bulk or film configurations, so a
trajectory containing many frames is not required.

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

From WSL:

```bash
cd "/mnt/c/Users/siteng/Documents/Simulation from Github/Analysis"
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

For the supplied V22 test system:

```bash
./z_profile \
    /home/siteng/makedata/V22/10/data.V22_PDMS_N32_10wt.npt_eq \
    /home/siteng/makedata/V22/10/V22_PDMS_N32_10wt.info
```

The default result is written in the current directory as:

```text
z_profile.V22_PDMS_N32_10wt.dat
```

Choose an explicit output path with:

```bash
./z_profile DATA_FILE INFO_FILE \
    --output my_profile.dat
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

The output is a headerless, whitespace-separated numeric table so it can be
loaded directly with MATLAB `load`. Each row is one sample point and contains
14 columns:

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

The current V22 and V35 beads all have the same mass, so these mass fractions
are also bead-number fractions. Reading masses from the data file keeps the
calculation valid for future fillers whose bead types may have different
masses.

An empty window has zero counts and zero relative density. Its four component
fractions are written as `NaN` because a `0/0` local composition is undefined.

## MATLAB example

```matlab
profile = load('z_profile.V22_PDMS_N32_10wt.dat');

z = profile(:,3);
relativeDensity = profile(:,6);
networkFraction = profile(:,11);
crosslinkerFraction = profile(:,12);
fillerFraction = profile(:,13);
moderatorFraction = profile(:,14);

figure;
plot(z, fillerFraction, 'LineWidth', 1.5);
xlabel('z (\AA)');
ylabel('Local filler mass fraction');
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
