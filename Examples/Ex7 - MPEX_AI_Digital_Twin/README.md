# Step 1 Magnetic Field Generation

This note explains how `Step_1_generate_Bfield.m` generates the Proto-MPEX
magnetic field and why complete elliptic integrals appear in the calculation.

## Files Involved

- `Step_1_generate_Bfield.m`: main script for selecting coil currents,
  computing the field, plotting flux surfaces, and writing output files.
- `CoilSetup_ProtoMPEX.xlsx`: coil geometry table. The script reads sheet
  `conf_G`.
- `../../Functions/CreateCoilStructure.m`: converts the spreadsheet geometry
  and selected power-supply currents into current-loop filament locations.
- `../../Functions/CalculateMagField.m`: sums the magnetic field from every
  current-loop filament.
- `../../Functions/bfield_circular_coil_analytic.m`: analytic magnetic field
  and vector potential for one circular current loop, evaluated with complete
  elliptic integrals.

## Overall Workflow

`Step_1_generate_Bfield.m` builds an axisymmetric vacuum magnetic field in the
`(r,z)` plane. The field has radial and axial components, `Br` and `Bz`; the
toroidal component `Bt` is written as zero.

The script follows this sequence:

1. Add the shared `Functions` directory to the MATLAB path.
2. Read the `conf_G` coil geometry from `CoilSetup_ProtoMPEX.xlsx`.
3. Define the field grid:

   ```matlab
   r1D = linspace(1e-6, 0.1, 400);
   z1D = linspace(z_Dump, z_Target, 500);
   ```

   The radial grid starts at `1e-6` m instead of zero because the circular-loop
   analytic formulas contain powers of `1/r`. The physical on-axis field has a
   finite limit, but this implementation avoids evaluating the singular
   algebraic form exactly at `r = 0`.

4. Define the coil currents for the selected operating case.
5. Build the coil filament model with `CreateCoilStructure`.
6. Compute `Br`, `Bz`, azimuthal vector potential `Atheta`, and magnetic flux
   `phi` with `CalculateMagField`.
7. Use contours of normalized flux to draw field-line/plasma-edge geometry.
8. Write the field to CSV and NetCDF outputs.

## Current Selection

The script contains many commented current sets. The active section currently
assigns `coilCurrents{1}` twice:

```matlab
coilCurrents{1}.TR1 = 530;
coilCurrents{1}.TR2 = 2300;
coilCurrents{1}.PS1 = 4500;
coilCurrents{1}.PS2 = 3500;
coilCurrents{1}.PS3 = 430;

coilCurrents{1}.TR1 = 530;
coilCurrents{1}.TR2 = 2300;
coilCurrents{1}.PS1 = 6800;
coilCurrents{1}.PS2 = 4000;
coilCurrents{1}.PS3 = 160;
```

Because both blocks use index `{1}`, the second block overwrites the first.
Therefore the effective current case is:

```matlab
TR1 = 530 A
TR2 = 2300 A
PS1 = 6800 A
PS2 = 4000 A
PS3 = 160 A
```

To compare multiple cases in one run, the cases need distinct indices, such as
`coilCurrents{1}`, `coilCurrents{2}`, etc.

## Coil Discretization

`CreateCoilStructure` turns each finite-size physical coil into an array of
circular current-loop filaments.

For each coil in the spreadsheet, the function reads:

- axial position `z`
- axial width `dz`
- inner and outer radii, `r_inner` and `r_outer`
- radial and axial winding counts, `layers_r` and `layers_z`
- power-supply name `ps`
- datum convention `datum`

It then computes:

```matlab
nfil  = layers_z .* layers_r;
dzfil = dz ./ layers_z;
drfil = (r_outer - r_inner) ./ layers_r;
```

For datum value `5`, the filament positions are placed across the rectangular
coil cross-section:

```matlab
coil{ii}.zfil = zfil_rel + z_ii - dz_ii/2;
coil{ii}.rfil = rfil_rel + r_ii - dr_ii/2;
```

As implemented, each filament is treated as one circular turn carrying the
assigned power-supply current. The finite coil field is then approximated by
superposition of all those circular-loop fields.

## Field Calculation

`CalculateMagField` first creates the requested 2D grid:

```matlab
[r2D, z2D] = meshgrid(r1D, z1D);
```

It initializes the total field:

```matlab
Br2D = 0;
Bz2D = 0;
Atheta2D = 0;
```

Then it loops over every coil and every filament in that coil:

```matlab
for ii = 1:numel(coil)
    current = coil{ii}.current;
    for jj = 1:coil{ii}.nfil
        [Br0, Bz0, Atheta0] = bfield_circular_coil_analytic( ...
            coil{ii}.rfil(jj), coil{ii}.zfil(jj), r2D, z2D);

        Br_n{ii}     = Br_n{ii}     + Br0 * current;
        Bz_n{ii}     = Bz_n{ii}     + Bz0 * current;
        Atheta_n{ii} = Atheta_n{ii} + Atheta0 * current;
    end

    Br2D     = Br2D     + Br_n{ii};
    Bz2D     = Bz2D     + Bz_n{ii};
    Atheta2D = Atheta2D + Atheta_n{ii};
end
```

`bfield_circular_coil_analytic` returns the field per ampere for one circular
loop. `CalculateMagField` multiplies that per-ampere result by the supply
current, then adds all contributions. This is the Biot-Savart superposition
principle applied to an axisymmetric collection of circular loops.

After the field is summed, the script computes magnetic flux from the vector
potential:

```matlab
phi2D = 2*pi*Atheta2D.*r2D;
```

For an axisymmetric field, `A_theta` defines the poloidal flux function. Contours
of `phi` are magnetic flux surfaces in the `(r,z)` plane. `Step_1` normalizes
this flux by the value at the limiting vessel location:

```matlab
xi = phi2D / phi0;
```

The contour `xi = 1` is then used as the plasma-edge or limiting flux surface.

## Role of Complete Elliptic Integrals

The core analytic problem is the magnetic field from a circular current loop at
an arbitrary off-axis point `(r,z)`. Starting from the Biot-Savart law, the
source current is integrated around the loop azimuth. For off-axis points, that
azimuthal integral cannot be reduced to elementary functions. It reduces to
complete elliptic integrals of the first and second kind:

```matlab
[K, E] = ellipke(u2);
```

Here `u2` is the elliptic-integral parameter:

```matlab
u2 = 4*r.*a ./ ((r+a).^2 + (z-d).^2);
u  = sqrt(u2);
```

where:

- `a` is the radius of the filament loop.
- `d` is the axial location of the filament loop.
- `r` and `z` are the evaluation coordinates.

`K` and `E` encode the geometry of the circular loop as seen from the field
point. They replace a numerical integration around each loop, making the field
evaluation faster and more accurate for a large grid.

The helper function uses these elliptic integrals to compute:

- `Brg`: radial magnetic field from one loop, per ampere.
- `Bzg`: axial magnetic field from one loop, per ampere.
- `Athetag`: azimuthal vector potential from one loop, per ampere.

The constant

```matlab
cnst = 1.e-7;
```

is `mu0/(4*pi)` in SI units. Therefore, with current in amperes and distances in
meters, the resulting magnetic field is in tesla.

## Output Files

`Step_1_generate_Bfield.m` writes:

- `z1D.csv`: axial grid.
- `r1D.csv`: radial grid.
- `Br2D.csv`: radial field.
- `Bz2D.csv`: axial field.
- `bfield_protoMPEX.nc`: NetCDF field file used by later steps.

The NetCDF file contains:

```text
r   radial coordinate, [nR]
z   axial coordinate, [nZ]
br  radial magnetic field, [nR x nZ]
bt  toroidal magnetic field, [nR x nZ], currently zero
bz  axial magnetic field, [nR x nZ]
```

Before writing NetCDF, the script checks whether the field arrays are in
`[z x r]` order and permutes them to `[r x z]` order if needed.

## Practical Notes

- The generated field is an axisymmetric vacuum field from coil currents only.
  Plasma currents are not included.
- The field is computed by superposition, so changing a coil current linearly
  changes that coil's contribution.
- The elliptic-integral formulas are singular at the exact location of a current
  filament. The evaluation grid used here is inside the plasma/vacuum region,
  away from the coil windings.
- `Bt` is explicitly set to zero in `Step_1`; downstream code should not expect
  a toroidal magnetic field from this file.
