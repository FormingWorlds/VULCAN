# Temperature-dependent UV cross sections

Most laboratory UV cross sections are measured at or below room temperature, which raises
reliability concerns for hot atmospheres. VULCAN can interpolate photoabsorption cross
sections in temperature, layer by layer (Section 2.5 of Tsai et al. 2021 [^tsai2021]).

## Motivation

Heays et al. (2017) [^heays2017] noted that for many molecules a temperature increase of a few
hundred kelvin only broadens the cross section slightly and leaves its wavelength integral
unchanged. However, for molecules with prominent transitions between excited vibrational
states (CO$_2$ in particular) both the absorption threshold and the cross section depend
strongly on temperature, which in turn affects the photolysis rate and the shielding of other
species.

## Available data

Temperature-dependent photoabsorption cross sections are included for:

| Species | Source |
|---|---|
| H$_2$O | ExoMol; with the above-200 nm measurement of Ranjan et al. (2020) and Schulz et al. (2002) above 1500 K |
| CO$_2$ | Venot et al. (2018); 1160 K from ExoMol |
| NH$_3$ | ExoMol |
| O$_2$ | Frederick & Mentall (1982); Vattulainen et al. (1997) |
| SH, H$_2$S, OCS, CS$_2$ | Gorman et al. (2019) |

For H$_2$O the noisy data above 216 nm are fitted log-linearly as a conservative estimate.

## Interpolation scheme

The cross section of a given species is allowed to vary across the atmosphere with the local
temperature. Interpolation is **linear in temperature** and **logarithmic in the cross
section**. Because data are sparse, linear-in-$T$ interpolation tends to *under*estimate the
cross section, so the implementation is a conservative lower bound on how strongly photolysis
increases with temperature.

## Implementation in VULCAN

Temperature dependence is activated by listing species in [`config.T_cross_sp`](../Reference/config.md#photochemistry-radiation) (empty by
default). In [`op.ReadRate.make_bins_read_cross`](../Reference/api/op.md#vulcan.op.ReadRate.make_bins_read_cross):

- For each such species, per-temperature files `<sp>_cross_<T>K.csv` are read from
  `thermo/photo_cross/<sp>/`, and the room-temperature file is registered as the 300 K entry.
- For every layer, the routine finds the two bracketing tabulated temperatures
  $T_\mathrm{low}$ and $T_\mathrm{high}$ around the local $T_z$, interpolates
  $\log_{10}\sigma$ linearly between them, and stores the result in
  `cross_T[sp]` (shape $n_z \times n_\mathrm{bin}$) and the per-branch `cross_J_T[(sp, i)]`.
- Outside the wavelength range covered by the T-dependent data, the standard room-temperature
  cross section `cross[sp]` is used. Temperatures below the lowest / above the highest
  tabulated value fall back to the nearest available cross section.

These layer-resolved arrays then feed `compute_tau` (optical depth) and `compute_J`
(photolysis rates), where the T-dependent branch is taken for any species in `T_cross_sp`.

!!! note "Caveat"
    For high temperatures ($T > 1000$ K) the assembled data only cover wavelengths longer
    than about 190 nm, so the effect of temperature dependence in the FUV is not captured and
    may be larger than the current model implies. The CO$_2$ shielding effect on other species
    is real but, in the cases studied, was found to be masked once sulfur species are included.

[^tsai2021]: Tsai, S.-M., Malik, M., Kitzmann, D., et al. (2021). A comparative study of atmospheric chemistry with VULCAN. *The Astrophysical Journal, 923*(2), 264. https://doi.org/10.3847/1538-4357/ac29bc

[^heays2017]: Heays, A. N., Bosman, A. D., & van Dishoeck, E. F. (2017). Photodissociation and photoionisation of atoms and molecules of astrophysical interest. *Astronomy & Astrophysics, 602*, A105. https://doi.org/10.1051/0004-6361/201628742