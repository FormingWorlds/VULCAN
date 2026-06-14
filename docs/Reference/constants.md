# Physical constants

Constants used throughout VULCAN, defined in [`phy_const.py`](api/phy_const.md). All values are in **cgs units**,
taken from the [astropy package](http://docs.astropy.org/en/stable/index.html).

| Symbol | Name | Value | Units |
|---|---|---|---|
| `kb` | Boltzmann constant | `1.38064852e-16` | erg K$^{-1}$ |
| `Navo` | Avogadro's number | `6.02214086e23` | mol$^{-1}$|
| `hc` | Planck constant × speed of light | `1.98644582e-9` | erg nm |
| `au` | Astronomical unit | `1.49597871e13` | cm |
| `r_sun` | Solar radius | `6.957e10` | cm |
| `r_jup` | Jupiter equatorial radius | `7.1492e9` | cm |
| `ag0` | Asymmetry factor (radiative transfer) | `0` | — |

!!! note "Gas constant $R$"
    The universal gas constant is $R = k_B \times N_\text{avo}$, and the specific gas constant for a species with molar mass $M$ is $R_\text{spec} = R / M$.
    `hc` is used to convert the stellar flux to actinic flux (photons cm$^{-2}$ s$^{-1}$ nm$^{-1}$).