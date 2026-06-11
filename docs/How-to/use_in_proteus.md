# Use VULCAN in PROTEUS

This page walks you through enabling and configuring the VULCAN atmospheric chemistry module in a PROTEUS config file. For a description of how VULCAN is wired into the framework, see the [PROTEUS coupling page](../Explanations/proteus_coupling.md). 

!!! info "Working PROTEUS + VULCAN installation"
    This page assumed you have installed PROTEUS by following its [installation guide](https://proteus-framework.org/PROTEUS/How-to/installation.html) and [installed VULCAN as an optional dependency](https://proteus-framework.org/PROTEUS/How-to/installation.html#atmospheric-chemistry-vulcan). 

---

## 1. Enable the module

VULCAN can be enabled in your PROTEUS input configuration file. In the TOML, set `module` to `"vulcan"` under `[atmos_chem]`:

```toml
[atmos_chem]
    module = "vulcan"
```

The other available values are `"none"` (chemistry disabled, the default) and `"dummy"` (placeholder for testing).

!!! warning "AGNI required"
    VULCAN only works with `atmos_clim.module = "agni"`. Using any other atmosphere module will disable chemistry with a warning.

---

## 2. Choose when chemistry runs

The `when` key controls the scheduling:

```toml
[atmos_chem]
    when = "offline"   # manually | offline | online
```

| Value | When it runs |
|---|---|
| `"manually" ` | Never; chemistry is skipped, and will be done manually.   |
| `"offline"`   | Once, after the main simulation loop completes. Good for post-processing a finished run. |
| `"online"`    | At every data-write snapshot during the simulation. Output files are labelled by simulation year. |

---

## 3. Choose a chemical network

```toml
[atmos_chem.vulcan]
    network = "SNCHO"   # CHO | NCHO | SNCHO
```

Pick the network that covers the elements relevant to your simulation:

| Value | Elements included | Use when |
|---|---|---|
| `"CHO"` | H, O, C | Carbon–hydrogen–oxygen atmospheres; fastest. |
| `"NCHO"` | H, O, C, N | Also tracking nitrogen (e.g. NH$_3$, N$_2$). |
| `"SNCHO"` | H, O, C, N, S | Full sulfur chemistry (SO$_2$, H$_2$S, S$_2$); always uses the photo network. |

---

## 4. Set the initial composition

```toml
[atmos_chem.vulcan]
    ini_mix = "profile"   # profile | outgas
```



`"profile"` initialises each species from the altitude-resolved VMR profile written by AGNI. Species absent from the profile are initialised to zero.`"outgas"` uses the bulk outgassed composition as a uniform starting point.

---

## 5. Configure transport processes

```toml
[atmos_chem]
    photo_on      = true    # photochemistry
    Kzz_on        = true    # eddy diffusion
    Kzz_const     = "none"  # constant Kzz [cm2/s]; "none" = use AGNI-derived profile
    moldiff_on    = true    # molecular diffusion
    updraft_const = 0.0     # constant updraft velocity [cm/s]
```

To use a fixed eddy diffusivity rather than the profile computed by AGNI, set `Kzz_const` to a value in cm$^2$ s$^{-1}$, for example:

```toml
    Kzz_const = 1e5
```

---

## 6. Set convergence and numerical options

```toml
[atmos_chem.vulcan]
    yconv_cri   = 0.05      # steady-state criterion on mixing ratio change
    slope_cri   = 0.0001    # steady-state criterion on rate of change
    clip_fl     = 1e-20     # stellar flux floor [erg/s/cm2/nm]
    clip_vmr    = 1e-10     # ignore species below this VMR
```

The defaults work for most cases. Tighten `yconv_cri` and `slope_cri` if you need stricter convergence; loosen them to speed up runs.

---

## 7. Optional: fix the surface boundary condition

To hold surface mixing ratios fixed throughout the chemistry integration (useful for sensitivity tests):

```toml
[atmos_chem.vulcan]
    fix_surf = true
```

---

## 8. Minimal working example

A complete `[atmos_chem]` block for an offline SNCHO run:

```toml
[atmos_chem]
    module        = "vulcan"
    when          = "offline"
    photo_on      = true
    Kzz_on        = true
    Kzz_const     = "none"
    moldiff_on    = true
    updraft_const = 0.0

    [atmos_chem.vulcan]
        network     = "SNCHO"
        ini_mix     = "profile"
        fix_surf    = false
        make_funs   = true
        yconv_cri   = 0.05
        slope_cri   = 0.0001
        clip_fl     = 1e-20
        clip_vmr    = 1e-10
        save_frames = false
```