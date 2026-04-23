
## Configuration File
<strong>All the settings and parameters, e.g. the atmospheric parameters, the elemental abundance etc, are prescribed in ```config.py```</strong>. Typically this is the only file you need to edit for each specific run. A useful cheatsheet describing what every parameter does can be found in ```cfg_examples/vulcan_cfg_Earth.py```. The configuration files used for the model validation in [Tsai et al. 2021](https://arxiv.org/abs/2108.01790) are also provided in the cfg_examples folder.

## Input Files
The key input files of VULCAN include the chemical network, atmospheric T-P profile, and stellar flux. ```NCHO_photo_network.txt``` is the deafult reaction network including nitrogen, carbon, hydrogen, and oxygen species. It is validated from ~ 500 to 3000 K with about 60 gaseous species and 700 reactions.
The rate coefficients A, B, C are written in A, B, C as in the Arrhenius formula k = A T^B exp(-C/T).
The input temperature-pressure(-Kzz) profile is required when Kzz_prof is set to 'file' in `config.py` and is placed in the `atm/` folder by default. The first line in the T-P file is commented for units, and the second line must specifies the column names: **Pressure	Temp** or **Pressure	Temp  	Kzz** (Kzz is optional). So the file consists of two columns without K<sub>zz</sub> and three columns with K<sub>zz</sub>.
See the included T-P files of HD 189733b and HD 209458b in `atm/` for example.
The stellar UV flux is stored in `atm/stellar_flux`, with the first column being weavelength in nm and the second column	being flux in ergs/cm**2/s/nm.
The thermodynamics data and cross sections are stored in `thermo/NASA9/` and `thermo/photo_cross/`, respectively. Change at your own risk!
If constant fluxes for certain species are used, the files are also placed in /atm, in the format of species, flux (cm-2 s-1), and deposite velocity (cm s-1).

## Editing or Using a different chemical network
VULCAN is developed in a flexible way that the chemical network is _not_ hard coded. Instead, ```make_chem_funs.py``` generates all the required funtions from the input chemical network (e.g. ```NCHO_photo_network.txt```) into ```chem_funs.py```.
You can edit the default netowrk, to remove or add reactions, to change rate constats, etc. You can also use a different chemical network, as long as it is in the same format as the defalut ones. That is, the reactions should be writen in the form of [ A + B -> C + D ], including the square brackets.
By default, ```make_chem_funs.py``` is always called prior to the main code to produce ```chem_funs.py``` based on the new chemical network . This step (which takes a few seconds) can be skiped by adding the agument ```-n```while running vulcan in the command line:
```
python vulcan.py -n
```
However, it is important NOT to skipping this step after making a change of the chemical network.

Noted that changing or using a different chemical network is not foolproof -- unrealistic values could lead to numerical issues. You can also see that only the forward reactions are listed, since VULCAN reverses the forward reactions to obtain the reverse reactions using the thermodynamic data.
So next, make sure all the species are included in the ```NASA9``` folder. If not, they need to be added manually by looking over ```nasa9_2002_E.txt``` or ```new_nasa9.txt```, which can also be found in ```NASA9```. Save the coefficients in a text file with the same name as used in the network (e.g. CO2.txt). The format of the NASA 9 polynomials is as follows
```
a1 a2 a3 a4 a5
a6 a7 0. a8 a9
```
Here, a7 and a8 are separated by 0. The first two rows are for low temperature (200 - 1000 K) and the last two rows are for high temperature (1000 - 6000 K).\

The reaction number, i.e. **id**, is irrelevent as it will be automatically generated (and writing into the network file) while calling ```make_chem_funs```. Three-body or dissociation reactions should also be separately listed after the comment line as the default network.
After changing the network, you can examine all the readable information, like the list of reactions and species in ```chem_funs.py```, being updated while running python vulcan.py (without -n argument).

## Boundary Conditions
If both use_topflux and use_botflux in `config.py` are set to False, it will use the default boundary condition -- zero flux boundary i.e. nothing coming in or out. When use_topflux = True, it reads the file prescribed in top_BC_flux_file as the incoming/outgoing flux at the top boundary. Similarly, when use_botflux = True, the file prescribed in bot_BC_flux_file is read in for the surface pressure and sinks at the bottom boundary. In addition, you can also use the dictionary use_fix_sp_bot to set fixed mole fraction at the surface. e.g. use_fix_sp_bot = {'CO2': 0.01} sets the surface CO<sub>2</sub> mixing ratio to 0.01.

## Coupling to a self-consistent climate model
VULCAN can be self-consistently coupled to a atmosphere climate model (AGNI) which solves for the atmospheric temperature profile given the composition calculated by VULCAN. AGNI uses correlated-k radiative transfer and mixing-length convection in order to determine realistic TP- and Kzz-profiles. More information on AGNI can be found [here](https://www.h-nicholls.space/AGNI/).

To enable chemistry-climate calculations, install AGNI on your machine and make it available inside the VULCAN folder (e.g. `VULCAN/AGNI/agni.jl`). Then set `agni_call_frq` to a value >0 in the VULCAN config object and run VULCAN like normal.

