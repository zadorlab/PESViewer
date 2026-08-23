# PESViewer
#### Potential energy surface visualizer

This code was originally developed by [Ruben Van de Vijver](https://github.com/rubenvdvijver). 
Our fork has been developed and improved on top of the original to maintain compatibility with KinBot and integrate new features.

PESViewer is a code to depict and analyze a potential energy surface 
characterized by wells, bimolecular products, transition states and barrierless reactions.
<!--
Their energy is needed to plot the potential energy surface.)
Written values of the energies and 2D plots (structural formulas) of the wells and products can be added to the figure.)
-->

To run PESViewer, you need python >= 3.7, matplotlib, numpy, and both RDKit and OpenBabel to create 2D plots.

While the code can be used as a standalone application, it is designed to work with KinBot, which automatically generates the input files for PESViewer.

## How to Install

PESViewer can be installed both in two different ways, from the conda-forge repo (`conda install`) or by cloning this github repo and then install it locally.

### conda-forge

    conda install -c conda-forge pesviewer

### From Github

If you want to have the very last version of PESViewer without waiting for a 
release or you want to modify PESViewer acccording to your needs you can clone the project 
from github:

    git clone git@github.com:zadorlab/PESViewer.git

and then, from within the PESViewer directory produced after cloning, type:

    pip install -e .

> **Note**
> Open Babel is a required dependency but is not declared in `pyproject.toml`,
> because its PyPI distribution is source-only and needs the Open Babel C++
> library and SWIG to build. Install it separately with conda:
>
>     conda install -c conda-forge openbabel

> **Note**
> If you want to modify PESViewer yourself it's better to fork the project 
> into your own repository and then clone it.

## INPUT

All the information needed for PESViewer to make plots and graph depictions of the PES is provided through a text input file.
In this file, stationary points, energies, names identifiers and options are listed in different sections. 
An example of an input file is provided in `input.txt`.

### Stationary Points
Each stationary point is written on a separate line, under the the relevant section. 

For wells, the name and energy of each unimolecular species must be specified under the `> <wells>` header, as follows. Additionally, a SMILES code can be provided alongside the name and energy. This can be convenient when no xyz is present, OpenBabel fails to interpret connectivity correctly or a particular resonant structure is desired:

    name energy [smiles]

Similarly, for bimolecular products each pair of species is written on a separate line under the `> <bimolec>` header:

    species1_species2 energy [smiles]

For reactions occurring via a transition state, the name, energy, reactant and product of the transition state must be specified under the `> <ts>` section. Additionally, a color can be provided to highlight each reaction:

    name energy reactant product [color]

For barrierless reactions, a similar approach is followed. The `> <barrierless>` header indicatess reactions occurring without a barrier. The energy is not needed in this case:

    name reactant product [color]

### Options
The plotting options (to be written in the input file) are listed here. 


| option | default | Description | Scope* |
| ------- | ------- | ------- | ------- |
| `title` | `0` | print a title (1) or not (0) on the PES. | T |
| `units` | `kJ/mol` | energy units in the input file. | TG | 
| `display_units` | same as the units of `units`  | energy units to be displayed. Allowed values: `kJ/mol`, `kcal/mol`, `eV` | TG | 
| `rounding` | `1`  | number of decimals for the energy values | TG | 
| `energy_shift` | `0.`  | shift energy scale by this amount measured in `units` | TG | 
| `use_xyz` | `1` | use xyz files to generate each species 2D depiction. The `*.xyz` files should be named the same as the name specified in the relevant line and placed inside a directory called `xyz` (0 switch it off). | TG | 
| `rescale` | `0` | which species is used to rescale all energies, the name of the relevant well or bimolecular species is needed. | TG |
| `fh` | `9.` | figure height. | T | 
| `fw` | `18.` | figure width. | T |
| `fs` | `1` | scale factor for species 2D depictions. The images generated under the `*_2d` directory need to be deleted for them to be regenerated with a new scaling. | T |
| `lw` | `1.5` | line width of the potential energy vs reaction coordinate plot. | T |
| `margin` | `0.2` | margin fraction on the x and y axis | T |
| `dpi` | `120` | Image resolution in dpi of the molecule figures | TG |
| `save` | `0` | does the plot need to be saved a an image (1) or displayed to be interactively modified (0). | T |
| `plot` | `1` | whether to plot (1) or not (0) the potential energy vs reaction coordinate representation of the PES. | T |
| `write_ts_values` | `1` | whether the TS energy values should be written in the plot. | T |
| `write_well_values` | `1` | whether the well and bimolecular energy values should be written in the plot. | T |
| `bimol_color` | `red` | color of the energy values for the bimolecular products. | T |
| `well_color` | `blue` | color of the energy values of the wells. | T |
| `ts_color` | `green` | color or the energy values of the ts, 'none' to use same color as line. | T |
| `show_images` | `1` | boolean tells whether the molecule images should be shown on the plot. | T |
| `rdkit4depict` | `1` | boolean that specifies which code to use for the 2D depiction. | TG |
| `axes_size` | `10` | font size of the axes. | T |
| `text_size` | `10` | font size of the energy values on the plot. | T |
| `linear_lines` | `0` | plot polynomials (0) or linear lines (1) between de stationary points. | T |
| `interpolation` | `hanning` | image interpolation method. | T |
| `graph_edge_color` | `black` | color of the graph edges, if set to `energy`, color is scaled accordingly. | G |
| `reso_2D` | `1` | enable (1) or disable (0) generation of 2D depictions for resonant structures. Additional images, one for each resonant structure named `_X` are generated under the `_2d` directory. | TG |
| `path_report` | `[]` | Compute the MEP between two species. Format: `chemid_start chemid_end` | TG |
| `search_cutoff` | `10` | Maximum length (in reactive steps) to search for of the path report. | TG |
| `node_size_diff` | `0` | Size difference of nodes in the graph depiction based on to their stability. More stable nodes are larger. 0 to make all nodes the same size. Reasonable values: 20-40.| G |

* This column shows whether the parameter impacts the traditional (T) Potential vs Reaction coordinate PES depiction, the graph (G) one, or both (TG).

The other input is a folder called `xyz/` containing the coordinates of the stationary points in xyz format (`name.xyz`)
(for bimolecular products, use several xyz coordinates files, `name_index.xyz`). These 


## RUN

With the input file `input.inp`, type:

    pesviewer input.inp

## OUTPUT

When the code runs and the depiction of the molecules is requested, first, all depictions are generated based on the content of the `/xyz` folder. On repeated runs the depictions are reused, unless missing. 

**Traditional PES depiction** 

The output is a modifiable matplotlib figure, which can be displayed and interactively arranged, or saved. 
The possible modifications are:
- modifing the x-position of a stationary point by draggning the energy value
- modifing the position of 2D structure images by dragging the image
- zooming in an out.

Helper files with `.txt` extension are also generated saving the positions of stationary points and images, which, in certain cases need to be deleted to recreate the plot on subsequent runs (e.g., when changing the `fs` parameter). However, the information about the adjustmens made are stored here.

By selecting a stationary point, all direct neighbors and pathways are lit up, and the others are dimmed to help navigation.

This is an example of a nicely arranged traditional PES plot from our recent paper:

![image](https://user-images.githubusercontent.com/40675474/227331800-373cf4b7-5d17-4f7a-8347-06544badc5b8.png)

* Martí, C., Michelsen, H. A., Najm, H. N., Zádor, J.: _Comprehensive kinetics on the C7H7 potential energy surface under combustion conditions._ J. Phys. Chem. A, **2023**, 127, 1941–1959. https://pubs.acs.org/doi/full/10.1021/acs.jpca.2c08035  

**Interactive graph representation**

Here, wells (black-border circles) and bimolecular products (blue-border circles) are shown as nodes of a graph. Their energies are displayed too, and the name of the species can be read by hovering over their depiction. Transition states are shown as edges. Edges representing saddle points are black (unless specified), barrierless reaction pathways are gray (unless specified). The thickness of the edge is inversely proportional to the absolute height of the barrier, and the energy can be read when hovering over the edge. Optionally, the color of the edges can be controlled using the 

This representation uses the [pyvis](https://pyvis.readthedocs.io/en/latest/) package, generates an `html` file, which can be opened in a browser `(file://<path/to/file.html>)`. Note the settings below the graph to help arrange the graph using intuitive physics analogies, which can also be turned off.

This is a static screenshot example from the same publication - follow this link https://pubs.acs.org/doi/full/10.1021/acs.jpca.2c08035 for the Supporting Information to download interactive plots like this.

![image](https://user-images.githubusercontent.com/40675474/227336672-c7448207-fc3b-42c3-ad89-7da45f84e985.png)

