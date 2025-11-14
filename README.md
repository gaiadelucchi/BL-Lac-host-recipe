# BL-Lac-host-recipe

This repository contains the code to reproduce the analysis discussed in the paper:

**BL Lac host galaxies: how to systematically characterise them in optical-NIR spectroscopy**
*Gaia Delucchi, Tullia Sbarrato, Giorgio Calderone, Chiara Righi, Silvano Tosi, and Boris Sbarufatti*
Submitted to ApJ.


## Workflow

- Download and install [Julia](https://julialang.org/downloads/);
- Download the current code [here](https://github.com/gaiadelucchi/BL-Lac-host-recipe/archive/refs/heads/main.zip);
- Download the input data for the analysis [here](https://drive.google.com/drive/folders/1l1-9qVmnqUeuIWh-vGAVdVeBzaGZUD65?usp=sharing)
- Type the following commands:
  - `unzip BL-Lac-host-recipe-main.zip`;
  - `cd BL-Lac-host-recipe`;
  - `julia --project=. --startup-file=no`;
- From within the Julia session install the dependencies with:
  - `using Pkg`
  - `Pkg.instantiate()`;

Note: the above steps are necessary only the first time the code is executed.
  
Finally, you can run the analysis by typing `include("QSFit_analysis_df.jl")` from within the Julia session.
