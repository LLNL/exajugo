## Overview
This project converts MATPOWER case files into the `.raw`, `.rop`, `.con`, and `.inl` formats required for running **ExaJuGo**.  
The provided script automatically generates these files starting from a MATPOWER `.m` case file.

For the **California Test System**, you can also create **wildfire-dependent `.con` files** by using wildfire data.  
Instructions, along with information about wildfire data types and sources, are available at:  
https://github.com/SLOPE-grid/sc-acopf-viz/tree/main

Once you have the wildfire data, you can generate `.con` files using the notebook located at:  
`preprocessing_4_ExaJuGo/Make_con_file.ipynb` in that repository.

## Requirements
- **MATLAB** 
- **MATPOWER** installed (https://matpower.org/)

**Important:** You must override the default `save2psse.m` function in MATPOWER by using the customized version provided in this directory.  
This ensures compatibility with the output format expected by ExaJuGo.

## Files and Directories
- **`example_m_files/`**  
  Contains example MATPOWER case files:
  - `CaliforniaTestSystem.m`: California test system used to generate ExaJuGo input files.
  - `case9_mod.m`: Modified version of the standard 9-bus system.

- **`California/`**  
  Directory where the generated `.raw`, `.rop`, `.con`, and `.inl` files will be stored.

- **`main_script.m`**  
  The main MATLAB script that:
  - Loads a case from `example_m_files/`
  - Converts and saves the output files into `California/`

## Moving Files to ExaJuGo
To run the ExaJuGo code with the generated California files, move the `California/` folder into the `examples/` directory of ExaJuGo.  
You can do this by running the following terminal command from the current directory:

```bash
mv ./California ./../examples
```

## Additional Examples
If you would like to use more MATPOWER case files, you can find many official examples here:  
https://github.com/MATPOWER/matpower/tree/master/data

You can download a `.m` case file and place it inside the `example_m_files/` directory, then update the script to use the new filename.

## Notes
- The `California/` folder contains the finalized input files needed to run ExaJuGo with the **California test system**.
- The `.con` and `.inl` files are created as **empty placeholders** (ExaJuGo requires them even if not populated).

## Quick Start

1. Install [MATPOWER](https://matpower.org/) and make sure it is added to your MATLAB path.
2. Replace the `save2psse.m` file in MATPOWER with the custom version from this directory.
3. Run the main script in MATLAB:

```matlab
create_ExaJuGo_files
```

4. Move the generated `California/` folder into the ExaJuGo `examples/` directory:

```bash
mv ./California ./../examples
```

5. Now you are ready to run ExaJuGo using the California test system!

### Optional Step: Create Wildfire-Dependent `.con` File

If you want to include wildfire-dependent data, follow these additional steps:

1. Visit the following link for instructions on wildfire data types and sources: [Wildfire Data in SC-ACOPF](https://github.com/SLOPE-grid/sc-acopf-viz/tree/main).
2. With the wildfire data, use the code in `preprocessing_4_ExaJuGo/Make_con_file.ipynb` to generate the `.con` file.
   
This will create a `.con` file that accounts for wildfire impacts on the California grid.

6. **Save the generated files** to the directory you moved into ExaJuGo (`./../examples/California/`).

---


