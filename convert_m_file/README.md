## Overview
This repository converts MATPOWER case files into the `.raw`, `.rop`, `.con`, and `.inl` formats required for running **ExaJuGO**.  
The provided script automatically generates these files starting from a MATPOWER `.m` case file.

- `.raw`: Power system network topology and parameters  
- `.rop`: Operating point information (e.g., voltage, generation, demand)  
- `.con`: Contingency definitions  
- `.inl`: Generator inclusion list (e.g., which generators are online)

For the **California Test System**, you can also create **wildfire-dependent `.con` files** by using wildfire data.  
Instructions, along with information about wildfire data types and sources, are available at:  
https://github.com/SLOPE-grid/sc-acopf-viz/tree/main

Once you have the wildfire data, you can generate `.con` files using the notebook located at:  
`preprocessing_4_ExaJuGO/Make_con_file.ipynb` in that repository.

## Requirements
- **MATLAB**
- **MATPOWER** installed ([https://matpower.org/](https://matpower.org/))

**Important:**  
You must use the customized `save2psse.m` function provided in this directory to ensure compatibility with ExaJuGO.  
You **do not** need to overwrite the version inside MATPOWER. Simply run the main script from the current directory, and MATLAB will use the local version automatically.

## Files and Directories
- **`example_m_files/`**  
  Contains example MATPOWER case files:
  - `CaliforniaTestSystem.m`: California test system from [CATS-CaliforniaTestSystem](https://github.com/WISPO-POP/CATS-CaliforniaTestSystem/tree/master).
  - `case9_mod.m`: Modified version of the standard 9-bus system.

- **`California/`**  
  - Directory where the generated `.raw`, `.rop`, `.con`, and `.inl` files for the California case are stored.
  - This directory name should reflect the test case name.
  - For example, if using `case9_mod.m`, you should name the directory something like `9bus/`.

- **`main_script.m`** (file: `create_ExaJuGO_files.m`)  
  The main MATLAB script that:
  - Loads a case from `example_m_files/`
  - Converts and saves the output files into the appropriate output folder

## Moving Files to ExaJuGO
To run the ExaJuGO code with the generated case files, move the output directory (e.g., `California/`) into the `examples/` directory of ExaJuGO.  
This is **required by ExaJuGO**, which expects test cases to be located in that folder.

You can do this by running the following terminal command from the current directory:

```bash
mv ./California ./../examples
```

If using a different test case (e.g., `case9_mod.m`), replace `California` with the appropriate folder name.

## Additional Examples
If you would like to use more MATPOWER case files, you can find many official examples here:  
https://github.com/MATPOWER/matpower/tree/master/data

To use a different case:
1. Download a `.m` case file and place it inside the `example_m_files/` directory.
2. Update `main_script.m` (`create_ExaJuGO_files.m`) to load the new file.
3. Adjust the output folder name to match the new test case name (e.g., `9bus/`).

## Notes
- The generated folder (e.g., `California/`) contains the finalized input files needed to run ExaJuGO.
- The `.con` and `.inl` files are created as **empty placeholders** (ExaJuGO requires them even if not populated).

## Quick Start

1. Install [MATPOWER](https://matpower.org/) and make sure it is added to your MATLAB path.
2. Run the main script `main_script.m` (`create_ExaJuGO_files.m`) in MATLAB:

```matlab
create_ExaJuGO_files
```

3. Move the generated output folder (e.g., `California/`) into the ExaJuGO `examples/` directory:

```bash
mv ./California ./../examples
```

4. Now you are ready to run ExaJuGO using the selected test system!

### Optional Step: Create Wildfire-Dependent `.con` File

If you want to include wildfire-dependent data **(only applicable to the California test system)**, follow these additional steps:

1. Visit the following link for instructions on wildfire data types and sources:  
   [Wildfire Data in SC-ACOPF](https://github.com/SLOPE-grid/sc-acopf-viz/tree/main)

2. With the wildfire data, use the notebook at:  
   `preprocessing_4_ExaJuGO/Make_con_file.ipynb`  
   to generate a wildfire-aware `.con` file for the California grid.

3. Replace the existing `.con` file in your output folder (e.g., `./../examples/California/`) with the wildfire-dependent version.
