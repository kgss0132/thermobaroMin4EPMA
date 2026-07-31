# thermobaroMin4EPMA (MATLAB) — Thermobarometric calculation for EPMA datasets

**thermobaroMin4EPMA** is a MATLAB-based, interactive workflow for applying published **geothermometers and geobarometers** to cation-calculated **electron probe microanalysis (EPMA) datasets**.
The program supports:

- geothermometry at a fixed pressure
- geobarometry at a fixed temperature
- geothermometry over a pressure range
- geobarometry over a temperature range

The workflow provides explicit control over:

- thermometer / barometer selection
- mineral or mineral-assemblage selection
- pressure and temperature units
- fixed-value and range calculations
- standardized pressure and temperature outputs
- reproducible Excel export and diagnostic plots for visually inspecting calculated pressure- or temperature-dependent trends

The calculation equations are implemented as modular MATLAB packages (`+thermo` and `+baro`), allowing individual calibrations to be inspected, debugged, replaced, and extended independently.

> Tested with **MATLAB R2024b**
>
> Input data should normally be prepared using **cationCalc4EPMA** or an equivalent cation-normalization workflow.

---

## Authors

- **Kazuki Matsuyama** (GSES, Nagoya University, Japan / Geoscience Montpellier, University of Montpellier, France)
- **Yumiko Harigane** (TUMSAT, Japan / GSJ, AIST, Japan)
- **Yoshihiro Nakamura** (GSJ, AIST, Japan)

---

## What this project does

Given a cation-calculated EPMA dataset, **thermobaroMin4EPMA**:

1. Imports an EPMA-derived dataset from an Excel or CSV file.
2. Asks the user to select one of four thermobarometric calculation modes.
3. Loads the available geothermometers or geobarometers from an Excel list.
4. Guides the user through interactive selection of:

    - thermometer or barometer group
    - published calibration
    - pressure or temperature unit
    - fixed value or calculation range
    - required mineral analyses or mineral assemblages

5. Converts user input to the internal standard units:

    - pressure: `kbar`
    - temperature: `degreeC`

6. Runs the selected thermometer or barometer module.
7. Standardizes the pressure and temperature variables required for downstream processing.
8. Adds calculation-condition metadata to the result table.
9. Exports the results to an Excel file.
10. Exports a diagnostic plot for visually inspecting the calculated pressure- or temperature-dependent trends in range calculations.

The main launcher, selection scripts, and individual calculation modules are separated so that published equations can be checked and updated independently.

---

## Why the `+thermo` and `+baro` folders?

Directories beginning with `+` are MATLAB packages. This design:

- enforces namespacing
  (`thermo.Pyroxene.BreyKohler1990`, `baro.Amphibole.RidolfiRenzulli2012`)
- prevents function-name conflicts
- groups calibrations by mineral or mineral assemblage
- separates geothermometers from geobarometers
- allows easy addition or replacement of individual equations

---

## Installation

### Requirements

- MATLAB R2024b (tested)
- Git (optional, for cloning)

No additional MATLAB toolboxes are required.

---

### Install via GitHub clone (recommended)

Clone the repository from GitHub:

```bash
git clone https://github.com/kgss0132/thermobaroMin4EPMA.git
```

Move into the project directory:

```bash
cd thermobaroMin4EPMA
```

---

### MATLAB setup

#### Option A: MATLAB Project (recommended)

If a `.prj` file is included:

1. Open MATLAB.
2. Double-click the `.prj` file, or select **Home → Project → Open**.
3. All paths are configured automatically.

This is the easiest and safest option for most users.

#### Option B: Manual path setup

If you are not using a MATLAB Project, add the project directory to your MATLAB path:

```matlab
addpath(genpath("path/to/thermobaroMin4EPMA"))
```

Then run:

```matlab
thermobaroMin4EPMA
```

---

## Inputs

### Cation-calculated EPMA dataset

Supported file formats are:

- Excel workbook (.xlsx)
- comma-separated file (.csv)

The input data should normally contain mineral compositions calculated as elemental cations or atoms per formula unit.
The exact required variables depend on the selected thermometer or barometer.

Typical data include:

- analysis identifier or dataCode
- Si, Ti, Al, Cr, Fe, Mn, Mg, Ca, Na, and K cations
- mineral-specific site occupancies or end-member components
- optional oxide wt% and metadata columns

The first column of each table is generally treated as the analysis identifier and is displayed during mineral selection.
Unique identifiers are recommended.

#### Excel input structure

For an Excel workbook, each worksheet is imported as a separate table and stored in rawdata_struct using the worksheet name.

For example:

```matlab
rawdata_struct.Olivine
```

Worksheet names should:

- be valid MATLAB structure-field names
- match the mineral names expected by the selected calculation module

An Excel workbook with separate mineral worksheets is recommended for multi-mineral thermometers and barometers.

#### CSV input structure

For a CSV file, the imported table is stored as:

```matlab
rawdata_struct.data
```

Because many calculation modules search for mineral-specific table names, CSV input is mainly suitable for methods that use a single table or explicitly support the data field.

---

### Thermometer and barometer list workbooks

The available calculation methods are controlled by:

- `functions/+thermo/Geothermometer_list.xlsx`
- `functions/+baro/Geobarometer_list.xlsx`

Each worksheet represents a thermometer or barometer group.
The worksheet name must match the corresponding MATLAB package folder.

Each method column uses the following four-row format:

- Row 1: Geothermometer / Geobarometer (Method name)
- Row 2: Reference (Published paper)
- Row 3: Type (Method classification)
- Row 4: MATLAB function name (Base filename without .m)

The selection dialog displays each method as:

`Reference (Method name)`

---

### Liquid-composition datasets

Some thermometers and barometers require the composition of an equilibrium liquid in addition to the cation-calculated mineral data.

Liquid-composition files are stored in:

`functions/+liquid/+liquidComp/`

A dataset based on **Kinzler, R.J. and Grove, T.L. (1992)** is included as a default example.

Users can add their own liquid-composition files to the same folder. When preparing a new file:

- use the included default file as a template
- preserve the required column names and data structure
- store the completed file in `functions/+liquid/+liquidComp/`

The required input format is fully compatible with files exported from **ExPetDB**, which incorporates the **Library of Experimental Phase Relations (LEPR)**, **TraceDs**, and **DCO-EPC** datasets. Files exported directly from ExPetDB, as well as independently prepared datasets that use the same column structure, can be placed in this directory and used in the calculations. This allows users to incorporate any required liquid-composition dataset into the program.

When the selected thermometer or barometer requires a liquid composition, the program automatically opens a selection dialog after the calculation module starts and asks the user which liquid-composition dataset should be used.

---

## Outputs

The program generates:

- an Excel file containing the complete calculation results
- a PNG diagnostic plot for range calculations, showing calculated temperature as a function of imposed pressure or calculated pressure as a function of imposed temperature

The result table may contain:

- selected analysis identifiers
- calculated temperature or pressure
- standardized `T_degreeC` and `P_kbar` variables
- converted `T_K` and `P_MPa` variables
- fixed-value or range metadata
- calibration-specific intermediate terms
- calibration-specific validity or warning information

### Calibration-specific warnings

Each thermometer and barometer module can implement validity checks based on the requirements and restrictions reported for that calibration. Depending on the method, these checks may include:

- mineral or liquid compositional limits
- imposed or calculated pressure and temperature ranges
- mineral–mineral or mineral–liquid equilibrium tests
- site-occupancy requirements
- water-content assumptions
- oxidation-state or ferric-iron assumptions
- restrictions associated with the experimental or empirical calibration dataset

Warning messages identify the affected analysis and the condition that was not satisfied. This allows users to distinguish calculations that satisfy the implemented criteria from those that approach or extend beyond a calibration limit and therefore require particular caution.

Because the relevant restrictions differ among published equations, warning criteria are defined separately within each thermometer or barometer module rather than by a single universal threshold. The warnings provide guidance on numerical and compositional applicability, but they do not demonstrate mineral equilibrium or establish the geological validity of a calibration for a particular sample.

---

### Output Excel file

The result table is exported to a worksheet named `Results` in an Excel workbook saved in the same directory as the input dataset.

The filename is constructed as:

`<input>_<mode>_<group>_<function>_<condition>.xlsx`

Example:

`Sample01_Geotherm_fixedP_Pyroxene_BreyKohler1990_15kbar.xlsx`

Decimal points in filename condition values are replaced with `p`.

Example:

`12.5 kbar` → `12p5kbar`

---

### Diagnostic plots for range calculations

Diagnostic plots are created only for range calculations.

Each diagnostic plot uses:

- horizontal axis: temperature in degreeC
- vertical axis: pressure in kbar

Each selected analysis or analysis combination is plotted as a separate curve.
Curves are distinguished using:

- changes in dataCode variables
- resets of the pressure or temperature range

The figure is exported as:

`<input>_<mode>_<group>_<function>_<condition>_TPplot.png`

These plots are intended for visual inspection of calculated pressure- or temperature-dependent trends. They do not combine independent thermometer and barometer solutions into a single equilibrium pressure–temperature diagram. If plot generation fails, the Excel result table is retained.

---

## How to run

### Default GUI mode

Simply call the main function:

```matlab
thermobaroMin4EPMA
```

The script will guide the user through:

- calculation-mode selection
- dataset selection
- thermometer or barometer selection
- pressure or temperature input
- mineral-analysis selection

No manual parameter editing is required for standard use.

---

### Calculation workflow

The main workflow is organized as follows:

- Stage 00: Initialize the program and add the `functions` folder to the MATLAB path
- Stage 01: Select the calculation mode
- Stage 02: Import the Excel or CSV dataset
- Stage 03: Load the thermometer or barometer list workbook
- Stage 04: Select the mineral group
- Stage 05: Select the published calibration
- Stage 06: Select the pressure or temperature unit
- Stage 07: Enter a fixed value or a minimum–maximum range
- Stage 08: Convert the input to `kbar` or `degreeC`
- Stage 09: Run the selected calculation module
- Stage 10: Standardize result variables and add metadata
- Stage 11: Export the result table
- Stage 12: Create and export the diagnostic trend plot for range calculations

---

### Calculation modes

#### 1. Geothermometer at fixed pressure

Launcher:

```matlab
startThermoCalc_fixedP
```

Available pressure units:

- `kbar`
- `MPa`

Conversion:

```text
1 kbar = 100 MPa
```

The selected module is called as:

```matlab
results = thermo.<group>.<mfile>(rawdata_struct, P_kbar)
```

The output table includes:

- `P_input_value`
- `P_input_unit`
- `P_MPa`
- `P_kbar`

#### 2. Geobarometer at fixed temperature

Launcher:

```matlab
startBaroCalc_fixedT
```

Available temperature units:

- `degreeC`
- `K`

Conversion:

```matlab
T_degreeC = T_K - 273.15;
```

The selected module is called as:

```matlab
results = baro.<group>.<mfile>(rawdata_struct, T_degreeC)
```

The output table includes:

- `T_input_value`
- `T_input_unit`
- `T_degreeC`
- `T_K`

#### 3. Geothermometer over a pressure range

Launcher:

```matlab
startThermoCalc_rangeP
```

The user enters minimum and maximum pressures in `kbar` or `MPa`.
The program creates 101 equally spaced pressure values:

```matlab
P_kbar_range = linspace(PMinKbar, PMaxKbar, 101).';
```

The complete pressure vector is passed to the selected thermometer module:

```matlab
results = thermo.<group>.<mfile>(rawdata_struct, P_kbar_range)
```

For diagnostic plotting, the output should contain:

- `T_degreeC`
- `P_kbar`

The launcher can also recognize selected alternative temperature-variable names and convert them to `T_degreeC`.

#### 4. Geobarometer over a temperature range

Launcher:

```matlab
startBaroCalc_rangeT
```

The user enters minimum and maximum temperatures in `degreeC` or `K`.
The program creates 101 equally spaced temperature values:

```matlab
T_degreeC_range = linspace(TMinDegreeC, TMaxDegreeC, 101).';
```

The complete temperature vector is passed to the selected barometer module:

```matlab
results = baro.<group>.<mfile>(rawdata_struct, T_degreeC_range)
```

For diagnostic plotting, the output should contain:

- `T_degreeC`
- `P_kbar`

The launcher can also recognize selected alternative pressure-variable names and convert them to `P_kbar`.

---

### Adding a new thermometer or barometer

#### New thermometer

Place the function in:

`functions/+thermo/+<group>/<mfile>.m`

Use the standard function interface:

```matlab
function results = <mfile>(rawdata_struct, P_kbar)
```

The function should:

- accept scalar pressure for fixed-pressure calculations
- accept vector pressure for pressure-range calculations
- return a MATLAB table
- include `T_degreeC` and `P_kbar`

Add the method to the corresponding worksheet in `Geothermometer_list.xlsx`.

#### New barometer

Place the function in:

`functions/+baro/+<group>/<mfile>.m`

Use the standard function interface:

```matlab
function results = <mfile>(rawdata_struct, T_degreeC)
```

The function should:

- accept scalar temperature for fixed-temperature calculations
- accept vector temperature for temperature-range calculations
- return a MATLAB table
- include `T_degreeC` and `P_kbar`

Add the method to the corresponding worksheet in `Geobarometer_list.xlsx`.

---

### Applicability and interpretation

**thermobaroMin4EPMA** performs numerical calculations but does not determine whether a calibration is geologically appropriate. Method-specific warnings indicate whether the implemented numerical, compositional, or experimental criteria are satisfied; they should not be interpreted as proof of mineral equilibrium or geological validity.

Before interpreting the results, confirm that:

- the required mineral phases are present
- the selected analyses represent a plausible equilibrium assemblage
- mineral zoning, exsolution, alteration, and retrogression have been evaluated
- the mineral compositions fall within the published calibration range
- the imposed and calculated pressure and temperature conditions lie within the experimental or empirical range
- ferric-iron, site-allocation, water-content, oxidation-state, and activity-model assumptions are appropriate
- calibration uncertainty is considered

Results outside the published calibration range should be treated as extrapolations, even when the software returns finite numerical values.

---

### Reproducibility

For each calculation, retain:

- the original EPMA dataset
- the cation-calculated input workbook
- the selected thermometer or barometer
- the selected mineral analyses
- the fixed pressure or temperature
- the specified pressure or temperature range
- the exported Excel result table
- the exported diagnostic plot
- MATLAB warnings and error messages
- the original published calibration

---

## Troubleshooting

### The main function is not found

Add the project directory to the MATLAB path:

```matlab
addpath(genpath("path/to/thermobaroMin4EPMA"))
```

Then run:

```matlab
thermobaroMin4EPMA
```

---

### The thermometer or barometer list workbook is not found

Confirm that these files exist:

- `functions/+thermo/Geothermometer_list.xlsx`
- `functions/+baro/Geobarometer_list.xlsx`

---

### A selected module cannot be found

Confirm that:

- the worksheet name matches the package-folder name
- the package folder begins with `+`
- Row 4 of the list workbook matches the MATLAB filename
- the MATLAB function declaration matches the filename

---

### A mineral table is missing

Confirm that:

- the Excel worksheet name matches the mineral name expected by the module
- the selected method supports the imported dataset structure
- the required mineral compositions were exported from `cationCalc4EPMA`

Read the header comments in the selected thermometer or barometer module for accepted table names and aliases.

---

### No results are exported

Possible causes include:

- file selection was canceled
- group or method selection was canceled
- pressure or temperature input was canceled
- no mineral combination was selected
- the selected module returned an empty table
- invalid input stopped the calculation

Check the MATLAB Command Window for the last displayed message.

---

### The Excel file is exported but the diagnostic plot is not created

Range plotting requires finite numeric variables named:

- `T_degreeC`
- `P_kbar`

The Excel result table is retained even when plot generation fails.

---

### Results contain `NaN`, `Inf`, or unrealistic values

Check:

- missing cation values
- negative or infinite input values
- zero denominators
- invalid logarithm arguments
- mineral-pair equilibrium
- mineral-site allocation
- pressure and temperature units
- calibration limits
- transcription of the original published equation

---

## Relationship to cationCalc4EPMA

### cationCalc4EPMA

Converts EPMA oxide wt% data into elemental cations and atoms per formula unit (a.p.f.u.).

### thermobaroMin4EPMA

Applies published thermometer and barometer equations to cation-calculated mineral compositions.

### Typical workflow

1. Raw EPMA oxide data
2. `cationCalc4EPMA`
3. Cation-calculated mineral tables
4. `thermobaroMin4EPMA`
5. Thermobarometric estimates

---

## Tests

### Running the test suite

The regression test is located in the `test/` directory.

Run the test from MATLAB as follows:

```matlab
cd test
test_cationCalc4epma
```

Because the calculation settings and mineral analyses are selected through interactive dialogs, the required selections must be made manually as described below.

---

### Input dataset

The test uses the following input workbook:

- `expected_output_Results_NaN2zero.xlsx`

This workbook is based on the `Results` output generated by the `cationCalc4EPMA` regression test. Before it is used by `thermobaroMin4EPMA`, all `NaN` values in the workbook are replaced with `0` for this test.

The workbook must be placed directly in the `test/` directory.

---

### Required user selections during testing

When the dialogs are displayed, select the following settings exactly.

#### 1. Calculation mode

Select: `Thermometer: pressure range`

This mode calculates temperature over a user-defined pressure interval.

#### 2. Input data file

Select: `expected_output_Results_NaN2zero.xlsx`

The file-selection dialog should normally open in the `test/` directory.

#### 3. Mineral group

Select: `Pyroxene`

This opens the list of geothermometers available for pyroxene compositions.

#### 4. Thermometer

Select: `Brey_and_Kohler_1990 (Opx Ca)`

Depending on how the thermometer name is displayed in the selection table, the corresponding entry may internally be associated with `BreyKohler1990Ca`. This thermometer calculates temperature from the Ca content of orthopyroxene.

#### 5. Pressure unit

Select: `kbar`

#### 6. Minimum pressure

Enter: `10`

#### 7. Maximum pressure

Enter: `30`

The calculation is therefore performed over a pressure range from 10 to 30 kbar.

The program converts this interval into the pressure values used by the pressure-range calculation and evaluates the selected thermometer at each pressure.

#### 8. Mineral data

Select: `test_003`

Only the `test_003` analysis should be selected for this regression test. Selecting additional analyses will generate additional calculation results and diagnostic curves, causing the output files to differ from the reference files.

---

### Expected output files

A successful calculation generates the following two files in the `test/` directory.

#### Calculation results

`expected_output_Results_NaN2zero_Geotherm_rangeP_Pyroxene_BreyKohler1990Ca_P10to30kbar.xlsx`

This workbook contains the temperature calculations for `test_003` over the pressure range from 10 to 30 kbar.

Due to the output structure of `thermobaroMin4EPMA`, the workbook includes:

- the selected mineral identifier
- pressure values
- calculated temperatures
- input compositional data
- calculation metadata
- thermometer-specific parameters

#### Diagnostic trend plot

`expected_output_Results_NaN2zero_Geotherm_rangeP_Pyroxene_BreyKohler1990Ca_P10to30kbar_TPplot.png`

This image shows the calculated temperature response to the imposed pressure range for the selected `test_003` orthopyroxene analysis.

---

### Successful test result

When both generated files match their reference files, MATLAB displays:

```matlab
==================================================================
All thermobaroMin4EPMA outputs are consistent with expected files!
==================================================================
```

This confirms that:

- the expected output files were generated
- the Excel workbook structure is unchanged
- the calculated values match the validated results
- the exported diagnostic plot matches the validated image

---
