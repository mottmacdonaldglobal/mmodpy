## JavaScript Script Module

This module is a collection of JavaScript files for automation in Oasys Primer or d3Plot.

## Usage
1. Install `mmodpy` library with pip (or download desired files from GitHub) <br />
`pip install mmodpy`
2. Open Oasys Primer and navigate to the `Script` window <br />
3. Select `Run` and open the desired JavaScript file when prompted <br />
4. Follow any additional prompts required by the script and inspect the results to ensure performance is as expected <br />

## Scripts
* `SSI_boundary_gen.js` 
  - Generates vertical and lateral boundary conditions for driving 3D SSI models with a compliant velocity field (appropriate for asymmetric or symmetric horizontal geometries)  <br />
* `interpolate_MAT.js` 
  - Assigns materials for elements in a 3D soil model based on soil unit type and nearest distance to materials in an aligned index soil column(s) model  <br />
* `symmetric_lateral_boundaries.js`
  - Creates lateral boundary conditions (not vertical) for a 3D soil domain with uniform layering in Primer <br />
* `asymmetric_lateral_boundaries.js`
  - Creates lateral boundary conditions (not vertical) for a 3D soil domain with potentially non-uniform layering in Primer <br />
* `create_lysmer_dampers.js` 
  - Creates a non-reflective bottom boundary condition with properly scaled LOAD_NODE_POINT cards for time history analysis of a 3D soil domain in Primer <br />
* `write_csv.js` 
  - Enables a csv file of your choice to be written from Primer or d3Plot <br />
* `setupSRA.js` 
  - Creates lateral boundary conditions for 1D SRA in Primer <br />
* `pile.js` 
  - Assists in the generation of embedded pile foundation components in Primer <br />
* `Beam-to-Shell.js` 
  - Converts parts composed of beam elements to a fine mesh of shell elements in Primer for more rigorous analysis in LS-DYNA <br />
* `staged_construction_viewer.js`
  - A tool for interactive visualization of construction stages defined in Primer <br />
