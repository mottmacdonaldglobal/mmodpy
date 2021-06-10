## Batch Script Module

This module is a collection of Windows .bat script to help with common tedious file manipulation tasks.

## Usage
1. Install `mmodpy` library with pip (or download desired files from GitHub) <br />
`pip install mmodpy`
2. Copy and paste the desired script into the appropriate directory <br />
3. Modify script parameters in a text editor as needed <br />
4. Double click to run the script <br />

## Scripts
* `cleanDirectory.bat` 
  - Removes all file in a given folder except for user-defined file-types (useful for prepping LS-DYNA analysis directories for new runs) <br />
* `dynaBatch.bat` 
  - Executes LS-DYNA analysis <br />
* `open_message_files.bat`
  - Walks through all files of the current directory and all subdirectories and opens `message` files output from LS-DYNA in Notepad++ <br />
* `thfBatch.bat` 
  - Executes Oasys T-HIS post-processes defined in a `FastTCF.inp` automation file (must be stored in the same directory) <br />
