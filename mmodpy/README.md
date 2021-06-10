# mmodpy

**M**ott**M**acDonald **O**asys/LS-**D**yna **Py**thon Toolkit

This project hosts a collection of Python functions that support LS-DYNA modeling and interaction with Oasys pre/post-processors.

Developer: @kevin.stanton

## Installation

1. Install with pip <br />
`pip install git+https://github.com/mottmacdonaldglobal/mmodpy`
2. Import the library <br />
`import mmodpy`

## Usage
`mmodpy.[module].[function](args)` or `mmodpy.[function](args)`

## Modules
* [beams](https://github.com/mottmacdonaldglobal/mmodpy/tree/main/mmodpy/beams)
  - Automatic beam material, section, and part card generation from batch data specified in text or excel format <br />
* [port](https://github.com/mottmacdonaldglobal/mmodpy/tree/main/mmodpy/port)
  - Transfer data from various software packages to LS-DYNA <br />
* [bat](https://github.com/mottmacdonaldglobal/mmodpy/tree/main/mmodpy/bat)
  - General Windows batch automation scripts <br />
* [js](https://github.com/mottmacdonaldglobal/mmodpy/tree/main/mmodpy/js)
  - Collection of JavaScript files for automation in Oasys Primer or d3Plot <br />

## General Functions
* `run_TCF` 
  - Executes an Oasys T-HIS FastTCF script on the local machine <br />
* `find_replace` 
  - Walks through a directory and all subdirectories and replaces specific text within 
    files matching a user-defined file type extension <br />
* `run_dyna`
  - Executes an LS-DYNA analysis on the local machine <br />
* `lin_scale`
  - Linearly scales seed motions to match a specified target response spectra. The scale factor is determined according to a weighting scheme either specified by the user or assumed to be constant for all periods. <br />
* `th_to_key` 
  - Writes input velocity time histories keyword for LS-DYNA using acceleration time-history data (m/sec2) from a csv file <br />
