# Port Module

This module facilitates LS-DYNA model generation from data originating from other programs.

<p>Developer: @Kevin.Stanton</p>


## Prerequisites

Required:
* Python 3
    
## Example Usage

1. Install `mmodpy` with pip <br />
`pip install git+https://github.com/mottmacdonaldglobal/mmodpy`
2. Import the `port` module <br />
`import mmodpy.port as port`
3. Run the `leapfrog_to_dyna` function <br />
`port.leapfrog_to_dyna(input_csv_path)`

## Workflow
![Diagram](https://github.com/mottmacdonaldglobal/mmodpy/blob/main/mmodpy/port/leapfrog_to_dyna.PNG)

## Options    

```python
def leapfrog_to_dyna(input_csv, output='model.key'):
    '''
    This function reads in data from LeapFrog and produces an LS-DYNA keyword
    file that contains PART, solid_data, and NODE data. Intended for
    porting 3D soil mesh data.
    
    Parameters
    ----------     
    input_csv : str
        Full file path for the input LeapFrog .csv data file
        Only filename required if located in cwd
    output : str
        Full file path desired for the output file (.key)
        Only filename required if located in cwd
        Default = 'model.key'
    
    Returns
    -------
    {file}.key : file
        LS-DYNA keyword file
    '''
```
