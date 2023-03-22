# Extended calibration of olivine-spinel Al exchange thermometers (OSAT)

Python scripts for the new models of the olivine-spinel Al exchange thermometers (OSAT). 

## Installation
```pip install https://github.com/eazzzon/olspthermo.git```

## Requirement
`numpy`

`pandas`

`scipy`

`periodictable`

## Uninstallation
```pip uninstall olspthermo```

## Example

To quick run, navigate to the `example` folder, call the python script as:

```python
python example_calculation.py
```

or run through the jupiter file.

**DON'T** change the header/column name of the example excel template files. 

## Monte Carlo simulation to progagation errors

We also provide a MC approach to progagate the analytical errors to the temperature estimation. an example is given by calling 

```
python example_mc.py
```

or run in the jupyter example.

## How to cite

Zhang Y, Namur O, Li W, Shorttle O, Gazel E, Jennings ES, Thy P, Grove TL, Charlier B (under review). An extended calibration of the olivine-spinel aluminum exchange thermometer: Application to the melting conditions and mantle lithology of large igneous provinces. 