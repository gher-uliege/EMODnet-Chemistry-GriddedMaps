### File and variable naming

`DIVA` and `DIVAnd` have been adapted to accept file and variable names with white spaces. However it is recommended to avoid white spaces as well as other characters (brackets etc) since they can produce issues with ERDDAP, for example.

It is suggested to replace whitespace with underscores in:
- file names
- directory names
- variable names.
```bash
Sea water temperature → Sea_water_temperature
Northeast Atlantic Ocean → Northeast_Atlantic_Ocean
```