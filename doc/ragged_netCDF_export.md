## Convert ODV netCDF to ragged netCDF

For large applications, the time needed to read the input file (either ODV spreadsheet or ODV netCDF) can be rather long. It is possible to speed up the process by re-writing the observations as vectors (instead of arrays of profiles).

### Procedure

1. Read the original netCDF ODV file
```julia
obsval, obslon, obslat, obsdepth, obstime, obsid =
    NCODV.load(Float64, ODVfile1, "Water body salinity");
```
2. Re-write the data
```julia
DIVAnd.saveobs(obsfile, "Water body salinity", obsval,
   (obslon, obslat, obsdepth, obstime), obsid)
```
3. Use the newly written files for the climatologies
```julia
obsval1, obslon1, obslat1, obsdepth1, obstime1, obsid1 =
   DIVAnd.loadobs(Float64,obsfile,"Water body salinity");
```