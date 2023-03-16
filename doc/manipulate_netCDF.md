## Manipulating netCDF files

There are plenty of ways to edit a netCDF file as the one obtained with `DIVAnd`.      
One frequently has to edit the attributes, modify the variable names etc.

### Tools
The most useful is probably [`nco`](http://nco.sourceforge.net/) (netCDF operator).           

### Mask an area (Credits: Ö. Bäck [SMHI])

The goal is to mask (i.e. set to land) an area after the DIVA analysis.      
The area to mask is specified by the indices in the coordinates.      
Note the quotes around the variable name (which contains spaces).
```bash
ncap2 -s "'ITS-90 water temperature'(:, :, 112:224, 0:15)=9.96921e36;" infile.nc outfile.nc
```
numbers are indices of the matrix that consist of dimensions (time, depth, lat, lon), so the above masks an area specified by `lat[112:224]` and `lon[0:15]` at all depths and time steps.

### Cut the netCDF
With this command, one can cut a piece of the field according to the coordinates.

```bash
ncks -d lon,10.,31. infile.nc outfile.nc
```
numbers here specify min and max of what you want, so anything outside 10 and 31 degrees are removed from all fields in the netCDF and also in the lon variable.
