## Duplicate removal

The duplicate detection and removal is a tedious task and might require to look individually to some profiles. We propose here a procedure that can help users to identify most of the duplicates.

1. Use `ODV` tool to remove the duplicates based on the metadata.
2. With `DIVAnd`, run the duplicate check with the function [`checkduplicates`](https://gher-uliege.github.io/DIVAnd.jl/stable/#DIVAnd.Quadtrees.checkduplicates).
```julia
dupl = DIVAnd.Quadtrees.checkduplicates(
    (obslon1,obslat1, obsdepth1, obstime1), obsval1,
    (obslon2,obslat2, obsdepth2, obstime1), obsval2,
    (Δlon, Δlat, Δdepth, Δtime), Δval);
```
with *reasonable* values for the Δ parameters (values below which two data points are considered as duplicates).      
For instance:
- Δlon = Δlat = 0.01°
- Δdepth = 2 meters
- Δtime = 1 hour (units should be the same as `obstime`)
- Δval = 0.1°C (if you work with temperature).
3. Check the percentage of supposed duplicates with respect to the 1st dataset
```julia
ndupl = length(findall(.!isempty.(dupl)));
pcdupl = round(ndupl / length(obslon) * 100; digits=2);
```
4. Return to step 2 and modify the Δ parameters to see how they affect the percentage of duplicates.
5. Once you are confident with the results, combine the 2 (or more)
data sets, as explained in the notebook [90-full-analysis](https://github.com/gher-uliege/Diva-Workshops/blob/master/notebooks/3-Analysis/90-full-analysis.ipynb).

  