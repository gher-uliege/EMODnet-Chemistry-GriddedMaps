## Computing weights "offline"

The weight computation is necessary when one wants to decrease the weight of observations according to their distance (in space and in time). This computation can be very expensive for large datasets, hence it is preferred to use perform it once, save the results in a file and read it whenever necessary.

### Example 

(With the obvious adaptation in the paths and file names)

```julia
using DIVAnd
using JLD
datadir = "/data/SeaDataCloud/NorthSea/"
varname = "Salinity"
obsfile = joinpath(datadir, "NorthSea_obs.nc")
netcdfODV = joinpath(datadir, "data_from_SDC_NS_DATA_DISCRETE_TS_V1b.nc")
@info("Reading data from the observation file")
@time obsval, obslon, obslat, obsdepth, obstime, obsid = DIVAnd.loadobs(Float64,obsfile,varname)
@info("Total number of data points: $(length(obsval))");

@time rdiag=1.0./DIVAnd.weight_RtimesOne((obslon,obslat),(0.03,0.03));
@show maximum(rdiag),mean(rdiag)
save("northsea_weights.jld", "rdiag", rdiag);
```