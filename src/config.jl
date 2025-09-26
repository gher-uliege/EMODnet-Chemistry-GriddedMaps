using GeoJSON
using PolygonOps
using OrderedCollections
using Dates
using ColorSchemes
using NCDatasets
using CairoMakie, GeoMakie
using GeoDatasets
using GeometryOps, GeoInterface

# Domain and resolution
deltalon = 0.25
deltalat = 0.25
lonr = -45.:deltalon:70.
latr = 24.:deltalat:83.

# Vertical levels
depthr = Float64[ 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

# Parameters for plotting
inch = 96
pt = 4 / 3
cm = inch / 2.54

# Files and directories
hostname = gethostname()
if hostname == "ctroupin-PRIME-A320M-C-R2-0"
    datadir = "/home/ctroupin/data/EMODnet-Chemistry/Eutrophication2024/netCDF"
else
    datadir = "/home/ctroupin/data/EMODnet/EMODnet-Chemistry/2024/netCDF"
end

figdir = "../figures/paper"
woadir = "/home/ctroupin/data/WOA/"

# Land/sea mask and coastline
lon_landsea, lat_landsea, landsea = GeoDatasets.landseamask(; resolution = 'i', grid = 5)
landsea[landsea.==2] .= 1;
landsea = Float64.(landsea)
landsea[landsea.==0] .= -999;

# Subset land/sea mask
goodlon = findall( (lon_landsea .>= first(lonr)) .& (lon_landsea .<= last(lonr)))
goodlat = findall( (lat_landsea .>= first(latr)) .& (lat_landsea .<= last(latr)))
lon_landsea = lon_landsea[goodlon]
lat_landsea = lat_landsea[goodlat]
landsea = landsea[goodlon, goodlat]

coordscoast = GeoDatasets.gshhg("i", 1);

domains = OrderedDict(
    "Arctic Ocean" => [-42, 54, 62, 83.5],
    "Northeast Atlantic Ocean" => [-42., -0.1, 25., 47.9],
    "Baltic Sea" => [9.4, 30.9, 53, 65.39],
    "Black Sea" => [26.5, 42., 40., 48.],
    "Mediterranean Sea" => [-7, 36.5, 30, 46],
    "North Sea" => [-20., 25., 47., 63.]
)

domainrivers = OrderedDict(
    "Loire River" => [-4., -1., 46.25, 48.],
    "Gulf of Riga" => [22.3, 25.0, 56.8, 58.4],
    "Po River" => [12., 14., 44., 46.],
    "Danube Delta" => [28.5, 30.5, 43.7, 45.6]
)

varlist = [
    "Water_body_ammonium",
    "Water_body_chlorophyll-a",
    "Water_body_dissolved_inorganic_nitrogen",
    "Water_body_silicate",
    "Water_body_phosphate",
    "Water_body_dissolved_oxygen_concentration",
]

varunits = Dict(
    "Water_body_ammonium" => "µmol/l",
    "Water_body_chlorophyll-a" => "mg/m{^3}",
    "Water_body_dissolved_inorganic_nitrogen" => "µmol/l",
    "Water_body_silicate" => "µmol/l",
    "Water_body_phosphate" => "µmol/l",
    "Water_body_dissolved_oxygen_concentration" => "µmol/l"
    )

# Acceptable range for the variables
# (used for the histograms)
varrange = Dict(
    "Water_body_ammonium" => [0., 30.],
    "Water_body_chlorophyll-a" => [0., 5.],
    "Water_body_dissolved_inorganic_nitrogen" => [0., 50.],
    "Water_body_silicate" => [0., 100.],
    "Water_body_phosphate" => [0., 5.],
    "Water_body_dissolved_oxygen_concentration" => [0., 400.]
)

regionnames = [
    "Arctic Ocean",
    "Northeast Atlantic Ocean",
    "Baltic Sea",
    "Black Sea",
    "Mediterranean Sea",
    "North Sea",
]

colorlist = ["#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f"]


domaincolors = OrderedDict(
    "Arctic Ocean" => "#e41a1c",
    "Northeast Atlantic Ocean" => "#377eb8",
    "Baltic Sea" => "#4daf4a",
    "Black Sea" => "#984ea3",
    "Mediterranean Sea" => "#ff7f00",
    "North Sea" => "#ffff33"
)

domainriverscolors = OrderedDict(
    "Loire River" => domaincolors["Northeast Atlantic Ocean"],
    "Gulf of Riga" => domaincolors["Baltic Sea"],
    "Po River" => domaincolors["Mediterranean Sea"],
    "Danube Delta" => domaincolors["Black Sea"]
)

"""
    get_unique_coordinates(filelist)

Read the coordinates from a list of netCDF files and keep the unique pairs of coordinates.
"""
function get_unique_coordinates(filelist::Vector{String}, varname::AbstractString)

    obslonvec = Float64[]
    obslatvec = Float64[]

    for (iii, datafile) in enumerate(filelist)
        @debug("Reading data from $(basename(datafile))")
        obsvalue, obslon, obslat, obsdepth, obstime, obsids =
            loadobs(Float64, datafile, varname)

        append!(obslonvec, obslon)
        append!(obslatvec, obslat)
       
    end
    coordinates = tuple.(obslonvec, obslatvec)
    unique!(coordinates)
    lon_u = first.(coordinates)
    lat_u = last.(coordinates)
    return lon_u::Array{Float64}, lat_u::Array{Float64}
end

"""
    make_scatter(datafilelist, varname)

Plot the locations of the `varname` observations by region using the list of file `datafilelist` 
"""
function make_scatter(datafilelist::Array, varname::String)

    local colorlist = collect(values(domaincolors))
    local domainlist = collect(keys(domaincolors))

    varname_ = replace(varname, "_" => " ")
    fig = Figure(size = (600, 600))
    ga = GeoAxis(
        fig[1, 1],
        title = "Observations of $(varname_)",
        dest = "+proj=ortho +lon_0=15 +lat_0=35",
    )
    heatmap!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea,
        colormap = Reverse(:greys),
        colorrange = [0, 2],
    )
    xlims!(-180, 180.0)
    ylims!(-90.0, 90.0)
    for (iii, datafile) in enumerate(datafilelist)
        @debug("Reading data from $(basename(datafile))")
        @time obsvalue, obslon, obslat, obsdepth, obstime, obsids =
            loadobs(Float64, datafile, varname)

        plt1 = plot!(
            ga,
            [NaN],
            [NaN],
            markersize = 5,
            color = colorlist[iii],
            label = domainlist[iii]
        )
        plot!(ga, obslon, obslat, markersize = 2, color = colorlist[iii])
    end
    leg = axislegend(
        ga,
        "Regions",
        framevisible = true,
        framecolor = :grey,
        position = :lb,
        backgroundcolor = :white,
        alpha = 0.85,
    )
    hidedecorations!(ga)
    return fig
end

"""
    make_hexbin(obslon, obslat, varname)

Create a hexbin plot using the coordinate vectors `obslon`, `obslat`.     
The variable `varname` is used for the figure title.
"""
function make_hexbin(obslonvec::Array{Float64}, obslatvec::Array{Float64}, varname::String="")
    
    if length(varname) > 0
        varname_ = replace(varname, "_" => " ")
        figtitle = "$(varname_)"
    else
        figtitle = ""
    end
    fig = Figure(size = (800, 600))
    ga = GeoAxis(
        fig[1, 1],
        title = figtitle,
        dest = "+proj=ortho +lon_0=15 +lat_0=35",
    )
    
    xlims!(-180, 180.0)
    ylims!(-90.0, 90.0)
    
    hb = hexbin!(ga, obslonvec, obslatvec, bins = 75, colormap = Reverse(:RdYlBu), strokewidth = 0.5,
        strokecolor = :gray50, threshold = 1, colorscale = log10)
    
    Colorbar(fig[1,2], hb, vertical = true, label = "Number \nof profiles\nper cell", labelrotation=0)
    heatmap!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea,
        colormap = Reverse(:greys),
        colorrange = [0, 2],
    )
    hidedecorations!(ga)
    return fig
end


function make_histogram(
    yearmin::Int64,
    yearmax::Int64,
    nobs_year::Vector{Int64},
    nobs_month::Vector{Int64},
    varname::String
)
    varname_ = replace(varname, "_" => " ")
    fig = Figure(size = (1200, 400))
    ax = Axis(fig[1, 1], title = "Number of $(varname_) observations per year", xticks=1960:10:2020)
    barplot!(
        ax,
        yearmin:yearmax,
        nobs_year,
        strokewidth = 0.5,
        strokecolor = (:black, 0.5),
        color = :grey,
    )
    xlims!(ax, 1960 - 0.5, 2024 + 0.5)
    hidespines!(ax, :t, :r)

    ax2 = Axis(
        fig[1, 2],
        xticks = 1:12,
        xtickformat = x -> Dates.monthname.(Int.(x)),
        xticklabelrotation = pi / 6,
        title = "Number of $(varname_) observations per month",
    )
    barplot!(
        ax2,
        1:12,
        nobs_month,
        strokewidth = 0.5,
        strokecolor = (:black, 0.5),
        color = :grey,
    )
    hidespines!(ax2, :t, :r)
    return fig
end


"""
    read_polygon_json(contourfile)

Read the coordinates as a list of tuples stored in the geoJSON file `contourfile`,
as downloaded from https://geojson.io

# Example
```julia-repl
julia> coordlist = read_polygon_json(contourfile)
```
"""
function read_polygon_json(contourfile::AbstractString)
    coordlist = []
    jsonbytes = read(contourfile);
    fc = GeoJSON.read(jsonbytes)
    for poly in fc
        coordinates = poly.geometry[1]
        push!(coordlist, coordinates)
    end
    return coordlist
end

"""
    edit_mask!(xi, yi, mask, coordinatelist)

Edit the land-sea mask (as read using `DIVAnd.load_bath`) by setting to zero (_land_ value)
the cells that fall inside the polygon(s) defined by `coordinatelist` (as read by using `read_polygon_json`).

The contour file can be downloaded from geojson.io. An example of such a file (`outsidemask.json`) is provided in 
the data directory.

# Example
```julia-repl
julia> coordinatelist = read_polygon_json(contourfile);
julia> edit_mask!(xi, yi, mask, coordinatelist)
```
"""
function edit_mask!(xi, yi, mask::BitMatrix, coordinatelist::Vector{Any})
    for (ii, xx) in enumerate(xi)
        for (jj, yy) in enumerate(yi)
            for coordinates in coordinatelist
                if PolygonOps.inpolygon((xx,yy), coordinates) != 0
                    mask[ii, jj] = 0.
                end
            end
        end
    end
    return nothing
end

"""
    draw_domain(ga, coords, thecolor, name, linestyle)

Draw the domain on the map.

# Example

"""
function draw_domain(ga::GeoAxis, coords::Vector{Float64}, thecolor, name::String="", linestyle=:solid)
    poly = GeoInterface.LineString([
        (coords[1], coords[3]),
        (coords[2], coords[3]),
        (coords[2], coords[4]),
        (coords[1], coords[4]),
        (coords[1], coords[3]),
    ])
    spoly = GeometryOps.segmentize(poly, max_distance = 1)
    lines!(ga, GeoMakie.geo2basic(spoly); linestyle = linestyle, color = thecolor, linewidth = 2, label = name)
    return nothing
end

"""
    read_results(resultfile, varname, depth2plot, month2plot)

Read the results obtained with `DIVAnd`
colors
# Example
```julia-repl
julia> lon, lat, depth, dates, field, error = read_results(resultfile, "Water_body_phosphate", 20, 4);
```
"""
function read_results(
    resultfile::String,
    varname::String,
    depth2plot::Int64,
    month2plot::Int64,
)
    NCDataset(resultfile) do ds
        lon = ds["lon"][:]
        lat = ds["lat"][:]
        depth = ds["depth"][:]
        dates = ds["time"][:]

        # Subsetting
        depthindex = findfirst(depth .== depth2plot)
        timeindex = findfirst(Dates.month.(dates) .== month2plot)
        # Read gridded field and error
        field = coalesce.(ds[varname][:, :, depthindex, timeindex], NaN)
        error = coalesce.(ds[varname*"_relerr"][:, :, depthindex, timeindex], NaN)

        return lon::Vector{Float64},
        lat::Vector{Float64},
        depth::Vector{Float64},
        dates::Vector{DateTime},
        field::Matrix{AbstractFloat},
        error::Matrix{AbstractFloat}
    end
end

"""
    read_woa(woafile, varname, depth2plot, lonr, latr)

Read data from the World Ocean Atlas

# Example
```julia-repl
julia> lonwoa, latwoa, depthwoa, dateswoa, fieldwoa, errorwoa =
    read_woa("woa23_all_i10_01.nc", "Water_body_phosphate", 20, lonr, latr)
```
"""
function read_woa(woafile::String, varname::String, depth2plot::Int64, lonr::StepRangeLen, latr::StepRangeLen)
    NCDataset(woafile, "r") do ds
        lon = ds["lon"][:]
        lat = ds["lat"][:]
        depth = ds["depth"][:]
        dates = ds["time"][:]

        # Subsetting
        goodlon = findall((lon .<= lonr[end]) .& (lon .>= lonr[1]))
        goodlat = findall((lat .<= latr[end]) .& (lat .>= latr[1]))
        depthindex = findfirst(depth .== depth2plot)
        timeindex = 1

        if varname == "Water_body_phosphate"
            varnamenc = "p_an"
            errnamenc = "p_sea"
        elseif varname == "Water_body_silicate"
            varnamenc = "i_an"
            errnamenc = "i_sea"
        elseif varname == "Water_body_dissolved_oxygen_concentration"
            varnamenc = "o_an"
            errnamenc = "o_sea"
        end

        # Read gridded field and error
        field = coalesce.(ds[varnamenc][goodlon, goodlat, depthindex, 1], NaN)
        error = coalesce.(ds[errnamenc][goodlon, goodlat, depthindex, 1], NaN)

        return lon[goodlon]::Vector{Float32},
        lat[goodlat]::Vector{Float32},
        depth::Vector{Float32},
        dates::Vector{DateTime},
        field::Matrix{AbstractFloat},
        error::Matrix{AbstractFloat}
    end
end

"""
    add_coast!(ga, coordscoast)

Add the coast read from GSHHG to the plot

# Example
```julia-repl
julia> coordscoast = GeoDatasets.gshhg("i", 1)
julia> fig = Figure(); ga = GeoAxis(fig[1,1])
julia> add_coast!(ga, coordscoast)
```
"""
function add_coast!(ga::GeoAxis, coordscoast::Vector{Tuple{Vector{Float64}, Vector{Float64}}})
    for iii = 1:length(coordscoast)
        lonc = coordscoast[iii][1]
        latc = coordscoast[iii][2]
        lonc[lonc.>=lonr[end]] .= NaN
        lonc[lonc.<=lonr[1]] .= NaN
        latc[latc.>=latr[end]] .= NaN
        latc[latc.<=latr[1]] .= NaN
        lines!(ga, lonc, latc, color = :black, linewidth = 0.5)
    end
    return nothing
end

"""
    plot_field_var(varname, lon, lat, field, depth2plot, month2plot)

Plot the 2D field corresponding to the coordinates `lon` and `lat`, 
at the depth level `depth2plot` and the time period `month2plot`.
"""
function plot_field_var(
    varname::String,
    lon,
    lat,
    field,
    depth2plot,
    month2plot,
    source::String = "DIVAnd",
    cmap = cgrad(:RdYlBu, rev = true);
    vmin::Float64=0.,
    vmax::Float64=1.
)
    varname_ = replace(varname, "_" => " ")

    fig = Figure()
    ga = GeoAxis(
        fig[1, 1],
        title = "$(source) $(varname_) field\nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    heatmap!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea,
        colormap = Reverse(:greys),
        colorrange = [0, 2],
    )

    hm = heatmap!(
        ga,
        lon,
        lat,
        field,
        colorrange = (vmin, vmax),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    # add_coast!(ga, coordscoast)
    contour!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea2,
        levels = [-0.1, 0.0],
        color = :black,
        linewidth = 0.5,
    )

    xlims!(ga, lonr[1], lonr[end])
    ylims!(ga, latr[1], latr[end])
    Colorbar(fig[1, 2], hm, label = varunits[varname], labelrotation = 0)
    return fig, ga, hm
end


"""
    plot_field_var_deepest(varname, lon, lat, field, depth2plot, month2plot)

Plot the 2D field corresponding to the coordinates `lon` and `lat`, 
at the depth level `depth2plot` and the time period `month2plot`.
"""
function plot_field_var_deepest(
    varname::String,
    lon,
    lat,
    field,
    month2plot,
    source::String = "DIVAnd",
    cmap = cgrad(:RdYlBu, rev = true);
    vmin::Float64=0.,
    vmax::Float64=1.
)
    varname_ = replace(varname, "_" => " ")

    fig = Figure()
    ga = GeoAxis(
        fig[1, 1],
        title = "$(source) $(varname_) field\nat deepest depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    heatmap!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea,
        colormap = Reverse(:greys),
        colorrange = [0, 2],
    )

    hm = heatmap!(
        ga,
        lon,
        lat,
        field,
        colorrange = (vmin, vmax),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    # add_coast!(ga, coordscoast)
    contour!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea2,
        levels = [-0.1, 0.0],
        color = :black,
        linewidth = 0.5,
    )

    xlims!(ga, lonr[1], lonr[end])
    ylims!(ga, latr[1], latr[end])
    Colorbar(fig[1, 2], hm, label = varunits[varname], labelrotation = 0)
    return fig, ga, hm
end

function plot_error_field(varname, lon, lat, fielderror, depth2plot, month2plot, source::String = "DIVAnd",
    cmaperror = cgrad(:RdYlGn_10, 10, rev = true, categorical = true);
    vmin::Float64=0.,
    vmax::Float64=1.
    )

    varname_ = replace(varname, "_" => " ")
    fig = Figure()

    ga = GeoAxis(
        fig[1, 1],
        title = "$(source) $(varname_) relative error\nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm = heatmap!(
        ga,
        lon,
        lat,
        fielderror * 100.0,
        colorrange = (0, 100.0),
        colormap = cmaperror
    )

    #add_coast!(ga, coordscoast)
    contour!(
        ga,
        lon_landsea,
        lat_landsea,
        landsea2,
        levels = [-0.1, 0.0],
        color = :black,
        linewidth = 0.5,
    )

    xlims!(ga, lonr[1], lonr[end])
    ylims!(ga, latr[1], latr[end])
    Colorbar(fig[1, 2], hm, label = "%", labelrotation = 0)
    return fig, ga, hm
end

"""
    plot_additional_field(varname, varunits, lon, lat, fielddeepest, fieldL2, fielderror, depth2plot, month2plot)

Plot the additional fields: masked gridded field, error field, field at deepest depth.
"""
function plot_additional_field(
    varname::String,
    varunits::String,
    lon,
    lat,
    fielddeepest, fieldL1, fieldL2, fielderror,
    depth2plot,
    month2plot,
    cmap = cgrad(:RdYlBu, rev = true),
    cmaperror = cgrad(:RdYlGn_10, 10, rev = true, categorical = true)
)
    varname_ = replace(varname, "_" => " ")

    fig = Figure(size=(1000, 650))

    ga1 = GeoAxis(
        fig[1, 1],
        title = "(a) $(varname_) masked\nusing a 50% relative error threshold\nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm1 = heatmap!(
        ga1,
        lon,
        lat,
        fieldL2,
        colorrange = (0, 1.0),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    add_coast!(ga1, coordscoast)

    xlims!(ga1, lonr[1], lonr[end])
    ylims!(ga1, latr[1], latr[end])
    Colorbar(fig[1, 2], hm1, label = varunits, labelrotation = 0)

    ga2 = GeoAxis(
        fig[1, 3],
        title = "(b) $(varname_) masked\nusing a 30% relative error threshold\nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm2 = heatmap!(
        ga2,
        lon,
        lat,
        fieldL1,
        colorrange = (0, 1.0),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    add_coast!(ga2, coordscoast)

    xlims!(ga2, lonr[1], lonr[end])
    ylims!(ga2, latr[1], latr[end])
    Colorbar(fig[1, 4], hm2, label = varunits, labelrotation = 0)

    ga3 = GeoAxis(
        fig[2, 1],
        title = "(c) Deepest value of $(varname_)\nin $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm3 = heatmap!(
        ga3,
        lon,
        lat,
        fielddeepest,
        colorrange = (0, 5.0),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    add_coast!(ga1, coordscoast)

    xlims!(ga3, lonr[1], lonr[end])
    ylims!(ga3, latr[1], latr[end])
    Colorbar(fig[2, 2], hm1, label = varunits, labelrotation = 0)
   

    ga4 = GeoAxis(
        fig[2, 3],
        title = "(d) $(varname_) relative error\nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm4 = heatmap!(
        ga4,
        lon,
        lat,
        fielderror * 100.0,
        colorrange = (0, 100.0),
        colormap = cmaperror,
        highclip = cmap.colors[end],
    )

    add_coast!(ga4, coordscoast)

    xlims!(ga4, lonr[1], lonr[end])
    ylims!(ga4, latr[1], latr[end])
    Colorbar(fig[2, 4], hm4, label = "%", labelrotation = 0)

    resize_to_layout!(fig)
    return fig

end



function plot_field_var_fast(
    varname::String,
    lon,
    lat,
    field,
    depth2plot,
    month2plot,
    source::String = "DIVAnd",
    cmap = cgrad(:RdYlBu, rev = true),
)
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1],
        title = "$(source) $(varname) field \nat $(Int64(depth2plot)) m depth in $(Dates.monthname(month2plot))\n\n",
        dest = "+proj=laea +lon_0=15 +lat_0=45",
        xticks = (-50:20.0:70),
        yticks = (20:10.0:85),
    )

    hm = heatmap!(
        ga,
        lon,
        lat,
        field,
        colorrange = (0., 1.0),
        colormap = cmap,
        highclip = cmap.colors[end],
    )

    xlims!(ga, lonr[1], lonr[end])
    ylims!(ga, latr[1], latr[end])
    Colorbar(fig[1, 2], hm)
    return fig, ga, hm
end