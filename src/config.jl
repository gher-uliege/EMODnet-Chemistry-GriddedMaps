deltalon = 0.25
deltalat = 0.25
lonr = -45.:deltalon:70.
latr = 24.:deltalat:83.


datadir = "/home/ctroupin/data/EMODnet-Chemistry/Eutrophication2024/netCDF"
figdir = "/home/ctroupin/Projects/EMODnet/EMODnet-Chemistry-GriddedMaps/figures/paper/"

varlist = [
    "Water_body_ammonium",
    "Water_body_chlorophyll-a",
    "Water_body_dissolved_inorganic_nitrogen",
    "Water_body_silicate",
    "Water_body_phosphate",
    "Water_body_dissolved_oxygen_concentration",
]

regionnames = [
    "Arctic Ocean",
    "Northeast Atlantic Ocean",
    "Baltic Sea",
    "Black Sea",
    "Mediterranean Sea",
    "North Sea",
]

colorlist = ["#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f"]
colorlist = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"]

function make_scatter(datafilelist::Array, varname::String)
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
            label = regionnames[iii],
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

function make_hexbin(datafilelist::Array, varname::String)
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

    obslonvec = Float64[]
    obslatvec = Float64[]

    for (iii, datafile) in enumerate(datafilelist)
        @debug("Reading data from $(basename(datafile))")
        obsvalue, obslon, obslat, obsdepth, obstime, obsids =
            loadobs(Float64, datafile, varname)

        append!(obslonvec, obslon)
        append!(obslatvec, obslat)
       
    end
    
    hexbin!(ga, obslonvec, obslatvec, bins = 50)
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
    fig = Figure(size = (1200, 400))
    ax = Axis(fig[1, 1], title = "Number of $(varname) observations per year", xticks=1960:10:2020)
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
        title = "Number of $(varname) observations per month",
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