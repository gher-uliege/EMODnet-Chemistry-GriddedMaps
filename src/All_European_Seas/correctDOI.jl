
using NCDatasets
using Glob
datadir = "/home/ctroupin/data/EMODnet-Chemistry/Eutrophication2024/Results/Monthly"
datafilelist = Glob.glob("*.nc", datadir)
const doiprefix = "https://doi.org/10.13120/" #

for datafile in datafilelist[1:1]
    @info("Working on file $(basename(datafile))")

    NCDataset(datafile, "a") do ds
        doisuffix = ds.attrib["product_id"]

        newdoi = joinpath(doiprefix, doisuffix)
        @info("New DOI = $(newdoi)")
        ds.attrib["doi"] = newdoi
        ds.attrib["Conventions"] = "CF-1.12" ;
        ds.attrib["source"] = "Observational data from SeaDataNet/EMODnet Chemistry Data Network" ;
        ds.attrib["bathymetry_source"] = "The GEBCO 30sec Digital Atlas published by the British Oceanographic Data Centre on behalf of IOC and IHO, 2003" 
        ds.attrib["bathymetry_source"] = "GEBCO Compilation Group (2020) GEBCO 2020 Grid (doi:10.5285/a29c5465-b138-234d-e053-6c86abc040b9)"
        ds.attrib["Acknowledgements"] = "Aggregated data products are generated by EMODnet Chemistry under the support of DG MARE Call for Tenders EASME/EMFF/2016/006-lot4, EASME/2019/OP/0003-lot4." ;
        ds.attrib["project"] = "EMODnet Chemistry - Phase 5" ;
        ds.attrib["citation"] = "Usage is subject to mandatory citation: \"This resource was generated in the framework of EMODnet Chemistry, under the support of DG MARE Call for Tender EASME/EMFF/2020/3.1.11/European Marine Observation and Data Network (EMODnet) - Lot 5 - Chemistry\""
        ds.attrib["DIVAnd_source"] = "https://github.com/gher-uliege/DIVAnd.jl" ;
		ds.attrib["DIVAnd_version"] = "2.7.12" ;
		ds.attrib["DIVAnd_code_doi"] = "https://doi.org/10.5281/zenodo.1303229" ;
        ds.attrib["DIVAnd_reference"] = "Barth, A.; Beckers, J.-M.; Troupin, C.; Alvera-Azcarate, A. & Vandenbulcke, L. divand-1.0: n-dimensional variational data analysis for ocean observations. Geoscientific Model Development, 2014, 7, 225-241. DOI:10.5194/gmd-7-225-2014" ;

    end
end

##20_sharloSextant_25