using DIVAnd
using Downloads
using Glob
using NCDatasets

datadir = "/home/ctroupin/data/EMODnet-Chemistry/Eutrophication2024/Results/Monthly"
outputdir = "/home/ctroupin/data/EMODnet-Chemistry/Eutrophication2024/XML"
datafilelist = Glob.glob("*.nc", datadir)
ignore_errors = true

const project = "EMODNET-chemistry"
cdilist = "./CDI-list-export.zip"

if !isfile(cdilist)
	Downloads.download("http://emodnet-chemistry.maris2.nl/download/export.zip", cdilist)
end

for datafile in datafilelist
	@info("Working on file $(basename(datafile))")

	ds = NCDataset(datafile, "r")
	varname = first(keys(ds))
	@info("Variable name: ")
	close(ds)

	xmlfilename = joinpath(outputdir, "Water_body_$(replace(varname," "=>"_")).xml")
	@info("Creating the XML file")
	divadoxml(datafile, varname, project, cdilist, xmlfilename, ignore_errors = ignore_errors)

end



