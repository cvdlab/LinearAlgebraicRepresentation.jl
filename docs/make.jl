push!(LOAD_PATH,"../src/")

using Documenter, LinearAlgebraicRepresentation

makedocs(
	format = :html,
	sitename = "LinearAlgebraicRepresentation.jl",
	pages = [
		"Home" => "index.md",
		"LAR" => "lar.md",
		"Interface" => "interface.md",
		"Arrangement" => "arrangement.md",
		"Mapper" => "mapper.md",
		"Assemblies" => "struct.md",
		"Cuboidal grids" => "largrid.md",
		"Simplicial grids" => "simplexn.md",
		"Domain integration" => "integr.md",
	]
)
