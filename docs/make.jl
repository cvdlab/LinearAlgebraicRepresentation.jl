push!(LOAD_PATH,"../src/")

using Documenter, LinearAlgebraicRepresentation

makedocs(
	format = :html,
	sitename = "LinearAlgebraicRepresentation.jl",
	assets = ["assets/lar.css", "assets/logo.png"]
	pages = [
		"Home" => "index.md",
		"L.A.R. Intro" => "lar.md",
		"Interface" => "interface.md",
		"Arrangement" => "arrangement.md",
		"Parametric primitives" => [
			"Mapper" => "mapper.md",
			"Assemblies" => "struct.md"
		],
		"Grid generation" => [
			"Cuboidal grids" => "largrid.md",
			"Simplicial grids" => "simplexn.md"
		],
		"Domain integration" => "integr.md",
	]
)
