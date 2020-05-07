push!(LOAD_PATH,"../src/")

using Documenter, LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD

makedocs(
	sitename = "CSG.jl ",
	assets = ["assets/lar.css", "assets/logo.png"],
	format = Documenter.HTML(
		prettyurls = get(ENV, "CI", nothing) == "true"
	),
	pages = [
		"Home" => "index.md",
		"Tutorial" => [
			"Model" => "model.md",
			"Boolean" => "boolean.md",
			"Examples" => "examples.md"
		],
		"Computational pipeline" => [
			"Spatial indices" => "spatial_indices.md",
			"Planar arrangement" => "planar_arrangement.md",
			"Splitting" => "splitting.md",
			"Congruence" => "congruence.md",
			"Topological Gift Wrapping" => "tgw.md",
			"Spatial_arrangement" => "spatial_arrangement.md",
			"Shells" => "shells.md"
		],
		"Utils" => "utils.md"
	]
)
