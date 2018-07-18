push!(LOAD_PATH,"../src/")

using Documenter, LARLIB

makedocs(
	format = :html,
	sitename = "LARLIB.jl",
	pages = [
		"Home" => "index.md",
		"LAR" => "lar.md",
		"Arrangement" => "arrangement.md",
		"Mapper" => "mapper.md",
		"Assemblies" => "struct.md",
		"Glossary" => "glossary.md"
	]
)
