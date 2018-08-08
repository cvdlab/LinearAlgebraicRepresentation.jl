try
   Pkg.installed("Triangle")
catch
   Pkg.clone("https://github.com/cvdlab/Triangle.jl.git")
   Pkg.build("Triangle")
end
