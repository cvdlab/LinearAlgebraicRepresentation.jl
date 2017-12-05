try
   Pkg.installed("TRIANGLE")
catch
   Pkg.clone("https://github.com/cvdlab/TRIANGLE.jl.git")
   Pkg.build("TRIANGLE")
end
