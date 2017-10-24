try
   Pkg.installed("TRIANGLE")
catch
   Pkg.clone("https://github.com/furio/TRIANGLE.jl.git")
   Pkg.build("TRIANGLE")
end
