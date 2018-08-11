try
    Pkg.installed("Triangle")
catch
    Pkg.clone("https://github.com/cvdlab/Triangle.jl.git")
    Pkg.checkout("Triangle", "ver-0.1.0") # Julia 0.6
    Pkg.build("Triangle")
end
