try
    checkTriangle = Pkg.installed("Triangle")
    if checkTriangle == nothing
        Pkg.clone("https://github.com/cvdlab/Triangle.jl.git")
        Pkg.checkout("Triangle", "ver-0.1.0") # Julia 0.6
        Pkg.build("Triangle")
    end
catch
    error("Cannot build dependencies")
end
