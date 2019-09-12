using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
import  Base.union


function boolops(assembly::Lar.Struct, op::Symbol)
    W, copEV, copFE, boolmatrix = Lar.bool2d(assembly)
    boolvars = [boolmatrix[:,k] for k=1:size(boolmatrix,2)]
    if eval(op) == -
        solution = boolvars[1] .& .!(.|(boolvars[2:end]...))
    else
        operator = eval(op)
        solution = operator.(boolvars...)
    end
    chain1d = sparse(copFE') * Int.(solution)
    EV = Lar.cop2lar(copEV)
    EVop = [ev for (k,ev) in enumerate(EV) if abs(chain1d[k])==1 ]
    V = convert(Lar.Points,W')
    return V,EVop
end


wall = Lar.svg2lar("test/svg/assembly/wall.svg", flag=false)
openings = Lar.svg2lar("test/svg/assembly/openings.svg", flag=false)
rectangle = Lar.svg2lar("test/svg/assembly/rectangle.svg", flag=false)
box = Lar.svg2lar("test/svg/assembly/box.svg", flag=false)

assembly = Lar.Struct([wall, openings, rectangle])
V,EV = boolops(assembly, :|)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

assembly = Lar.Struct([wall, openings, rectangle])
V,EV = boolops(assembly, :-)
diff = Lar.Struct([ (V,EV) ])
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

assembly = Lar.Struct([wall, openings, box])
V,EV = boolops(assembly, :&)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

V,EV = boolops( Lar.Struct([ box,diff ]), :-)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
