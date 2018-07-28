"""Module for facet extraction, extrusion and simplicial grids"""

"""
	extrudeSimplicial(model, pattern)
	
Algorithm for multimensional extrusion of a simplicial complex. A full documentation 
may be found in the ``Simple_X^n`` module, coded as `simplexn`. `model` is a LAR model, i.e. a pair (vertices,cells) to be extruded, whereas pattern is an array of `Int64`, to be used as lateral measures of the *extruded* model. pattern elements are assumed as either *solid* or *empty* measures, according to (+/-) sign.

# Example
```julia

```
"""
function extrudeSimplicial(model, pattern)
    V, FV = model
    d, m = length(FV[1]), length(pattern)
    coords = collect(cumsum(append!([0], abs.(pattern))))
    offset, outcells, rangelimit = length(V), [], d*m
    for cell in FV  
        tube = [v+k*offset for k in range(0, m+1) for v in cell]
        cellTube = [tube[k:k+d] for k in range(1, rangelimit)]
        outcells = vcat(outcells, reshape(cellTube, d, m))
    end
    cellGroups = []
    for i in 1:size(outcells, 2)
        if pattern[i]>0
            cellGroups = vcat(cellGroups, outcells[:, i])
        end
    end
    outVertices = [vcat(v, [z]) for z in coords for v in V]
    cellGroups = convert(Array{Array{Int64, 1}, 1}, cellGroups)  ## aggiunta mia
    outModel = outVertices, cellGroups
end

"""
	simplexGrid(shape)
	
Generate a simplicial complex decomposition of a cubical grid of ``d``-cuboids, where ``d`` is the length of `shape` array. Vertices (0-cells) of the grid have `Int64` coordinates.

# Examples
```julia

julia> simplexGrid([0]) # 0-dimensional complex
# output
([0], Array{Int64,1}[])

julia> V,EV = simplexGrid([1]) # 1-dimensional complex
# output
([0 1], Array{Int64,1}[[1, 2]])

julia> V,FV = simplexGrid([1,1]) # 2-dimensional complex
# output
([0 1 0 1; 0 0 1 1], Array{Int64,1}[[1, 2, 3], [2, 3, 4]])

julia> V,CV = simplexGrid([10,10,1]) # 3-dimensional complex
# output
([0 1 … 9 10; 0 0 … 10 10; 0 0 … 1 1], Array{Int64,1}[[1, 2, 12, 122], [2, 12, 122, 123], [12, 122, 123, 133], [2, 12, 13, 123], [12, 13, 123, 133], [13, 123, 133, 134], [2, 3, 13, 123], [3, 13, 123, 124], [13, 123, 124, 134], [3, 13, 14, 124]  …  [119, 229, 230, 240], [109, 119, 120, 230], [119, 120, 230, 240], [120, 230, 240, 241], [109, 110, 120, 230], [110, 120, 230, 231], [120, 230, 231, 241], [110, 120, 121, 231], [120, 121, 231, 241], [121, 231, 241, 242]])

julia> V
# output
3×242 Array{Int64,2}:
 0  1  2  3  4  5  6  7  8  9  10  0  1  2  3  …   1   2   3   4   5   6   7   8   9  10
 0  0  0  0  0  0  0  0  0  0   0  1  1  1  1     10  10  10  10  10  10  10  10  10  10
 0  0  0  0  0  0  0  0  0  0   0  0  0  0  0      1   1   1   1   1   1   1   1   1   1

julia> using LARVIEW

julia> hpc = LARVIEW.lar2exploded_hpc(V,CV) # exploded visualization of the grid

julia> view(hpc)

julia> V,HV = simplexGrid([1,1,1,1]) # 4-dim cellular complex from the 4D simplex
# output
([0 1 … 0 1; 0 0 … 1 1; 0 0 … 1 1; 0 0 … 1 1], Array{Int64,1}[[1, 2, 3, 5, 9], [2, 3, 5, 9, 10], [3, 5, 9, 10, 11], [5, 9, 10, 11, 13], [2, 3, 5, 6, 10], [3, 5, 6, 10, 11], [5, 6, 10, 11, 13], [6, 10, 11, 13, 14], [3, 5, 6, 7, 11], [5, 6, 7, 11, 13]  …  [4, 6, 10, 11, 12], [6, 10, 11, 12, 14], [3, 4, 6, 7, 11], [4, 6, 7, 11, 12], [6, 7, 11, 12, 14], [7, 11, 12, 14, 15], [4, 6, 7, 8, 12], [6, 7, 8, 12, 14], [7, 8, 12, 14, 15], [8, 12, 14, 15, 16]])
```
"""
function simplexGrid(shape)
    model = [[]], [[1]]
    for item in shape
        model = extrudeSimplicial(model, fill(1, item))
    end
    V, CV = model
    V = hcat(V...)
    V = convert(Array{Float64,2}, V)
    return V, CV
end


VOID = V0,CV0 = [[]],[[0]]    # the empty simplicial model

function cumsum(iterable)
    # cumulative addition list(cumsum(range(4))) => [0, 1, 3, 6]
    iterable = iter(iterable)
    s = iterable.next()
    yield s
    for c in iterable
        s = s + c
        yield s
    end
end

function larExtrude1(model,pattern)
    V, FV = model
    d, m = len(FV[0]), len(pattern)
    coords = list(cumsum([0]+(AA(ABS)(pattern))))
    offset, outcells, rangelimit = len(V), [], d*m
    for cell in FV
        tube = [v + k*offset for k in range(m+1) for v in cell]  
        cellTube = [tube[kk+d+1] for k in range(rangelimit)]
        outcells += [reshape(cellTube, newshape=(m,d,d+1)).tolist()]      
    end
        
    outcells = AA(CAT)(TRANS(outcells))
    cellGroups = [group for k,group in enumerate(outcells) if pattern[k]>0 ]
    outVertices = [v+[z] for z in coords for v in V]
    outModel = outVertices, CAT(cellGroups)
    return outModel
end

function larSimplexGrid1(shape)
    model = VOID
    for item in shape
        model = larExtrude1(model,item*[1])
    end
    return model
end

function larSimplexFacets(simplices)
    out = []
    d = len(simplices[0])
    for simplex in simplices
        out += AA(sorted)([simplex[0k]+simplex[k+1d] for k in range(d)])
    end
    out = set(AA(tuple)(out))
    return  sorted(out)
end

""" Transformation to triangles by sorting circularly the vertices of faces """
function quads2tria(model)
   V,FV = model
   out = []
   nverts = len(V)-1
   for face in FV
      centroid = CCOMB([V[v] for v in face])
      V += [centroid] 
      nverts += 1
      
      v1, v2 = DIFF([V[face[0]],centroid]), DIFF([V[face[1]],centroid])
      v3 = VECTPROD([v1,v2])
      if ABS(VECTNORM(v3)) < 10**3
         v1, v2 = DIFF([V[face[0]],centroid]), DIFF([V[face[2]],centroid])
         v3 = VECTPROD([v1,v2])
      transf = mat(INV([v1,v2,v3]))
      verts = [(V[v]*transf).tolist()[0][-1]  for v in face]

      tcentroid = CCOMB(verts)
      tverts = [DIFF([v,tcentroid]) for v in verts]   
      rverts = sorted([[ATAN2(vert),v] for vert,v in zip(tverts,face)])
      ord = [pair[1] for pair in rverts]
      ord = ord + [ord[0]]
      edges = [[n,ord[k+1]] for k,n in enumerate(ord[-1])]
      triangles = [[nverts] + edge for edge in edges]
      out += triangles
   return V,out
end

if __name__ == "__main__"
   # example 1
   V = [[0,0],[1,0],[2,0],[0,1],[1,1],[2,1],[0,2],[1,2],[2,2]]
   FV = [[0,1,3],[1,2,4],[2,4,5],[3,4,6],[4,6,7],[5,7,8]]
   model = larExtrude1((V,FV),4*[1,2,-3])
   VIEW(EXPLODE(1,1,1.2)(MKPOLS(model)))
   
   # example 2
   model = larExtrude1( VOID, 10*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude1( model, 10*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude1( model, 10*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   
   # example 3
   model = larExtrude1( VOID, 10*[1,-1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   model = larExtrude1( model, 10*[1] )
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(model)))
   
   grid_2d = larSimplexGrid1([3,3])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_2d)))
   
   grid_3d = larSimplexGrid1([2,3,4])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(grid_3d)))
   
   V,CV = larSimplexGrid1([1,1,1])
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS((V,CV))))
   SK2 = (V,larSimplexFacets(CV))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK2)))
   SK1 = (V,larSimplexFacets(SK2[1]))
   VIEW(EXPLODE(1.5,1.5,1.5)(MKPOLS(SK1)))
   
