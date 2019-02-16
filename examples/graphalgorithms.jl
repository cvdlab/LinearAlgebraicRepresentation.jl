using LinearAlgebraicRepresentation
using Plasm, DataStructures
using PyCall
Lar = LinearAlgebraicRepresentation
p = PyCall.pyimport("pyplasm")

# Example 1

EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],
[6,7],[4,7],[2,7],[5,8],[1,8],[3,8]]

VV = Lar.verts2verts(EV)
V = hcat([[2.,1],[2,2],[1,2],[0,3],[0,0],[3,0],[3,3],[1,1]]...)

spanningtree, fronds = Lar.depth_first_search(EV)

hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])

# Example 2

spanningtree, fronds = Lar.depth_first_search(EV,5)

hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])


"""
	edge_biconnect(EV::Lar.Cells)
	
Compute maximal 2-edge-connected components of an undirected graph.

# Example

```julia
V = hcat([[1.,2],[1,3],[0,3],[1,4],[2,3],[2,2],[1,1],[1,0],[0,1]]...)
EV = [[3,4],[2,4],[2,3],[2,5],[1,2],[5,6],[1,6],[1,7],[7,9],[8,9],[7,8]]
using Plasm
Plasm.view( Plasm.numbering(1.3)((V,[[[k] for k=1:size(V,2)], EV])) )

julia> edge_biconnect(EV::Lar.Cells)
4-element Array{Array{Array{Int64,1},1},1}:
 [[4, 2], [3, 4], [2, 3]]        
 [[6, 1], [5, 6], [2, 5], [1, 2]]
 [[9, 7], [8, 9], [7, 8]]        
 [[1, 7]]                        
```
"""
function edge_biconnect(EV::Lar.Cells)
    VV = Lar.verts2verts(EV)
	number = zeros(Int, length(VV))
	lowpt = zeros(Int, length(VV))
	components = Array{Array{Int,1},1}[]

	function biconnect(v, u)
		number[v] = i; i +=  1
		lowpt[v] = number[v]
		for w in VV[v] 
			if number[w] == 0  # w is not yet numbered 
				push!(stack, [v, w]) # add [v, w] to edges stack
				
				biconnect(w, v)
				
				lowpt[v] = min(lowpt[v], lowpt[w])
				if lowpt[w] ≥ number[v] 
					# start new biconnected component
					component = Array{Int,1}[]
					while number[stack[end][1]] ≥ number[w]
						u1,u2 = pop!(stack) 
						push!(component, [u1,u2])
					end
					v,w = pop!(stack)
					push!(component, [v,w])
					push!(components, component)
				end
			else 
				if number[w] < number[v] && w ≠ u
					push!(stack, [v,w]) # add (v, w) to edge stack;
					lowpt[v] = min(lowpt[v], number[w]);
				end
			end
		end
	end 
	
	i = 1
	stack = Array{Int64,1}[]
	for w=1:length(VV) 
		if number[w] == 0 # w is not yet numbered  
			biconnect(w, 0)
		end
	end
	return components
end
