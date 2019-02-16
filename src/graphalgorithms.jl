"""

```julia
EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],
[6,7],[4,7],[2,7],[5,8],[1,8],[3,8]]

```

"""
function verts2verts(EV::Lar.Cells)::Lar.Cells
    cscEV = Lar.characteristicMatrix(EV)
    cscVE = convert(Lar.ChainOp, cscEV')
    cscVV = cscVE * cscEV   
    rows,cols,data = SparseArrays.findnz(cscVV)
    VV = [ findnz(cscVV[k,:])[1] for k=1:size(cscVV,2) ]
    return [setdiff(vs,[k]) for (k,vs) in enumerate(VV)]
end




"""
	depth_first_search(EV::Lar.Cells [, start::Int=1])::Tuple{Lar.Cells, Lar.Cells}

*Depth-First-Search algorithm* (Tarjan, 1972) to compute a *spanning tree* of an undirected graph. The graph is entered as an array of edges, given as an array `EV` of 1-cells. Complexity linear in the number of nodes and edges: O(V, E).

Return a partition of `EV` into two arrays `spanningtree` and `fronds` which define a ``palm tree``.  The algorithm may `start` from any vertex index.

# Example

```
julia> EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],
       [6,7],[4,7],[2,7],[5,8],[1,8],[3,8]];

julia> depth_first_search(EV,8)
(Array{Int64,1}[[8, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7]], Array{Int64,1}[[6, 1], [7, 2], [7, 4], [5, 8], [3, 8]])

julia> depth_first_search(EV)
(Array{Int64,1}[[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [5, 8]], Array{Int64,1}[[6, 1], [7, 2], [7, 4], [8, 1], [8, 3]])
```
"""
function depth_first_search(EV::Lar.Cells, start::Int=1)::
		Tuple{Lar.Cells, Lar.Cells}
	# initialize the data structures and start DFS
    VV = Lar.verts2verts(EV)
    spanningtree = Array{Int,1}[];  
    fronds = Array{Int,1}[]
	number = zeros(Int, length(VV))

	function DFS(v::Int, u::Int) 
	# vertex u is the father of
	# vertex v in the spanning tree being constructed
		number[v] = i; i += 1
		for w in VV[v] 
			if number[w] == 0 
				# w is not yet numbered  
				push!(spanningtree, [v,w])
				DFS(w, v)
			elseif (number[w] < number[v]) && (w != u)
				push!(fronds, [v,w])
			end
		end
	end
	
	# initialization of DFS algorithm
	i = 1
	DFS(start, 1)
	return spanningtree::Lar.Cells, fronds::Lar.Cells
end


"""
	edge_biconnect(EV::Lar.Cells)
	
Compute maximal 2-edge-connected components of an undirected graph.

# Examples

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

```julia
julia> V = hcat([[0,4],[4,4],[4,0],[2,0],[0,0],[2,2],[1,2],[3,2],[1,3],[3,3]]...)
julia> EV = [[1,2],[9,10],[1,5],[7,9],[8,10],[2,3],[6,7],[6,8],[4,6],[4,5],[3,4]]
julia> Plasm.view( Plasm.numbering(1.3)((V,[[[k] for k=1:size(V,2)], EV])) )

julia> EVs = edge_biconnect(EV::Lar.Cells)
3-element Array{Array{Array{Int64,1},1},1}:
 [[8, 6], [10, 8], [9, 10], [7, 9], [6, 7]]
 [[4, 6]]                                  
 [[5, 1], [4, 5], [3, 4], [2, 3], [1, 2]]  

julia> map(Lar.depth_first_search, EVs)
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



