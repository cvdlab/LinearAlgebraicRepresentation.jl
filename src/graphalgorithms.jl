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





