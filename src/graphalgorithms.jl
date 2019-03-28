
#---------------------------------------------------------------------
#	Biconnected components
#---------------------------------------------------------------------

""" 
	verts2verts(EV::Lar.Cells)::Lar.Cells

Adjacency lists of vertices of a cellular 1-complex.

# Example

```julia 
julia> V,(VV,EV,FV) = Lar.cuboidGrid([3,3],true);

julia> verts2verts(EV::Lar.Cells)::Lar.Cells
16-element Array{Array{Int64,1},1}:
 [2, 5]         
 [1, 3, 6]      
 [2, 4, 7]      
 [3, 8]         
 [1, 6, 9]      
 [2, 5, 7, 10]  
 [3, 6, 8, 11]  
 [4, 7, 12]     
 [5, 10, 13]    
 [6, 9, 11, 14] 
 [7, 10, 12, 15]
 [8, 11, 16]    
 [9, 14]        
 [10, 13, 15]   
 [11, 14, 16]   
 [12, 15]       
```
"""
function verts2verts(EV::Lar.Cells)::Lar.Cells
    cscEV = Lar.characteristicMatrix(EV)
    cscVE = convert(Lar.ChainOp, cscEV')
    cscVV = cscVE * cscEV   
    rows,cols,data = findnz(cscVV)
    VV = [ findnz(cscVV[k,:])[1] for k=1:size(cscVV,2) ]
    return [setdiff(vs,[k]) for (k,vs) in enumerate(VV)] # remove reflexive pairs
end



""" 
	DFV_visit( VV::cells, out::Array, count::Int, visited::Array, parent::Array, d::Array, low::Array, stack::Array, u::Int )::Array

Hopcroft-Tarjan algorithm. Depth-First-Visit recursive algorithm.
"""
function DFV_visit( VV::Lar.Cells, out::Array, count::Int, visited::Array, parent::Array, d::Array, low::Array, stack::Array, u::Int )::Array
		
    visited[u] = true
    count += 1
    d[u] = count
    low[u] = d[u]
    for v in VV[u]
        if ! visited[v]
            push!(stack, [(u,v)])
            parent[v] = u
            DFV_visit( VV,out,count,visited,parent,d,low,stack, v )
            if low[v] >= d[u]
                push!(out, [outputComp(stack,u,v)])
            end
            low[u] = min( low[u], low[v] )
        else
            if ! (parent[u]==v) && (d[v] < d[u])
                push!(stack, [(u,v)])
            end
            low[u] = min( low[u], d[v] )
        end
    end
    out
end

visited=[]; count=0; out=[]; parent=[]; d=[];low=[]; stack=[]; u=1

"""
	outputComp(stack::Array, u::Int, v::Int)::Array

Internal output of biconnected components. Utility function for DFV_visit graph algorithm (Hopcrof, Tarjan (1973)).
"""
function outputComp(stack::Array, u::Int, v::Int)::Array
    out = []
    while true
        e = pop!(stack)[1]
        push!(out,e)
        if e == (u,v) 
        	break
        end
    end
    return [out] 
end


""" 
	biconnectedComponent(model)

Main procedure for computation of biconnected components of a graph 
represented by a LAR unsigned pair (Points, Cells).

# Example

```julia
julia> V =  [3.0  6.0  6.0  6.0  6.0  0.0  3.0  10.0  10.0  10.0  0.0  3.0  0.0;
 			 2.0  2.0  0.0  4.0  8.0  4.0  4.0   4.0   0.0   8.0  8.0  0.0  0.0]

julia> EV = [[1, 2], [2, 3], [2, 4], [5, 10], [3, 12], [5, 11], [6, 7], [8, 9], 
			 [6, 13], [12, 13], [1, 7], [1, 12], [4, 7], [4, 8], [4, 5], [8, 10]]
			 
julia> VV = [[k] for k=1:size(V,2)]

julia> Plasm.view( Plasm.numbering(3)((V,[VV, EV])) )
```

```julia
model = (V,EV)
V,EVs = biconnectedComponent(model)
EW = convert(Lar.Cells, cat(EVs))
Plasm.view( Plasm.numbering(3)((V,[VV,EW])) )
```

Another example:

```julia
V = [0. 1  0  0 -1;
	 0  0  1 -1  0]
EV = [[1,2],[2,3],[1,4],[4,5],[5,1],[1,3]  ,[2,4]]
model = (V,EV)
Plasm.view( Plasm.numbering(1)((V,[[[k] for k=1:size(V,2)],EV])) )
V,EVs = Lar.biconnectedComponent(model)

cscEV = Lar.characteristicMatrix(EV)
biconcomp = Lar.Arrangement.biconnected_components(cscEV)

Matrix(Lar.characteristicMatrix(EV))

V,EVs = biconnectedComponent(model)
EW = convert(Lar.Cells, cat(EVs))
Plasm.view( Plasm.numbering(.3)((V,[VV,EW])) )
```

"""
function biconnectedComponent(model)
    W,EV = model
    V = collect(1:size(W,2))
    count = 0
    stack,out = [],[]
    visited = [false for v in V]
    parent = Union{Int, Array{Any,1}}[[] for v in V]
    d = Any[0 for v in V]
    low = Any[0 for v in V]    
    VV = Lar.verts2verts(EV)
    out = Any[]
    for u in V 
        if ! visited[u] 
            out = DFV_visit( VV,out,count,visited,parent,d,low,stack, u )
        end
    end
    out = [component for component in out if length(component) >= 1]
    EVs = [[map(sort∘collect,edges) for edges in cat(comp) if length(edges)>1] 
    		for comp in out]
    EVs = cat(filter(x->!isempty(x), EVs))
    EVs = sort(EVs, lt=(x,y)->length(x)>length(y))
    return W, EVs
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



