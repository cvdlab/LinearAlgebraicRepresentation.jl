using LinearAlgebraicRepresentation
using LinearAlgebra, SparseArrays, DataStructures
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL



"""
	numbering(::Number)(::Lar.Points,::Vector{OrderedDict},::GL.COLORS,::Float)
			::Vector{ViewerGL.Viewer}
	
Cell numbering of (cuboidal) complexes 3D and 2D.

## Examples

The two examples given here use the first method of the function 
`GL.numbering()()`, given in the file `src/GLText.jl`.
An use example of the below method is provided by the
`function viewsubcomplexes`. 

```
model = Lar.cuboidGrid([3,4,2], true);
GL.VIEW(GL.numbering()(model));

model = Lar.cuboidGrid([10,10], true);
meshes = GL.numbering(1.5)(model);
GL.VIEW(meshes)
```
"""
function numbering(sizeScaling=1.)
	function numbering0(V,skelDict,color=COLORS[1],alpha=0.2)
		cells = skelDict #cells are dictionaries
		meshes = []
		if length(cells)>2
@show cells
			background = GL.GLHulls(V, cells[3], color, alpha)
		end
		if size(V,1)==2
			V = GL.embed(1)(model)[1]
		end
		ns = sizeScaling
		gcode = GL.textWithAttributes("centre", 0, 0.1ns, 0.2ns, 0.025ns)
		push!(meshes,GL.GLLines(V,collect(keys(cells[2])),color))

		colors = GL.COLORS[3], GL.COLORS[7], GL.COLORS[5], GL.COLORS[8]
		
		for (h,skel) in zip(1:length(skelDict),skelDict)
			nums = []
			#for (k,cell) in zip(1:length(skel),skel)
			for key in keys(skelDict[h])
				k,cell = (skelDict[h][key],key)
				center = sum([V[:,v] for v in cell])/length(cell)
				code = GL.embed(1)( gcode(string(k)) )
				scaling = (0.6+0.1h,0.6+0.1h,1)
				push!(nums, Lar.struct2lar( Lar.Struct([
						Lar.t(center...), Lar.s(scaling...), code ]) ))
			end
			for num in nums
				mesh = GL.GLLines(num[1],num[2],colors[h])
				push!( meshes, mesh )
			end
		end
		if length(cells)>2 push!( meshes, background ) end
		return meshes
	end
	return numbering0
end


"""
	makesubsets(Model,model)

Given a `Model` of a cellular complex, and the `model` of any *p-subcomplex of it*, 
return a  `Vector` of *p+1* dictionaries **0 ≤ p ≤ d**, having the `p`-cells as keys, 
and as values the ordinal numbers of the subcomplex cells within the original complex.

"""
function makesubsets(Model,model)
	# dictionary of global unit chains
	SkelDict = []
	for (k,skel) in enumerate(Model[2])
		push!(SkelDict, DataStructures.OrderedDict(zip(skel,1:length(skel))))
	end
	# dictionary of local unit chains
	skeldict = []
	for (k,skel) in enumerate(model[2])
		push!(skeldict, DataStructures.OrderedDict(zip(skel,1:length(skel))))
	end
	# dictionary of local unit chains with global keys
	skelDict = []
	for (k,skel) in enumerate(model[2])
		vect = []
		for key in collect(keys(skeldict[k]))
			push!(vect, (key,SkelDict[k][key]))
		end
		push!(skelDict, DataStructures.OrderedDict(vect))
	end
	return skelDict
end


"""
A visualization function for all the subcomplexes of Model generated by ff, 
in this case ff is a specialization of the FF relation, but in general may by any 
expression generating a complete and well-defined sub complex of the Model parameter.
"""
function viewsubcomplexes(Model,ff,scaling)
		V, (VV,EV,FV,FE) = Model
		ffE = [union(FE[f]...) for f in ff]
		
		function viewsubcomplex(range)
			for k in range
				fVV = [[v] for v in [union(FV[f]...) for f in ff][k]]
				fEV = [EV[e] for e in ffE[k]]
				fFV = [FV[f] for f in ff[k]]
				#model = ( ([1 0 0.3; 0 1 0.2; 0 0 1] * V)[1:2,:], Lar.Cells[fVV,fEV,fFV]);
				model = (Lar.Points(V), Lar.Cells[fVV,fEV,fFV]); 
				skelDict = makesubsets(Model,model)
				GL.VIEW(push!(numbering(scaling)(V, skelDict, GL.COLORS[1], 0.5),GL.GLFrame2));
			end
		end
	
	return viewsubcomplex
end



