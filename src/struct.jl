using LinearAlgebra

"""
	t(args::Array{Number,1}...)::Matrix

Return an *affine transformation Matrix* in homogeneous coordinates. Such `translation` Matrix has ``d+1`` rows and ``d+1`` columns, where ``d`` is the number of translation parameters in the `args` array.

# Examples

```julia
julia> t(1,2)			# 2D translation
# return
3×3 Array{Float64,2}:
 1.0  0.0  1.0
 0.0  1.0  2.0
 0.0  0.0  1.0

julia> t(1.,2,3)		# 3D translation
# return
4×4 Array{Float64,2}:
 1.0  0.0  0.0  1.0
 0.0  1.0  0.0  2.0
 0.0  0.0  1.0  3.0
 0.0  0.0  0.0  1.0
```
"""
function t(args...)
	d = length(args)
	mat = Matrix{Float64}(I, d+1, d+1)
	for k in range(1, length=d)
        	mat[k,d+1]=args[k]
	end
	return mat
end

"""
	s(args::Array{Number,1}...)::Matrix


Return an *affine transformation Matrix* in homogeneous coordinates. Such `scaling` Matrix has ``d+1`` rows and ``d+1`` columns, where ``d`` is the number of scaling parameters in the `args` array.

# Examples

```julia
julia> s(2,3)			# 2D scaling
# return
3×3 Array{Float64,2}:
 2.0  0.0  0.0
 0.0  3.0  0.0
 0.0  0.0  1.0

julia> s(2.,3.,4.)		# 3D scaling
# return
4×4 Array{Float64,2}:
 2.0  0.0  0.0  0.0
 0.0  3.0  0.0  0.0
 0.0  0.0  4.0  0.0
 0.0  0.0  0.0  1.0
```
"""
function s(args...)
	d = length(args)
	mat = Matrix{Float64}(I, d+1, d+1)
	for k in range(1, length=d)
		mat[k,k]=args[k]
	end
	return mat
end

"""
	r(args...)

Return an *affine transformation Matrix* in homogeneous coordinates. Such `Rotation` Matrix has *dimension* either equal to 3 or to 4, for 2D and 3D rotation, respectively.
The `{Number,1}` of `args` either contain a single `angle` parameter in *radiants*, or a vector with three elements, whose `norm` is the *rotation angle* in 3D and whose `normalized value` gives the direction of the *rotation axis* in 3D.

# Examples

```julia
julia> r(pi/6)				# 2D rotation of ``π/6`` angle
# return
3×3 Array{Float64,2}:
 0.866025  -0.5       0.0
 0.5        0.866025  0.0
 0.0        0.0       1.0

julia> r(0,0,pi/4)
# return
4×4 Array{Float64,2}:		# 3D rotation about the ``z`` axis, with ``π/6`` angle
 0.707107  -0.707107  0.0  0.0
 0.707107   0.707107  0.0  0.0
 0.0        0.0       1.0  0.0
 0.0        0.0       0.0  1.0
 
julia> r(1,1,1)		# 3D rotation about the ``x=y=z`` axis, with angle ``1.7320508`` angle
# return
4×4 Array{Float64,2}:
  0.226296  -0.183008   0.956712  0.0
  0.956712   0.226296  -0.183008  0.0
 -0.183008   0.956712   1.21332   0.0
  0.0        0.0        0.0       1.0

```
"""
function r(args...)
    args = collect(args)
    n = length(args)
    if n == 1 # rotation in 2D
        angle = args[1]; COS = cos(angle); SIN = sin(angle)
        mat = Matrix{Float64}(I, 3, 3)
        mat[1,1] = COS;    mat[1,2] = -SIN;
        mat[2,1] = SIN;    mat[2,2] = COS;
    end

     if n == 3 # rotation in 3D
        mat = Matrix{Float64}(I, 4, 4)
        angle = norm(args); 
        if norm(args) != 0.0
			axis = normalize(args)
			COS = cos(angle); SIN= sin(angle)
			if axis[2]==axis[3]==0.0    # rotation about x
				mat[2,2] = COS;    mat[2,3] = -SIN;
				mat[3,2] = SIN;    mat[3,3] = COS;
			elseif axis[1]==axis[3]==0.0   # rotation about y
				mat[1,1] = COS;    mat[1,3] = SIN;
				mat[3,1] = -SIN;    mat[3,3] = COS;
			elseif axis[1]==axis[2]==0.0    # rotation about z
				mat[1,1] = COS;    mat[1,2] = -SIN;
				mat[2,1] = SIN;    mat[2,2] = COS;
			else
				I = Matrix{Float64}(I, 3, 3); u = axis
				Ux=[0 -u[3] u[2] ; u[3] 0 -u[1] ;  -u[2] u[1] 1]
				UU =[u[1]*u[1]    u[1]*u[2]   u[1]*u[3];
					 u[2]*u[1]    u[2]*u[2]   u[2]*u[3];
					 u[3]*u[1]    u[3]*u[2]   u[3]*u[3]]
				mat[1:3,1:3]=COS*I+SIN*Ux+(1.0-COS)*UU
			end
		end
	end
	return mat
end

"""
	removeDups(CW::Cells)::Cells

Remove dublicate `cells` from `Cells` object. Then put `Cells` in *canonical form*, i.e. with *sorted indices* of vertices in each (unique) `Cells` Array element.
"""
function removeDups(CW::Cells)::Cells
	CW = collect(Set(CW))
	CWs = collect(map(sort,CW))
	return CWs
end

"""
	Struct
	
A *container* of geometrical objects is defined by applying the function `Struct` to
the array of contained objects. Each value is defined in local coordinates and may be transformed by affine transformation tensors.

The value returned from the application of `Struct` to an `Array{Any, 1}` of `LAR` values, `matrices`, and `Struct` values is a value of 
`Struct type`.  The coordinate system of this value is the one associated with the first object of the `Struct` arguments.  Also,
the resulting geometrical value is often associated with a variable name.

The generation of containers may continue hierarchically by suitably applying `Struct`. Notice that each LAR object in a `Struct` container is transformed by each matrix before it *within the container*, going from right to left. The action of a transformation (tensor) extends to each object on its right within its own container. Whereas,  the action of a tensor does not extend outside its container, according to the semantics of *PHIGS* structures.

# Example

```julia
julia> assembly = L.Struct([L.sphere()(), L.t(3,0,-1), L.cylinder()()])
# return
LARLIB.Struct(Any[([0.0 -0.173648 … -0.336824 -0.17101; 0.0 0.0 … 0.0593912 0.0301537;
-1.0 -0.984808 … 0.939693 0.984808], Array{Int64,1}[[2, 3, 1], [4, 2, 3], [4, 3, 5], [4,
5, 6], [7, 5, 6], [7, 8, 6], [7, 9, 8], … , [1.0 0.0 0.0 3.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0
-1.0; 0.0 0.0 0.0 1.0], ([0.5 0.5 … 0.492404 0.492404; 0.0 0.0 … -0.0868241 -0.0868241;
0.0 2.0 … 0.0 2.0], Array{Int64,1}[[4, 2, 3, 1], [4, 3, 5, 6], [7, 5, 8, 6], [7, 9, 10,
8], [9, 10, 11, 12], [13, 14, 11, 12], … , [68, 66, 67, 65], [68, 69, 67, 70], [71, 69,
72, 70], [71, 2, 72, 1]])], Array{Float64,2}[[-1.0; -1.0; -1.0], [3.5; 1.0; 1.0]],
"14417445522513259426", 3, "feature")

julia> assembly.name = "simple example"
# return
"simple example"

julia> assembly
# return
LARLIB.Struct(Any[([0.0 -0.173648 … -0.336824 -0.17101; 0.0 0.0 … 0.0593912 0.0301537;
-1.0 -0.984808 … 0.939693 0.984808], Array{Int64,1}[[2, 3, 1], [4, 2, 3], [4, 3, 5], [4,
5, 6], [7, 5, 6], [7, 8, 6], … , [71, 2, 72, 1]])], Array{Float64,2}[[-1.0; -1.0; -1.0],
[3.5; 1.0; 1.0]], "simple example", 3, "feature")

julia> using LARVIEW

julia> LARVIEW.view(assembly)
```
"""
mutable struct Struct
	body::Array
	box
	name::AbstractString
	dim
	category::AbstractString
	
	function Struct()
		self = new([],Nullable{Any},"new",Nullable{Any},"feature")
		self.name = string(object_id(self))
		return self

	end

	function Struct(data::Array)
		self = Struct()
		self.body = data
		self.box = box(self)
		self.dim = length(self.box[1])
		return self
	end
	
	function Struct(data::Array,name)
		self = Struct()
		self.body=[item for item in data]
		self.box = box(self)
		self.dim = length(self.box[1])
		self.name = string(name)
		return self
	end

	function Struct(data::Array,name,category)
		self = Struct()
		self.body = [item for item in data]
		self.box = box(self)
		self.dim = length(self.box[1])
		self.name = string(name)
		self.category = string(category)
		return self
	end
	
end

	function name(self::Struct)
		return self.name
	end
	function category(self::Struct)
		return self.category
	end
	
	function len(self::Struct)
		return length(self.body)
	end
	function getitem(self::Struct,i::Int)
		return self.body[i]
	end
	function setitem(self::Struct,i,value)
		self.body[i]=value
	end
	function pprint(self::Struct)
		return "<Struct name: $(self.name)"
	end
	function set_name(self::Struct,name)
		self.name = string(name)
	end
	function clone(self::Struct,i=0)
		newObj = deepcopy(self)
		if i!=0
			newObj.name="$(self.name)_$(string(i))"
		end
		return newObj
	end
	function set_category(self::Struct,category)
		self.category = string(category)
	end

"""
	struct2lar(structure::Struct)::Union{LAR,LARmodel}

"""
function struct2lar(structure)
	listOfModels = LARLIB.evalStruct(structure)
	vertDict= Dict()
	index,defaultValue,W,FW,EW = 0,0,Array{Float64,1}[],Array{Int64,1}[],Array{Int64,1}[]
	
	for model in listOfModels
		if  length(model)==2
			V,FV = model
		elseif length(model)==3
			V,FV,EV = model
		end
		
		for incell in FV
			outcell=[]
			for v in incell
				key = map(LARLIB.approxVal(7), V[:,v])
				if get(vertDict,key,defaultValue)==defaultValue
					index += 1
                   	vertDict[key]=index
					append!(outcell,index)
					append!(W,[key])                   
				else
					append!(outcell,vertDict[key])
				end
			end
			append!(FW,[outcell])
		end
	end
	if length(listOfModels[1])==3
		for model in listOfModels
			V,FV,EV = model
			for incell in EV
				outcell=[]
				for v in incell
					key = map(LARLIB.approxVal(7), V[:,v])
					if get(vertDict,key,defaultValue)==defaultValue
						index += 1
						vertDict[key]=index
						append!(outcell,[index])
						append!(W,[key])                   
					else
						append!(outcell,vertDict[key])
					end
				end
				append!(EW,[outcell])
			end
			
		end
	end
	
	topology = listOfModels[end]
	if length(topology)==2
		#FW = removeDups(FW)
		larmodel = hcat(W...),FW
		return larmodel
	elseif length(topology)==3
		#FW = removeDups(FW)
		#EW = removeDups(EW)
		larmodel = hcat(W...),FW,EW
		return larmodel
	end
end

"""
	embedTraversal(cloned::Struct,obj::Struct,n::Int,suffix::String)

# TODO:  debug embedTraversal
"""
function embedTraversal(cloned::Struct,obj::Struct,n::Int,suffix::String)
	for i=1:length(obj.body)
		if isa(obj.body[i],Matrix)
			mat = obj.body[i]
			d,d = size(mat)
			newMat = Matrix{Float64}(I, d+n, d+n)
			for h in range(1, length=d)
				for k in range(1, length=d)
					newMat[h,k]=mat[h,k]
				end
			end
			push!(cloned.body,[newMat])
		elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==3 
			V,FV,EV = deepcopy(obj.body[i])
			dimadd = n
			ncols = size(V,2)
			nmat = zeros(dimadd,ncols)
			V = [V;zeros(dimadd,ncols)]
			push!(cloned.body,[(V,FV,EV)])
		elseif  (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==2 
			V,EV = deepcopy(obj.body[i])
			dimadd = n
			ncols = size(V,2)
			nmat = zeros(dimadd,ncols)
			V = [V;zeros(dimadd,ncols)]
			push!(cloned.body,[(V,EV)])
		elseif isa(obj.body[i],Struct)
			newObj = Struct()
			newObj.box = [ [obj.body[i].box[1];zeros(dimadd)],
							[obj.body[i].box[2];zeros(dimadd)] ]
			newObj.category = (obj.body[i]).category
			push!(cloned.body,[embedTraversal(newObj,obj.body[i],obj.dim+n,suffix)])
		end
	end
	return cloned
end

"""
	embedStruct(n::Int)(self::Struct,suffix::String="New")

# TODO:  debug embedStruct
"""
function embedStruct(n::Int)
	function embedStruct0(self::Struct,suffix::String="New")
		if n==0
			return self, length(self.box[1])
		end
		cloned = Struct()
		cloned.box = [ [self.body[i].box[1];zeros(dimadd)],
						[self.body[i].box[2];zeros(dimadd)] ]
		cloned.name = self.name*suffix
		cloned.category = self.category
		cloned.dim = self.dim+n
		cloned = embedTraversal(cloned,self,n,suffix)
		return cloned
	end
	return embedStruct0
end

"""
	box(model)

"""
function box(model)
	if isa(model,Matrix)
		return nothing
	elseif isa(model,Struct)
		listOfModels = evalStruct(model)
		#dim = checkStruct(listOfModels)
		theMin,theMax = box(listOfModels[1])
		for theModel in listOfModels[2:end]
			modelMin,modelMax= box(theModel)
			for (k,val) in enumerate(modelMin)
				if val < theMin[k]
					theMin[k]=val
				end
			end
			for (k,val) in enumerate(modelMax)
				if val > theMax[k]
					theMax[k]=val
				end
			end
		end
		return [theMin,theMax]

	elseif (isa(model,Tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
		V = model[1]
		theMin = minimum(V, 2)
		theMax = maximum(V, 2)
	end

	return [theMin,theMax]
end

"""
	apply(affineMatrix::Array{Float64,2}, larmodel::Union{LAR,LARmodel})

"""
function apply(affineMatrix, larmodel)
	data = collect(larmodel)
	V = data[1]
	
	m,n = size(V)
	W = [V; fill(1.0, (1,n))]
	V = (affineMatrix * W)[1:m,1:n]

	data[1] = V	
	larmodel = Tuple(data)
	return larmodel
end

"""
	checkStruct(lst)

"""
function checkStruct(lst)
	obj = lst[1]
	if isa(obj,Matrix)
		dim = size(obj)[1]-1
	elseif (isa(obj,Tuple) || isa(obj,Array))
		dim = length(obj[1][:,1])
	
	elseif isa(obj,Struct)
		dim = length(obj.box[1])
	end
	return dim
end

"""
	traversal(CTM,stack,obj,scene=[])

"""
function traversal(CTM::Matrix, stack, obj, scene=[])
	for i = 1:length(obj.body)
		if isa(obj.body[i],Matrix)
			CTM = CTM*obj.body[i]
		elseif (isa(obj.body[i],Tuple) || isa(obj.body[i],Array)) && 
			(length(obj.body[i])==2 || length(obj.body[i])==3)
			l = apply(CTM, obj.body[i])
			push!(scene,l)
		elseif isa(obj.body[i],Struct)
			push!(stack,CTM)	
			traversal(CTM,stack,obj.body[i],scene)
			CTM = pop!(stack)
		end
	end
	return scene
end

"""
	evalStruct(self)

"""
function evalStruct(self::Struct)
	dim = checkStruct(self.body)
   	CTM, stack = Matrix{Float64}(I, dim+1, dim+1), []
   	scene = traversal(CTM, stack, self, []) 
	return scene
end