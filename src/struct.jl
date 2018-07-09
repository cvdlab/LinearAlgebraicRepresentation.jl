


"""
	approxVal(PRECISION)(value)
	
Transform the float `value` to get a `PRECISION` number of significant digits.
"""
function approxVal(PRECISION)
	PRECISION = 7
    function approxVal0(value)
        out = round(value*(10^(PRECISION)))/10^(PRECISION)
        if out == -0.0
            out = 0.0
        end
        return out 
    end
    return approxVal0
end



"""
	t(args...)

"""
function t(args...)
	d = length(args)
	mat = eye(d+1)
	for k in range(1,d)
        	mat[k,d+1]=args[k]
	end
	return mat
end



"""
	s(args...)

"""
function s(args...)
	d = length(args)
	mat = eye(d+1)
	for k in range(1,d)
		mat[k,k]=args[k]
	end
	return mat
end




"""
	r(args...)

"""

function r(args...)
    args = collect(args)
    n = length(args)
    if n == 1 # rotation in 2D
        angle = args[1]; COS = cos(angle); SIN = sin(angle)
        mat = eye(3)
        mat[1,1] = COS;    mat[1,2] = -SIN;
        mat[2,1] = SIN;    mat[2,2] = COS;
    end

     if n == 3 # rotation in 3D
        mat = eye(4)
        angle = norm(args); axis = normalize(args)
        COS = cos(angle); SIN= sin(angle)
        if axis[2]==axis[3]==0.0    # rotation about x
            mat[2,2] = COS;    mat[2,3] = -SIN;
            mat[3,2] = SIN;    mat[3,3] = COS;

        elseif axis[1]==axis[3]==0.0   # rotation about y
            mat[1,1] = COS;    mat[1,3] = SIN;
            mat[3,1] = -SIN;    mat[3,3] = COS;
        elseif axis[1]==axis[2]==0.0    # rotation about z
            mat[1,1] = SIN;    mat[1,2] = -SIN;
            mat[2,1] = COS;    mat[2,2] = COS;
        
        else
	    I = eye(3); u = axis
	    Ux=[0 -u[3] u[2] ; u[3] 0 -u[1] ;  -u[2] u[1] 1]
	    UU =[u[1]*u[1]    u[1]*u[2]   u[1]*u[3];
             u[2]*u[1]    u[2]*u[2]   u[2]*u[3];
             u[3]*u[1]    u[3]*u[2]   u[3]*u[3]]
	    mat[1:3,1:3]=COS*I+SIN*Ux+(1.0-COS)*UU
		end
	end

return mat

end



"""
	removeDups(CW)

"""
function removeDups(CW)
	CW = collect(Set(CW))
	CWs = collect(map(sort,CW))
	return CWs
end 



"""
	Struct

"""

type Struct
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
	struct2lar(structure)

"""
function struct2lar(structure)
	listOfModels = evalStruct(structure)
	vertDict= Dict()
	index,defaultValue,W,CW,FW = 0,0,Array{Float64,1}[],Array{Int64,1}[],Array{Int64,1}[]
	
	for model in listOfModels
		if  length(model)==2
			V,FV = model
		elseif length(model)==3
			V,FV,EV = model
		end
		
		for (k,incell) in enumerate(FV)
			outcell=[]
			for v in incell
				key = map(approxVal(7), V[:,v])
				if get(vertDict,key,defaultValue)==defaultValue
					index += 1
                   	vertDict[key]=index
					append!(outcell,index)
					append!(W,[key])                   
				else
					append!(outcell,vertDict[key])
				end
			end
			append!(CW,[outcell])
		end
		
		if length(model)==3
		
			for (k,incell) in enumerate(EV)
				outcell=[]
				for v in incell
					key = map(approxVal(7), V[:,v])
					if get(vertDict,key,defaultValue)==defaultValue
						index += 1
						vertDict[key]=index
						append!(outcell,[index])
						append!(W,[key])                   
					else
						append!(outcell,vertDict[key])
					end
				end
				append!(FW,[outcell])
			end
			
		end
	end
	
	if length(listOfModels[end])==2
		CW = removeDups(CW)
		return hcat(W...),CW
	end
	
	if length(listOfModels[end])==3
		FW = removeDups(FW)
		CW = removeDups(CW)
		return hcat(W...),CW,FW
	end
end





"""
	embedTraversal(cloned,obj,n,suffix)

"""

function embedTraversal(cloned,obj,n,suffix)

	for i in range(1,len(obj))
		if isa(obj.body[i],Matrix)
			mat = obj.body[i]
			d,d = size(mat)
			newMat = eye(d+n*1)
			for h in range(1,d-1)
				for k in range(1,d-1)
					newMat[h,k]=mat[h,k]
				end
				newMat[h,d-1+n*1]=mat[h,d-1]
			end
			append!(cloned.body,newMat)
		elseif (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==3 
			V,FV,EV = obj.body[i]
			dimadd = fill([0.0],n)
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,[(V,FV,EV)])
		elseif  (isa(obj.body[i],Tuple) ||isa(obj.body[i],Array))&& length(obj.body[i])==2 
			V,EV = deepcopy(obj.body[i])
			dimadd = fill([0.0],n)
			for k in dimadd
				for v in V
					append!(v,k)
				end
			end
			append!(cloned.body,[(V,EV)])
		elseif isa(obj.body[i],Struct)
			newObj = Struct()
			newObj.box = hcat((obj.body[i].box,[fill([0],n),fill([0],n)]))
			newObj.category = obj.body[i].category
			append!(cloned.body,embedTraversal(newObj,obj.body[i],n,suffix))
		end
	end
	return cloned
end



"""
	embedStruct(n)

"""

function embedStruct(n)
	function embedStruct0(self,suffix="New")
		if n==0
			return self, length(self.box[1])
		end
		cloned = Struct()
		cloned.box = hcat((self.box,[fill([0],n),fill([0],n)]))	
		cloned.name = string(object_id(cloned))
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
		return Array[theMin,theMax]

	elseif (isa(model,Tuple) ||isa(model,Array))&& (length(model)==2 || length(model)==3)
		V = model[1]
		theMin = minimum(V, 2)
		theMax = maximum(V, 2)
	end

	return [theMin,theMax]
end

 


"""
	apply(affineMatrix)(larmodel)

"""
function apply(affineMatrix)
	function apply0(larmodel)
		data = collect(larmodel)
		V = data[1]
		
		m,n = size(V)
		W = [V; fill(1.0, (1,n))]
		V = (affineMatrix * W)[1:m,1:n]

		data[1] = V	
		larmodel = Tuple(data)
	end 
	return apply0
end

#function apply(affineMatrix)
#	function apply0(larmodel::Struct)
#		larmodel = Struct([affineMatrix,larmodel])
#	end 
#	return apply0
#end




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

function traversal(CTM,stack,obj,scene=[])
	for i in 1:len(obj)
		if isa(obj.body[i],Matrix)
			CTM = CTM*obj.body[i]
		elseif (isa(obj.body[i],Tuple) || isa(obj.body[i],Array)) && (length(obj.body[i])==2 || length(obj.body[i])==3)
			l = apply(CTM)(obj.body[i])
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
function evalStruct(self)
	dim = checkStruct(self.body)
   	CTM, stack = eye(dim+1), []
   	scene = traversal(CTM, stack, self, []) 
return scene
end




