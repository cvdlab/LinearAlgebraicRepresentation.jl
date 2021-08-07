using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
import Base.cat
using DataStructures


"""
	apply(fun::Function, params::Array)

"""
function apply(fun::Function, params::Array)
	return fun(params)
end
function apply(fun::Function)
	function apply0(params::Array)
		return fun(params)
	end
	return apply0
end

"""
	apply(affineMatrix::Matrix)(larmodel::LAR)::LAR

Apply the `affineMatrix` parameter to the vertices of `larmodel`.

# Example

```
julia> square = Lar.cuboid([1,1])
([0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0], Array{Int64,1}[[1, 2, 3, 4]])

julia> Lar.apply(Lar.t(1,2))(square)
([1.0 1.0 2.0 2.0; 2.0 3.0 2.0 3.0], Array{Int64,1}[[1, 2, 3, 4]])
```
"""
function apply(affineMatrix)
	function apply0(larmodel)
		return Lar.struct2lar(Lar.Struct([ affineMatrix,larmodel ]))
	end
	return apply0
end



"""
	comp(funs::Array)

Standard mathematical composition.

Pipe execution from right to left on application to actual parameter.
"""
function comp(funs)
    function compose(f,g)
	  return x -> f(g(x))
	end
    id = x->x
    return reduce(compose, funs; init=id)
end



"""
	cons(funs::Array)(x::Any)::Array

*Construction* functional of FL and PLaSM languages.

Provides a *vector* functional that returns the array of
applications of component functions to actual parameter.

# Example

```
julia> Lar.cons([cos,sin])(0)
2-element Array{Float64,1}:
 1.0
 0.0
```
"""
function cons(funs)
	return x -> [f(x) for f in funs]
end



"""
	k(Any)(x)

*Constant* functional of FL and PLaSM languages.

Gives a constant functional that always returns the actual parameter
when applied to another parameter.

#	Examples

```
julia> Lar.k(10)(100)
10

julia> Lar.k(sin)(cos)
sin
```
"""
function k(Any)
	x->Any
end


"""
	aa(fun::Function)(args::Array)::Array

AA applies fun to each element of the args sequence

# Example

```
julia> Lar.aa(sqrt)([1,4,9,16])
4-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
```
"""
function aa(fun)
	function aa1(args::Array)
		map(fun,args)
	end
	return aa1
end



"""
	id(x::Anytype)

Identity function.  Return the argument.

"""
id = x->x




"""
	distr(args::Union{Tuple,Array})::Array

Distribute right. The parameter `args` must contain a `list` and an element `x`.
Return the `pair` array with the elements of `args` coupled with `x`

# Example

```
julia> Lar.distr(([1,2,3],10))
3-element Array{Array{Int64,1},1}:
 [1, 10]
 [2, 10]
 [3, 10]
```
"""
function distr(args)
	list,element = args
	return [ [e,element] for e in list ]
end




"""
	distl(args::Union{Tuple,Array})::Array

Distribute right. The parameter `args` must contain an element `x` and a `list`.
Return the `pair` array with `x` coupled with the elements of `args`.

# Example

```
julia> Lar.distl((10, [1,2,3]))
3-element Array{Array{Int64,1},1}:
 [10, 1]
 [10, 2]
 [10, 3]
```
"""
function distl(args)
	element, list = args
	return [ [element, e] for e in list ]
end
