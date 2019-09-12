
"""
	cubicbezier2D(curvePts::Array{Array{Float,1},1})

Produce the two ``coordinate functions`` for a ``cubic`` *Bézier curve* in ``2D``.
The input is the array `curvePts` of 4 control points with 2 coordinates.

# Example

```
julia> curvePts = [[1232.24, 1340.84],[1259.53, 1119.86],
	[1417.91, 1015.16],[1133.47, 1010.63]]

julia> bx,by = cubicbezier2D(curvePts);

julia> (bx(0), by(0))
(1232.24, 1340.84)

julia> (bx(1), by(1))
(1133.47, 1010.63)
```
"""
function cubicbezier2D( curvePts::Array{Array{Float64,1},1} )
	b1(u) = (1 - u)^3
	b2(u) = 3*u*(1 - u)^2
	b3(u) = 3*u^2*(1 - u)
	b4(u) = u^3
	cntrlverts = hcat(curvePts...)
	x = cntrlverts[1,:]
	y = cntrlverts[2,:]
	Bx = u -> x[1]*b1(u) + x[2]*b2(u) + x[3]*b3(u) + x[4]*b4(u)
	By = u -> y[1]*b1(u) + y[2]*b2(u) + y[3]*b3(u) + y[4]*b4(u)
	return Bx,By
end



"""
	const cmdsplit

Regex for splitting the `d` element of a `<path` element into
a sequence of graphics commnds.
"""
const cmdsplit = r"\s*([mMzZlLhHvVcCsSqQtTaA])\s*"


"""
	const digitRegEx

Regex for extracting numeric params from from a sequence of graphics
commnds of a `<path` element.
"""
const digitRegEx = r"\s*(-?[0-9]*\.?\d+)\s*"

"""
	pathparse(elements)

Parse graphics commands in SVG `<path` tagged element.

# Example

```julia

```
"""
function pathparse(data, cnrtlpolygon=false)
	tokens = split(data,'\"')
	pathstring = tokens[2]

	iterator = eachmatch(cmdsplit, pathstring)
	offsets = [m.offset for m in iterator]
	commands = [m.match for m in iterator]

	substrings = [pathstring[offsets[k]+1:offsets[k+1]-1]
		for k=1:length(offsets)-1]
	push!(substrings, pathstring[offsets[end]+1:end])

	params = Array{Float64,1}[]
	for substring in substrings
		numbers = []
		for match in eachmatch(digitRegEx, substring)
			push!(numbers,String(match.match))
		end
		numbers = map(x->parse(Float64,x),numbers)
		push!(params,numbers)
	end

	lines = Array{Float64,1}[]
	global startpoint, endpoint = 0.0, 0.0
	for (command, args) in zip(commands, params)
		if command == "M"
			global startpoint = args
		elseif command == "L"
			endpoint = args
			line = vcat([startpoint,endpoint]...)
			startpoint = endpoint
			push!(lines, line)
		elseif command == "C"
			contrlpt1 = args[1:2]
			contrlpt2 = args[3:4]
			endpoint = args[5:6]
			curvePts = [startpoint,contrlpt1,contrlpt2,endpoint]
			#println(curvePts)

			Bx,By = cubicbezier2D(curvePts)
			pts = [[Bx(u),By(u)] for u=0:.1:1]
			curvelines = [ vcat([pts[k],pts[k+1]]...) for k=1:length(pts)-1 ]
			for line in curvelines
				push!(lines, line)
			end

 			if cnrtlpolygon
				push!(lines,vcat(curvePts[1:2]...),vcat(curvePts[2:3]...),
					vcat(curvePts[3:4]...))
			end
			startpoint = endpoint
		end

	end
	return lines
end



"""
	lines2lar(lines)

LAR model construction from array of float quadruples.
Each `line` in input array stands for `x1,y1,x2,y2`.
"""
function lines2lar(lines)
	vertdict = OrderedDict{Array{Float64,1}, Int64}()
	EV = Array{Int64,1}[]
	idx = 0
	for h=1:size(lines,2)
		x1,y1,x2,y2 = lines[:,h]

		if ! haskey(vertdict, [x1,y1])
			idx += 1
			vertdict[[x1,y1]] = idx
		end
		if ! haskey(vertdict, [x2,y2])
			idx += 1
			vertdict[[x2,y2]] = idx
		end
		v1,v2 = vertdict[[x1,y1]],vertdict[[x2,y2]]
		push!(EV, [v1,v2])
	end
	V = hcat(collect(keys(vertdict))...)
	return V,EV
end



"""
	normalize(V::Lar.Points; flag=true::Bool)::Lar.Points

2D normalization transformation (isomorphic by defaults) of model
vertices to normalized coordinates ``[0,1]^2``. Used with SVG importing.
"""
function normalize(V::Lar.Points; flag=true)
	m,n = size(V)
	if m > n # V by rows
		V = convert(Lar.Points, V')
	end

	xmin = minimum(V[1,:]); ymin = minimum(V[2,:]);
	xmax = maximum(V[1,:]); ymax = maximum(V[2,:]);
	box = [[xmin; ymin] [xmax; ymax]]	# containment box
	aspectratio = (xmax-xmin)/(ymax-ymin)
	if flag
		if aspectratio > 1
			umin = 0; umax = 1
			vmin = 0; vmax = 1/aspectratio; ty = vmax
		elseif aspectratio < 1
			umin = 0; umax = aspectratio
			vmin = 0; vmax = 1; ty = vmax
		end
		T = Lar.t(0,ty) * Lar.s(1,-1) * Lar.s((umax-umin), (vmax-vmin)) *
			Lar.s(1/(xmax-xmin),1/(ymax-ymin)) * Lar.t(-xmin,-ymin)
	else
		# T = Lar.t(0, ymax-ymin) * Lar.s(1,-1)
		T = Lar.s(1,-1)
	end
	dim = size(V,1)
	W = T[1:dim,:] * [V;ones(1,size(V,2))]
	#V = map( x->round(x,digits=8), W )
	V = map(Lar.approxVal(8), W)

	if m > n # V by rows
		V = convert(Lar.Points, V')
	end
	return V
end



"""
	svg2lar(filename::String; flag=true)::Lar.LAR

Parse a SVG file to a `LAR` model `(V,EV)`.
Only  `<line >` and `<rect >` and `<path >` SVG primitives are currently translated.
TODO:  interpretation of transformations.
"""
function svg2lar(filename::String; flag=true)::Lar.LAR
	outlines = Array{Float64,1}[]
	matchall(r::Regex, s::AbstractString; overlap::Bool=false) =
		collect(( m.match for m=eachmatch(r,s,overlap=overlap) ));

	for line in eachline(filename)
		parts = split(lstrip(line), ' ')
		elements = [part for part in parts if part≠""]
		tag = elements[1]
		# SVG <line > primitives
		if tag == "<line"
			r = r"(.)(.)="
			regex = r"([0-9]*?\.[0-9]*)"
			prefixes = matchall(r,line)
			values = matchall(regex,line)
			outline = map(Meta.parse, values)
			#exprs = [ eval(Meta.parse(*(pre,val))) for (pre,val) in zip(prefixes, values) ]
			push!(outlines, outline)
		# SVG <rect > primitives
		elseif tag == "<rect"
			r = r"(.)(.)="
			regex = r"([0-9]*?\.[0-9]*)"
			prefixes = matchall(r,line)
			values = matchall(regex,line)
			x, y, width, height = map(Meta.parse, values)
			# regex = r"""(<rect x=")(.+?)(" y=")(.+?)(" )(.*?)( width=")(.+?)(" height=")(.+?)("/>)"""
			# coords = collect(match( regex , line)[k] for k in (4,6,8,10))
			# x, y, width, height = [ parse(Float64, string) for string in coords ]
			line1 = [ x, y, x+width, y ]
			line2 = [ x, y, x, y+height ]
			line3 = [ x+width, y, x+width, y+height ]
			line4 = [ x, y+height, x+width, y+height ]
			push!(outlines, line1, line2, line3, line4)
		# SVG <path  > primitives
		# see https://github.com/chebfun/chebfun/issues/1617
		elseif tag == "<path"
			dataregex = r"""d=(".*?")(.*?)"""
			data = string(match( dataregex , line)[1])
			polyline = pathparse(data)
			append!(outlines, polyline)
		end
	end
	lines = hcat(outlines...)
	lines = map( x->round(x,sigdigits=8), lines )
	# LAR model construction
	V,EV = lines2lar(lines)
	# normalization
	V = normalize(V,flag=flag)
	return V,EV
end
