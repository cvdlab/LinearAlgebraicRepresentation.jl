#####################################################
### From Julia 0.6 to Julia 1.0
### --- SIMPLEXN - speedup ---
### Authors: Fabio Fatelli, Emanuele Loprevite
#####################################################

include("../src/simplexn-serial.jl")
include("../src/psimplexn.jl")

using Statistics
using Plots
Plots.scalefontsizes(0.7)
gr() # loading backend

VOID = [Int64[]],[[0]] # the empty simplicial model
nt, N = 7, 5

# Compute the execution mean time
function timing(f::Function,x,n::Int64)
	t = Array{Float64}(undef,n)
	f(x...)
	for k in 1:n
		t[k] = @elapsed f(x...)
	end
	return mean(t)
end

# Create and save the plots in the specified folder
function plotting(name::String,timeS::Array{Float64,1},timeP::Array{Float64,1},flag::Int64)
	l = max(length(timeS),length(timeP))+1
	pathName = "./"*name
	s, p, xlb, ylb, DPI = "Serial", "Parallel", "Input", "Time (seconds)", 150
	plot(timeS,xlims=(1,l),title=name,label=s,xlabel=xlb,ylabel=ylb,dpi=DPI,legend=:topleft)
	savefig(pathName*"serial")
	plot(timeP,xlims=(1,l),title=name,label=p,xlabel=xlb,ylabel=ylb,dpi=DPI,legend=:topleft)
	savefig(pathName*"parallel")
	plot([timeS,timeP],xlims=(1,l),title=name,label=[s,p],xlabel=xlb,ylabel=ylb,dpi=DPI,legend=:topleft)
	savefig(pathName*"both")
end

# Create the plots
function plotting(name::String,timeS::Array{Float64,1},timeP::Array{Float64,1})
	l = max(length(timeS),length(timeP))+1
	s, p, xlb, ylb, DPI = "Serial", "Parallel", "Input", "Time (seconds)", 150
	p1 = plot(timeS,label=s,title=name)
	p2 = plot(timeP,label=p)
	p3 = plot([timeS,timeP],label=[s,p])
	plot(p1,p2,p3,xlims=(1,l),xlabel=xlb,ylabel=ylb,dpi=1.5*DPI,legend=:topleft)
end


timeSer = [timing(larExtrude1,[VOID,repeat([1,2,3,4],k^4)],nt) for k in 1:2*N]
timePar = [timing(plarExtrude1,[VOID,repeat([1,2,3,4],k^4)],nt) for k in 1:2*N]
plotting("larExtrude1",timeSer,timePar)


simp = [collect(1:500)]
sLen = length(simp[1])
for k in 2:2*N^2
	push!(simp,simp[end].+sLen)
end
timeSer = [timing(larSimplexFacets,[simp[1:k]],nt) for k in 1:2*N^2]
timePar = [timing(plarSimplexFacets,[simp[1:k]],nt) for k in 1:2*N^2]
plotting("larSimplexFacets",timeSer,timePar)


verts, quads = [[0,0,0],[0,1,0],[1,0,0],[1,1,0],[2,2,0],[2,3,0],[3,2,0],[3,3,0]], [[0,1,2,3],[4,5,6,7]]
len = length(quads[1])
for k in 2:2*N^2
	append!(verts,map(x->x.+1,verts[end-len+1:end]))
	append!(quads,[quads[end].+len])
end
timeSer = [timing(quads2tria,[(verts[1:len*k],quads[1:k])],nt) for k in 1:2*N^2]
timePar = [timing(pquads2tria,[(verts[1:len*k],quads[1:k])],nt) for k in 1:2*N^2]
plotting("quads2tria",timeSer,timePar)
