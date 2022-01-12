#
#	This file contains some LinearAlgebraicRepresentation function that
#	 have not been updated to the last version yet. The problem is raised
#	 the `cat` acronym, here renamed `pcat`.
#	This file will be removed as soon as the appropriate changes are made.
#

function chainbasis2polygons(V,copEV,copFE)
	FE = [findnz(copFE[k,:])[1] for k=1:copFE.m]
	EV = [findnz(copEV[k,:])[1] for k=1:copEV.m]
    
	FEs = Array{Int64,1}[]
	EVs = Array{Array{Int64,1},1}[]
	FVs = Array{Int64,1}[]
	for f=1:copFE.m
		push!( FEs, collect(Set(CAGD.pcat([e for e in FE[f]]))) )
		# edges in EVs are aggregated by face, in order to answer point-classifications
		push!( EVs, [EV[e] for e in FE[f]]  )
		push!( FVs, collect(Set(CAGD.pcat([EV[e] for e in FE[f]])))  )
	end
	polygons = collect(zip(EVs,FVs,FEs)) # (EV,FV,FFE)s for polygon
	W = convert(Lar.Points,V')
	return W,polygons,FE
end


function internalpoints2d(W,copEV,copFE) # W by rows
	# transform each 2-cell in a solid (via Lar model)
	#-------------------------------------------------------------------------------
	U,pols,FE = chainbasis2polygons(W,copEV,copFE) # V by rows
	# compute, for each `pol` (3-cell) in `pols`, one `internalpoint`.
	#-------------------------------------------------------------------------------
	internalpoints = []
	for f=1:length(pols)
		(EV,FV,FE) = pols[f]
		#GL.VIEW([ GL.GLFrame, GL.GLLines(V,EV) ]);
		internalpoint = Lar.getinternalpoint2d(W,EV,FV, f, copEV,copFE)
		push!(internalpoints,internalpoint)
	end
	return internalpoints
end