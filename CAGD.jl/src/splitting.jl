Lar = LinearAlgebraicRepresentation

using IntervalTrees
using SparseArrays
using NearestNeighbors
using DataStructures


#-------------------------------------------------------------------------------
#   PAIRWISE DECOMPOSITION
#-------------------------------------------------------------------------------

"""
	facesplitting(model::CAGD.Model)::CAGD.Model

Performs the pairwise decomposition of all faces with the possibly intersecting.
"""

function facesplitting(model; atol = 10e-8)
    sp_idx = CAGD.spaceIndex(model, 2)
    de_models = [
        CAGD.face_decomposition(model, face_idx, sp_idx[face_idx], atol = atol)
        for face_idx = 1 : size(model, 2, 1)
    ]
    return CAGD.uniteMultipleModels(de_models)
end

function face_decomposition(
		model::CAGD.Model,
		face_idx::Int,
		sp_idx_face::Array{Int,1};
		atol = 1e-8
	)::CAGD.Model

	function face_intersection(													# Inserire il modello in depth
			V::Lar.Points,
			EV::Lar.ChainOp,
			face::Lar.Cell
		)::CAGD.Model

		# Get indices of `face` vertices in traversal order
	    Vidx = Lar.buildFV(EV, face)
	    retV = Lar.Points(undef, 2, 0)

	    visited_verts = []
		for i in 1:length(Vidx)
			# vector to vertex
			o = V[:, Vidx[i]]
			# next vertex (circular ordering)
			j = i < length(Vidx) ? i+1 : 1
			# vector from i-th to j-th
	        d = V[:, Vidx[j]] - o

			err = atol
			
			# if the point is coplanar, it should be added
			if (-err < d[3] < err)
				#retV = [retV V[1:2, Vidx[i]]]                                  Da aggiungere solo se la faccia è tutta parallela
				#retV = [retV V[1:2, Vidx[j]]]
			else

				# check whether the z-ratio between o and d is ∈ [-1, 0]
				# then there is an intersextion point `p = o + α * d`
	            alpha = -o[3] / d[3]

	            if -err <= alpha <= 1+err
	                p = o + alpha*d

	                if -err < alpha < err || 1-err < alpha < 1+err
	                    if !(Lar.vin(p, visited_verts))
	                        push!(visited_verts, p)
	                        retV = [retV reshape(p, 3, 1)[1:2, 1]]
	                    end
	                else
	                    retV = [retV reshape(p, 3, 1)[1:2, 1]]
	                end
	            end
	        end

	    end

	    vnum = size(retV, 2)

	    if vnum == 1
	        vnum = 0
	        retV = Lar.Points(undef, 2, 0)
	    end
	    enum = (÷)(vnum, 2)
	    retEV = spzeros(Int8, enum, vnum)

	    for i in 1:enum
	        retEV[i, 2*i-1:2*i] = [-1, 1]
	    end

	    retModel = CAGD.Model(retV)
		CAGD.addModelCells!(retModel, 1, retEV)
		return retModel
	end

	#== THE FUNCTION ==#

	Gface, Gfaceidx = CAGD.getModelCellGeometry(model, 2, face_idx, true)
	M = CAGD.build_projection_matrix(Gface)
	PG = (M * [model.G; ones(1, size(model, 0, 2))])[1:3, :]

	# Generate face Model (in 2D)
	Pmodel = CAGD.Model(PG[1:2, Gfaceidx])
	CAGD.addModelCells!(Pmodel, 1, model.T[1][CAGD.getModelLoCell(model, 2, face_idx), Gfaceidx])

	# Each other face adds new pieces
	for f in sp_idx_face
		# println("Intersecting $face_idx with $f")
		f_FE = model.T[2][f, :]
		f_EV = model.T[1][CAGD.getModelLoCell(model, 2, f), :]
		_, vs, _ = findnz(f_EV);
		if *((abs.(PG[3, vs]) .< atol)...)
			newVidx = sort(unique(vs));
			oldVidx = SparseArrays.spzeros(Int, size(PG, 2), 1)
			for i = 1 : length(newVidx)  oldVidx[newVidx[i]] = i  end
			newModel = CAGD.Model(PG[1:2, newVidx])
			I, J, Val = findnz(f_EV)
			f_EV = SparseArrays.sparse(I, map(x -> oldVidx[x], J), Val)#, max(I...), f_EV.n)
			CAGD.addModelCells!(newModel, 1, f_EV)
			#newModel = CAGD.getSubModel(model, 2, f)
		else
			newModel = face_intersection(PG, model.T[1], f_FE)
		end
		CAGD.uniteModels!(Pmodel, newModel)
	end

	Pmodel = CAGD.mergemodel(Pmodel, signed_merge = false)

	model = CAGD.planar_arrangement(Pmodel, sparsevec(ones(Int8, length(Gfaceidx))))
	@assert !isnothing(model) "UNEXPECTED ERROR: a face should be mapped to itself"

	n = size(model, 0, 2)
	retModel = CAGD.Model((inv(M)*[model.G; zeros(1, n); ones(1, n)])[1:3, :])
	CAGD.addModelCells!(retModel, 1, model.T[1])
	CAGD.addModelCells!(retModel, 2, model.T[2])

	return retModel
end
