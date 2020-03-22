Lar = LinearAlgebraicRepresentation

using IntervalTrees
using SparseArrays
using NearestNeighbors
using DataStructures


#-------------------------------------------------------------------------------
#   PAIRWISE DECOMPOSITION
#-------------------------------------------------------------------------------

"""
	pairwise_decomposition(model::CAGD.Model)::CAGD.Model

Performs the pairwise decomposition of all faces with the possibly intersecting.
"""

function pairwise_decomposition(model)
    sp_idx = CAGD.spaceIndex(model, 2)
    de_models = [
        CAGD.face_decomposition(model, face_idx, sp_idx[face_idx])
        for face_idx = 1 : size(model, 2, 1)
    ]
    return CAGD.uniteMultipleModels(de_models)
end

function face_decomposition(
		model::CAGD.Model,
		face_idx::Int,
		sp_idx_face::Array{Int,1}
	)::CAGD.Model

	function build_projection_matrix(vs::Lar.Points)::Matrix{Float64}
    	u1 = vs[:, 2] - vs[:, 1]
    	u2 = vs[:, 3] - vs[:, 1]
    	u3 = Lar.cross(u1, u2)
    	T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    	T[1:3, 4] = - vs[:,1]
    	M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    	M[1:3, 1:3] = [u1 u2 u3]
    	return M' * T
	end

	function face_intersection(													# Inserire il modello in depth
			V::Lar.Points,
			EV::Lar.ChainOp,
			face::Lar.Cell
		)::CAGD.Model

	    vs = Lar.buildFV(EV, face)
	    retV = Lar.Points(undef, 2, 0)

	    visited_verts = []
	    for i in 1:length(vs)
	        o = V[:, vs[i]]
	        j = i < length(vs) ? i+1 : 1
	        d = V[:, vs[j]] - o

	        err = 10e-8
	        if !(-err < d[3] < err)

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
	M = build_projection_matrix(Gface)
	PG = (M * [model.G; ones(1, size(model, 0, 2))])[1:3, :]

	# Generate face Model (in 2D)
	Pmodel = CAGD.Model(PG[1:2, Gfaceidx])
	CAGD.addModelCells!(Pmodel, 1, model.T[1][CAGD.getModelLoCell(model, 2, face_idx), Gfaceidx])

	# Each other face adds new pieces
	for f in sp_idx_face
		# println("Intersecting $face_idx with $f")
		newModel = face_intersection(PG, model.T[1], model.T[2][f, :])
		CAGD.uniteModels!(Pmodel, newModel)
	end

	Pmodel = CAGD.mergeModelVertices(Pmodel, signed = false)

	model = CAGD.planar_arrangement(Pmodel, sparsevec(ones(Int8, length(Gfaceidx))))
	@assert !isnothing(model) "UNEXPECTED ERROR: a face should be mapped to itself"

	n = size(model, 0, 2)
	retModel = CAGD.Model((inv(M)*[model.G; zeros(1, n); ones(1, n)])[1:3, :])
	CAGD.addModelCells!(retModel, 1, model.T[1])
	CAGD.addModelCells!(retModel, 2, model.T[2])

	return retModel
end
