using SparseArrays

function tgw(model, dim)

    function select_sigma()
        for i = 1 : lo_num  if visited[i] == 1  return i  end  end
        for i = 1 : lo_num  if visited[i] == 0  return i  end  end
        return -1
    end

    function get_closer_petal(τ, stem, sign)
        if ord_angles[τ] == undef  ord_angles[τ] = CAGD.eval_ord_angle(model, dim - 1, τ)  end
        stem_idx = findfirst(isequal(stem), ord_angles[τ])
        return ord_angles[τ][(stem_idx + sign - 1, length(ord_angles[τ])) + 1]
    end

    # Prepearing tgw model
    tgw_model = CAGD.Model(model.G)
    for d = 1 : dim - 1  CAGD.addModelCells!(tgw_model, d, model.T[d])  end
    lo_cls = model.T[dim - 1]
    lo_num, llo_num = size(lo_cls)

    # Preprocessing angles between d-1 cells w.r.t. d-2 cells
    # ord_angles has for every d-2 cell the circular ordering between d-1 cells
    ord_angles = Array{Array{Int64, 1}, 1}(undef, llo_num)

    # Either 0 (not visited), 1 (once visited) or 2 (twice visited)
    visited = zeros(Int8, lo_num)
    # If visited[σ] == 1 then other_sign[σ] is the sign for the second visit
    other_sign = zeros(Int8, lo_num)

    # Components accumulator
    components = [];
    comp_idx = 0;

    # Until all dim-1 cells have been visited twice, continue
    while (σ = select_sigma()) > 0
        # Build new 3-cell
        c = SparseArrays.sparsevec(zeros(Int8, 1, lo_num))
        if visited[σ] == 0
            # Set sign in representation
            c[σ] = 1
            other_sign[σ] = -1
            # Set start to a new component
            components.push!([])
            comp_idx += 1
        else
            c[σ] = other_sign[σ]
        end
        components[comp_idx].push!(σ)
        visited[σ] += 1

        # Build other faces via corolla method
        # boundary = d-2 exposed cells of currently analysed d cell
        while (boundary = c * model.T[dim - 1]).nzind != []
            # corolla are the next d-1 cells to be added
            corolla = zeros(Int8, lo_num)
            for τ in boundary
                # petals are incident d-1 cells over τ d-2 cell
                petals = model.T[dim - 1][:, τ]
                # stem is the cell generating the cell generating the petals
                #TODO it should be one only ?!?
                stem = intersect(c.nzind, petals.nzind)[1]
                petal = CAGD.get_closer_petal(τ, stem, sign(boundary[τ]))          # -sign(...)
                # the new petal and the stem must travel τ in opposit direction
                if model.T[dim - 1][stem, τ] == - model.t[dim - 1][petal, τ]
                    petal_sign = c[stem]
                elseif model.T[dim - 1][stem, τ] == model.t[dim - 1][petal, τ]
                    petal_sign = -c[stem]
                else
                    throw(ArgumentError("Something Weird Happened"))
                end

                if corolla[petal] == 0 && c[petal] == 0
                    # if petal is new, it is simply added
                    corolla[petal] = petal_sign
                elseif corolla[petal] ==
                    (corolla[petal] == petal_sign) ||
                        throw(ArgumentError("Incoherency detected"))
                end
            end

        end


     end
end
