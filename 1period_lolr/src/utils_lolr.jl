using LinearAlgebra

@inline function u(c::Real, gamma::Real)
    # NOTE: This assumes gamma != 1.0 (CRRA). If gamma == 1.0, replace with log utility.
    if c <= 0
        return -Inf
    end
    return c^(1 - gamma) / (1 - gamma)
end

function nearest_index(grid::AbstractVector, values::AbstractArray)
    idx = similar(values, Int)
    n = length(grid)
    for i in eachindex(values)
        v = values[i]
        j = searchsortedfirst(grid, v)
        if j <= 1
            idx[i] = 1
        elseif j > n
            idx[i] = n
        else
            idx[i] = (abs(grid[j] - v) < abs(grid[j - 1] - v)) ? j : (j - 1)
        end
    end
    return idx
end

# Expectation over (g', eps') with iid eps (independent of current eps).
function expected_next_iid(v::AbstractArray{<:Real,3}, P_g::AbstractMatrix, pi_eps::AbstractVector)
    Nb, Ng, Ne = size(v)
    @assert size(P_g, 1) == Ng && size(P_g, 2) == Ng
    @assert length(pi_eps) == Ne
    v_mat = reshape(v, Nb * Ng, Ne)
    v_eps = v_mat * pi_eps
    v_eps = reshape(v_eps, Nb, Ng)
    return v_eps * P_g'
end

# Expectation over (g', eps') with iid eps for arrays with LOLR state.
function expected_next_iid(v::AbstractArray{<:Real,4}, P_g::AbstractMatrix, pi_eps::AbstractVector)
    Nb, Nl, Ng, Ne = size(v)
    @assert size(P_g, 1) == Ng && size(P_g, 2) == Ng
    @assert length(pi_eps) == Ne

    v_mat = reshape(v, Nb * Nl * Ng, Ne)
    v_eps = v_mat * pi_eps
    v_eps = reshape(v_eps, Nb * Nl, Ng)
    vE_mat = v_eps * P_g'
    return reshape(vE_mat, Nb, Nl, Ng)
end

# In-place expectation over (g', eps') with iid eps for 3D arrays.
function expected_next_iid!(out::AbstractMatrix, v::AbstractArray{<:Real,3},
    P_g::AbstractMatrix, pi_eps::AbstractVector, tmp_vec::AbstractVector)
    Nb, Ng, Ne = size(v)
    @assert size(out, 1) == Nb && size(out, 2) == Ng
    @assert size(P_g, 1) == Ng && size(P_g, 2) == Ng
    @assert length(pi_eps) == Ne
    @assert length(tmp_vec) == Nb * Ng

    v_mat = reshape(v, Nb * Ng, Ne)
    mul!(tmp_vec, v_mat, pi_eps)
    v_eps = reshape(tmp_vec, Nb, Ng)
    mul!(out, v_eps, P_g')
    return out
end
