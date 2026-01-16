using LinearAlgebra

function u(c::Real, gamma::Real)
    if c <= 0
        return -Inf
    end
    if gamma == 1.0
        return log(c)
    end
    return c^(1 - gamma) / (1 - gamma)
end

function u(c::AbstractArray, gamma::Real)
    out = similar(c)
    pos = c .> 0
    if gamma == 1.0
        out[pos] = log.(c[pos])
    else
        out[pos] = c[pos].^(1 - gamma) ./ (1 - gamma)
    end
    out[.!pos] .= -Inf
    return out
end

function build_transition(P_g::AbstractMatrix, P_eps::AbstractMatrix)
    return kron(P_eps, P_g)
end

function build_transition(P_g::AbstractMatrix, pi_eps::AbstractVector)
    Ne = length(pi_eps)
    P_eps = repeat(reshape(pi_eps, 1, Ne), Ne, 1)
    return kron(P_eps, P_g)
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

# Gather v[bidx, g, e] for each g,e using linear indexing (no explicit loops).
function gather_b(v::Array{Float64,3}, bidx::Array{Int,2})
    Nb, Ng, Ne = size(v)
    idx3 = reshape(bidx, Nb, Ng, 1)
    base_g = reshape((0:Ng - 1) .* Nb, 1, Ng, 1)
    base_e = reshape((0:Ne - 1) .* Nb * Ng, 1, 1, Ne)
    lin = idx3 .+ base_g .+ base_e
    return v[lin]
end

# Expectation over (g', eps') conditional on current (g, eps).
function expected_next(v::Array{Float64,3}, P_ge::Matrix{Float64})
    Nb, Ng, Ne = size(v)
    v_mat = reshape(v, Nb, Ng * Ne)
    vE_mat = v_mat * P_ge'
    return reshape(vE_mat, Nb, Ng, Ne)
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
