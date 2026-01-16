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

# Expectation over (g', s', eps') with iid eps and Markov sunspot.
function expected_next_iid(v::AbstractArray{<:Real,4}, P_g::AbstractMatrix, P_s::AbstractMatrix, pi_eps::AbstractVector)
    Nb, Ng, Ns, Ne = size(v)
    @assert size(P_g, 1) == Ng && size(P_g, 2) == Ng
    @assert size(P_s, 1) == Ns && size(P_s, 2) == Ns
    @assert length(pi_eps) == Ne

    v_mat = reshape(v, Nb * Ng * Ns, Ne)
    v_eps = v_mat * pi_eps
    v_eps = reshape(v_eps, Nb, Ng, Ns)

    vE = similar(v_eps, Float64)
    Ps_t = P_s'
    for b in 1:Nb
        vE[b, :, :] = P_g * v_eps[b, :, :] * Ps_t
    end
    return vE
end
