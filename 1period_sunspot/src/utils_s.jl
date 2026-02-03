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
function expected_next_iid(v::AbstractArray{<:Real,4}, K::AbstractMatrix, pi_eps::AbstractVector)
    Nb, Ng, Ns, Ne = size(v)
    @assert size(K, 1) == Ng * Ns && size(K, 2) == Ng * Ns
    @assert length(pi_eps) == Ne

    v_mat = reshape(v, Nb * Ng * Ns, Ne)
    v_eps = v_mat * pi_eps
    v_eps = reshape(v_eps, Nb, Ng, Ns)

    v_eps_mat = reshape(v_eps, Nb, Ng * Ns)
    vE_mat = v_eps_mat * K'
    return reshape(vE_mat, Nb, Ng, Ns)
end
