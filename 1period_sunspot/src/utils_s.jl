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

function fortran_initial_guess(model::Model; debt_weight::Float64 = 0.10, c_floor::Float64 = 1e-8)
    Nb, Ng, Ns, Ne = model.Nb, model.Ng, model.Ns, model.Ne
    y = reshape(model.g, Ng, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, Ne))
    imf_net = (1 - model.R_l) .* model.l_policy

    vnd = Array{Float64}(undef, Nb, Ng, Ns, Ne)
    for gi in 1:Ng, ei in 1:Ne
        c = max.(y[gi, ei] .+ imf_net[gi, ei] .- debt_weight .* model.b, c_floor)
        vnd_ge = u(c, model.gamma) ./ (1 - model.beta)
        for si in 1:Ns
            vnd[:, gi, si, ei] .= vnd_ge
        end
    end

    vd = zeros(Float64, Nb, Ng, Ns, Ne)
    d = falses(Nb, Ng, Ns, Ne)
    e = trues(Nb, Ng, Ns, Ne)
    Q = ones(Float64, Nb, Ng, Ns, Ne)
    X = ones(Float64, Nb, Ng, Ns, Ne)
    R = ones(Float64, Nb, Ng, Ns)
    n = zeros(Float64, Nb, Ng, Ns)
    schedule_mask = trues(Nb, Ng, Ns)
    pdefault = zeros(Float64, Nb, Ng, Ns)
    b_policy_idx = repeat(reshape(collect(1:Nb), Nb, 1, 1, 1), 1, Ng, Ns, Ne)

    return Solution(
        vnd,
        vd,
        d,
        e,
        Q,
        X,
        R,
        n,
        schedule_mask,
        pdefault,
        b_policy_idx,
        0,
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
    )
end
