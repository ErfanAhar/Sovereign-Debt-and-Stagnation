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

function fortran_initial_guess(model::Model; debt_weight::Float64 = 0.10, c_floor::Float64 = 1e-8)
    Nb, Ng, Ne = model.Nb, model.Ng, model.Ne
    y = reshape(model.g, Ng, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, Ne))
    imf_net = (1 - model.R_l) .* model.l_policy

    vnd = Array{Float64}(undef, Nb, Ng, Ne)
    for gi in 1:Ng, ei in 1:Ne
        c = max.(y[gi, ei] .+ imf_net[gi, ei] .- debt_weight .* model.b, c_floor)
        vnd[:, gi, ei] .= u(c, model.gamma) ./ (1 - model.beta)
    end

    vd = zeros(Float64, Nb, Ng, Ne)
    d = falses(Nb, Ng, Ne)
    e = trues(Nb, Ng, Ne)
    Q = ones(Float64, Nb, Ng, Ne)
    X = ones(Float64, Nb, Ng, Ne)
    R = ones(Float64, Nb, Ng)
    n = zeros(Float64, Nb, Ng)
    pdefault = zeros(Float64, Nb, Ng)
    b_policy_idx = repeat(reshape(collect(1:Nb), Nb, 1, 1), 1, Ng, Ne)

    return Solution(
        vnd,
        vd,
        d,
        e,
        Q,
        X,
        R,
        n,
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
