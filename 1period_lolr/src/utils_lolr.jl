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

function fortran_initial_guess(model::Model; debt_weight::Float64 = 0.10, c_floor::Float64 = 1e-8)
    Nb, Nl, Ng, Ne = model.Nb, model.Nl, model.Ng, model.Ne
    y = reshape(model.g, Ng, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, Ne))
    l0_idx = findmin(abs.(model.l))[2]

    vnd = Array{Float64}(undef, Nb, Nl, Ng, Ne)
    for gi in 1:Ng, ei in 1:Ne, li in 1:Nl
        c = max.(y[gi, ei] .- debt_weight .* model.b .- model.l[li], c_floor)
        vnd[:, li, gi, ei] .= u.(c, model.gamma) ./ (1 - model.beta)
    end

    vd = zeros(Float64, Nb, Nl, Ng, Ne)
    d = falses(Nb, Nl, Ng, Ne)
    e = trues(Nb, Nl, Ng, Ne)
    Q = ones(Float64, Nb, Nl, Ng, Ne)
    X = ones(Float64, Nb, Nl, Ng, Ne)
    R = ones(Float64, Nb, Nl, Ng)
    n = zeros(Float64, Nb, Nl, Ng)
    pdefault = zeros(Float64, Nb, Nl, Ng)
    b_policy_idx = repeat(reshape(collect(1:Nb), Nb, 1, 1, 1), 1, Nl, Ng, Ne)
    l_policy_idx = fill(l0_idx, Nb, Nl, Ng, Ne)
    l_policy_idx_d = fill(l0_idx, Nb, Nl, Ng, Ne)

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
        l_policy_idx,
        l_policy_idx_d,
        0,
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
    )
end
