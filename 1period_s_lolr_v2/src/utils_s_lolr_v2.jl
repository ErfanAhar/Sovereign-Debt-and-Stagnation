using QuantEcon
using LinearAlgebra

function init_model(m::Model = Model())
    sigma_u = 1.0
    mc = QuantEcon.tauchen(m.Ne, 0.0, sigma_u, 0.0, m.n_std)
    eps = collect(mc.state_values)
    P_eps = mc.p
    pi_eps = vec((P_eps^10_000)[1, :])
    pi_eps ./= sum(pi_eps)

    pi_g = fill(1.0 / length(m.g), length(m.g))
    for _ in 1:10_000
        pi_g = vec(pi_g' * m.P_g)
    end
    g_bar = dot(pi_g, m.g)

    lbar_dgs = m.Δ_dgs .* g_bar
    l_max = max(m.l_max, maximum(lbar_dgs))
    l = collect(range(m.l_min, l_max, length = m.Nl))

    return Model(
        Nb = m.Nb,
        Nl = m.Nl,
        Ng = length(m.g),
        Ns = size(m.P_s, 1),
        Ne = m.Ne,
        b_min = m.b_min,
        b_max = m.b_max,
        l_min = m.l_min,
        l_max = l_max,
        b = collect(range(m.b_min, m.b_max, length = m.Nb)),
        l = l,
        g = copy(m.g),
        eps = eps,
        P_g = copy(m.P_g),
        P_s = copy(m.P_s),
        pi_eps = pi_eps,
        K = kron(m.P_s, m.P_g),
        beta = m.beta,
        gamma = m.gamma,
        sigma_eps = m.sigma_eps,
        n_std = m.n_std,
        rstar = m.rstar,
        theta = m.theta,
        kappa = m.kappa,
        phi_g = copy(m.phi_g),
        nbar_g = copy(m.nbar_g),
        R_l = m.R_l,
        Δ_dgs = copy(m.Δ_dgs),
        lbar_dgs = lbar_dgs,
        pub = m.pub,
        max_iter = m.max_iter,
        max_iter_vd = m.max_iter_vd,
        max_iter_x = m.max_iter_x,
        tol = m.tol,
        tol_vd = m.tol_vd,
        tol_x = m.tol_x,
    )
end

@inline function u(c::Real, gamma::Real)
    if c <= 0
        return -Inf
    elseif gamma == 1.0
        return log(c)
    else
        return c^(1 - gamma) / (1 - gamma)
    end
end

function nearest_index(grid::AbstractVector, values::AbstractArray)
    idx = similar(values, Int)
    n = length(grid)
    for i in eachindex(values)
        @inbounds j = searchsortedfirst(grid, values[i])
        if j <= 1
            @inbounds idx[i] = 1
        elseif j > n
            @inbounds idx[i] = n
        else
            @inbounds idx[i] = abs(grid[j] - values[i]) < abs(grid[j - 1] - values[i]) ? j : j - 1
        end
    end
    return idx
end

function expected_next!(out::AbstractArray{Float64,3}, v::AbstractArray{<:Real,4},
    K::AbstractMatrix, pi_eps::AbstractVector, tmp_vec::Vector{Float64})
    Nb, Ng, Ns, Ne = size(v)
    @assert size(out) == (Nb, Ng, Ns)
    @assert length(tmp_vec) == Nb * Ng * Ns

    v_mat = reshape(v, Nb * Ng * Ns, Ne)
    mul!(tmp_vec, v_mat, pi_eps)
    v_eps = reshape(tmp_vec, Nb, Ng * Ns)
    out_mat = reshape(out, Nb, Ng * Ns)
    mul!(out_mat, v_eps, K')
    return out
end

function expected_next!(out::AbstractArray{Float64,4}, v::AbstractArray{<:Real,5},
    K::AbstractMatrix, pi_eps::AbstractVector, tmp_vec::Vector{Float64})
    Nb, Nl, Ng, Ns, Ne = size(v)
    @assert size(out) == (Nb, Nl, Ng, Ns)
    @assert length(tmp_vec) == Nb * Nl * Ng * Ns

    v_mat = reshape(v, Nb * Nl * Ng * Ns, Ne)
    mul!(tmp_vec, v_mat, pi_eps)
    v_eps = reshape(tmp_vec, Nb * Nl, Ng * Ns)
    out_mat = reshape(out, Nb * Nl, Ng * Ns)
    mul!(out_mat, v_eps, K')
    return out
end

function max_abs_diff_safe(new::AbstractArray{<:Real}, old::AbstractArray{<:Real})
    @assert size(new) == size(old)
    err = 0.0
    for i in eachindex(new, old)
        @inbounds a = new[i]
        @inbounds b = old[i]
        if isfinite(a) && isfinite(b)
            err = max(err, abs(a - b))
        elseif a == b
            continue
        else
            err = Inf
        end
    end
    return err
end

function damp_update_safe!(dest::AbstractArray{Float64}, old::AbstractArray{Float64}, new::AbstractArray{Float64}, damp::Float64)
    @assert size(dest) == size(old) == size(new)
    for i in eachindex(dest, old, new)
        @inbounds a = old[i]
        @inbounds b = new[i]
        if isfinite(a) && isfinite(b)
            @inbounds dest[i] = (1.0 - damp) * a + damp * b
        else
            @inbounds dest[i] = b
        end
    end
    return dest
end

function initial_solution(model::Model)
    Nb, Nl, Ng, Ns, Ne = model.Nb, model.Nl, model.Ng, model.Ns, model.Ne
    y = reshape(model.g, 1, 1, Ng, 1, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, 1, 1, 1, Ne))
    b = reshape(model.b, Nb, 1, 1, 1, 1)
    l = reshape(model.l, 1, Nl, 1, 1, 1)
    c0 = repeat(max.(y .- 0.10 .* b .- l, 1e-8), 1, 1, 1, Ns, 1)
    vnd = u.(c0, model.gamma) ./ (1 - model.beta)
    vd = vnd .- 1.0
    d = falses(Nb, Nl, Ng, Ns, Ne)
    e = trues(Nb, Nl, Ng, Ns, Ne)
    Q = ones(Nb, Nl, Ng, Ns, Ne)
    X = ones(Nb, Nl, Ng, Ns, Ne)
    R = ones(Nb, Nl, Ng, Ns)
    n = zeros(Nb, Nl, Ng, Ns)
    pdefault = zeros(Nb, Nl, Ng, Ns)
    schedule_mask = trues(Nb, Nl, Ng, Ns)
    b_policy_idx = repeat(reshape(collect(1:Nb), Nb, 1, 1, 1, 1), 1, Nl, Ng, Ns, Ne)
    l_policy_idx = repeat(reshape(collect(1:Nl), 1, Nl, 1, 1, 1), Nb, 1, Ng, Ns, Ne)
    l_policy_idx_d = copy(l_policy_idx)
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
        schedule_mask,
        b_policy_idx,
        l_policy_idx,
        l_policy_idx_d,
        Float64[],
        Float64[],
        Float64[],
    )
end
