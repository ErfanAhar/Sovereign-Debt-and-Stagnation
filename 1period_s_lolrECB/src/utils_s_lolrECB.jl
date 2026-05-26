using LinearAlgebra
using QuantEcon

# Build the debt and shock grids and combine growth and sunspot transitions.
# This performs Step 0 of the solution algorithm in the note.
function init_model(m::Model = Model())
    mc = QuantEcon.tauchen(m.Ne, 0.0, 1.0, 0.0, m.n_std)
    eps = collect(mc.state_values)
    pi_eps = vec((mc.p^10_000)[1, :])
    pi_eps ./= sum(pi_eps)

    return Model(
        Nb = m.Nb,
        Ng = length(m.g),
        Ns = size(m.P_s, 1),
        Ne = m.Ne,
        b_min = m.b_min,
        b_max = m.b_max,
        b = collect(range(m.b_min, m.b_max, length = m.Nb)),
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
        ecb_spread = m.ecb_spread,
        intervention_prob = m.intervention_prob,
        pub = m.pub,
        max_iter = m.max_iter,
        max_iter_vd = m.max_iter_vd,
        max_iter_x = m.max_iter_x,
        tol = m.tol,
        tol_vd = m.tol_vd,
        tol_x = m.tol_x,
    )
end

# CRRA flow utility for one consumption value; infeasible consumption gets -Inf.
@inline function u(c::Real, gamma::Real)
    if c <= 0
        return -Inf
    elseif gamma == 1.0
        return log(c)
    end
    return c^(1 - gamma) / (1 - gamma)
end

# Apply CRRA utility elementwise to candidate consumption values.
function u(c::AbstractArray, gamma::Real)
    out = similar(c, Float64)
    for i in eachindex(c)
        @inbounds out[i] = u(c[i], gamma)
    end
    return out
end

# Map off-grid debt choices to their closest debt-grid positions.
function nearest_index(grid::AbstractVector, values::AbstractArray)
    idx = similar(values, Int)
    for i in eachindex(values)
        value = values[i]
        j = searchsortedfirst(grid, value)
        if j <= 1
            idx[i] = 1
        elseif j > length(grid)
            idx[i] = length(grid)
        else
            idx[i] = abs(grid[j] - value) < abs(grid[j - 1] - value) ? j : j - 1
        end
    end
    return idx
end

# Integrate a value or payoff array over next-period growth, sunspot, and shock.
function expected_next(v::AbstractArray{<:Real,4}, model::Model)
    Nb, Ng, Ns, Ne = size(v)
    @assert (Ng, Ns, Ne) == (model.Ng, model.Ns, model.Ne)
    v_eps = reshape(v, Nb * Ng * Ns, Ne) * model.pi_eps
    out = reshape(v_eps, Nb, Ng * Ns) * model.K'
    return reshape(out, Nb, Ng, Ns)
end

# Measure the largest finite fixed-point change while preserving infeasible values.
function max_abs_diff_safe(new::AbstractArray{<:Real}, old::AbstractArray{<:Real})
    @assert size(new) == size(old)
    err = 0.0
    for i in eachindex(new, old)
        a = new[i]
        b = old[i]
        if isfinite(a) && isfinite(b)
            err = max(err, abs(a - b))
        elseif a != b
            return Inf
        end
    end
    return err
end

# Damp the repayment-value update to stabilize the outer fixed point.
function damp_update_safe!(old::AbstractArray{Float64}, new::AbstractArray{Float64}, damp::Float64)
    for i in eachindex(old, new)
        a = old[i]
        b = new[i]
        old[i] = isfinite(a) && isfinite(b) ? (1 - damp) * a + damp * b : b
    end
    return old
end

# Return the ECB purchase price P(d, g, s) in Section 2.1.
# Purchases occur only during repayment, low growth, and a bad sunspot.
@inline function ecb_purchase_price(model::Model, gi::Int, si::Int)
    active_state = gi == argmin(model.g) && si == 1
    return active_state ? 1 / (1 + model.rstar + model.ecb_spread) : 0.0
end

# Return the exogenous probability that a lender sells a bond to the ECB.
@inline function ecb_sale_probability(model::Model, gi::Int, si::Int)
    active_state = gi == argmin(model.g) && si == 1
    return active_state ? model.intervention_prob : 0.0
end

# Construct initial values, prices, schedules, and policy indexes.
# This supplies the Step 1 guess for repayment value vnd.
function initial_solution(model::Model)
    Nb, Ng, Ns, Ne = model.Nb, model.Ng, model.Ns, model.Ne
    y = reshape(model.g, 1, Ng, 1, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, 1, 1, Ne))
    b = reshape(model.b, Nb, 1, 1, 1)
    c0 = repeat(max.(y .- 0.10 .* b, 1e-8), 1, 1, Ns, 1)
    vnd = u(c0, model.gamma) ./ (1 - model.beta)
    vd = vnd .- 1.0
    d = falses(Nb, Ng, Ns, Ne)
    e = trues(Nb, Ng, Ns, Ne)
    Q = ones(Nb, Ng, Ns, Ne)
    X = ones(Nb, Ng, Ns, Ne)
    H = ones(Nb, Ng, Ns) ./ (1 + model.rstar)
    L = copy(H)
    R = fill(1 + model.rstar, Nb, Ng, Ns)
    n = zeros(Nb, Ng, Ns)
    pdefault = zeros(Nb, Ng, Ns)
    schedule_mask = trues(Nb, Ng, Ns)
    ecb_used = falses(Nb, Ng, Ns)
    b_policy_idx = repeat(reshape(collect(1:Nb), Nb, 1, 1, 1), 1, Ng, Ns, Ne)
    return Solution(vnd, vd, d, e, Q, X, H, L, R, n, pdefault, schedule_mask,
        ecb_used, b_policy_idx, Float64[], Float64[], Float64[])
end
