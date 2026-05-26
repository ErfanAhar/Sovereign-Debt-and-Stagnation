using Random
using Statistics

# Draw a discrete state index from a cumulative probability vector.
@inline function _draw_index(cdf::AbstractVector, draw::Real)
    min(searchsortedfirst(cdf, draw), length(cdf))
end

# Simulate debt, output, default, borrowing rates, and ECB support.
function simulate(model::Model, sol::Solution; T::Int = 200_000, b0_idx::Int = 1,
    g0_idx::Int = 1, s0_idx::Int = 2, eps0_idx::Int = Int(cld(model.Ne, 2)),
    seed::Union{Nothing,Int} = nothing)
    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)
    bprime_idx = nearest_index(model.b, model.b ./ reshape(model.g, 1, model.Ng))
    kbprime_idx = nearest_index(model.b, model.kappa .* model.b ./ reshape(model.g, 1, model.Ng))
    P_g_cdf = cumsum(model.P_g, dims = 2)
    P_s_cdf = cumsum(model.P_s, dims = 2)
    pi_eps_cdf = cumsum(model.pi_eps)

    b_idx, g_idx, s_idx, eps_idx = b0_idx, g0_idx, s0_idx, eps0_idx
    in_default = false
    b_path = Vector{Float64}(undef, T)
    y_path = Vector{Float64}(undef, T)
    n_path = zeros(T)
    R_path = fill(NaN, T)
    default_path = falses(T)
    ecb_used_path = falses(T)
    ecb_sale_prob_path = zeros(T)
    g_idx_path = Vector{Int}(undef, T)
    s_idx_path = Vector{Int}(undef, T)
    eps_idx_path = Vector{Int}(undef, T)

    for t in 1:T
        y_base = model.g[g_idx] * exp(model.sigma_eps * model.eps[eps_idx])
        if !in_default && sol.d[b_idx, g_idx, s_idx, eps_idx]
            in_default = true
        end

        # Record current states and outcomes before choosing next-period debt.
        b_path[t] = model.b[b_idx]
        y_path[t] = in_default ? model.phi_g[g_idx] * y_base : y_base
        default_path[t] = in_default
        g_idx_path[t], s_idx_path[t], eps_idx_path[t] = g_idx, s_idx, eps_idx

        if in_default
            # Default prevents new private issuance; inherited debt evolves mechanically.
            b_next_default = bprime_idx[b_idx, g_idx]
            kb_next_default = kbprime_idx[b_idx, g_idx]
        else
            # Repayment selects bprime and may induce ECB secondary-market support.
            b_next = sol.b_policy_idx[b_idx, g_idx, s_idx, eps_idx]
            n_path[t] = sol.n[b_next, g_idx, s_idx]
            R_path[t] = sol.R[b_next, g_idx, s_idx]
            ecb_used_path[t] = sol.ecb_used[b_next, g_idx, s_idx]
            if ecb_used_path[t]
                ecb_sale_prob_path[t] = ecb_sale_probability(model, g_idx, s_idx)
            end
        end

        # Draw next-period exogenous states.
        gp = _draw_index(view(P_g_cdf, g_idx, :), rand(rng))
        sp = _draw_index(view(P_s_cdf, s_idx, :), rand(rng))
        ep = _draw_index(pi_eps_cdf, rand(rng))

        # Update debt and default status, allowing re-entry with a haircut.
        if !in_default
            b_idx = b_next
        elseif rand(rng) < model.theta
            b_idx = b_next_default
        elseif sol.vnd[kb_next_default, gp, sp, ep] >= sol.vd[b_next_default, gp, sp, ep]
            in_default = false
            b_idx = kb_next_default
        else
            b_idx = b_next_default
        end
        g_idx, s_idx, eps_idx = gp, sp, ep
    end

    return (b = b_path, y = y_path, n = n_path, R = R_path, default = default_path,
        ecb_used = ecb_used_path, ecb_sale_prob = ecb_sale_prob_path, g_idx = g_idx_path,
        s_idx = s_idx_path, eps_idx = eps_idx_path)
end

# Compute a conditional mean when the conditioning set is nonempty.
@inline _safe_mean(x, mask) = any(mask) ? mean(x[mask]) : NaN

# Compute a conditional frequency when the conditioning set is nonempty.
@inline _safe_rate(numer, denom) = denom > 0 ? numer / denom : NaN

# Summarize simulated default, borrowing, and ECB-intervention moments after burn-in.
function summarize_simulation(sim, model::Model; burnin::Int = 10_000)
    @assert 0 <= burnin < length(sim.b)
    idx = (burnin + 1):length(sim.b)
    default = sim.default[idx]
    nondefault = .!default
    starts = copy(default)
    starts[2:end] .= default[2:end] .& .!default[1:end-1]
    if burnin > 0 && sim.default[burnin]
        starts[1] = false
    end

    spread = sim.R[idx] .- (1 + model.rstar)
    valid_spread = nondefault .& isfinite.(spread)
    crisis = (sim.g_idx[idx] .== argmin(model.g)) .& (sim.s_idx[idx] .== 1) .& nondefault
    ecb_used = sim.ecb_used[idx] .& nondefault

    return (
        default_rate = sum(starts) / length(idx),
        mean_b_to_gdp = mean(sim.b[idx] ./ sim.y[idx]),
        mean_n_to_gdp = mean(sim.n[idx] ./ sim.y[idx]),
        mean_credit_spread = _safe_mean(spread, valid_spread),
        crisis_share = mean(crisis),
        ecb_support_share = mean(ecb_used),
        ecb_support_share_in_crisis = _safe_rate(sum(ecb_used .& crisis), sum(crisis)),
        mean_ecb_sale_probability = _safe_mean(sim.ecb_sale_prob[idx], ecb_used),
        sample_size = length(idx),
    )
end
