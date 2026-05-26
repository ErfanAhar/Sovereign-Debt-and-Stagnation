using Random
using Statistics

@inline function _draw_index(cdf::AbstractVector, u::Real, n::Int)
    idx = searchsortedfirst(cdf, u)
    return idx > n ? n : idx
end

@inline _safe_mean(x, mask) = any(mask) ? mean(x[mask]) : NaN
@inline _safe_rate(numer, denom) = denom > 0 ? numer / denom : NaN

function _episode_starts(default_path::AbstractVector{Bool}; prev_default::Bool = false)
    starts = falses(length(default_path))
    prev = prev_default
    for t in eachindex(default_path)
        starts[t] = default_path[t] && !prev
        prev = default_path[t]
    end
    return starts
end

function simulate(model::Model, sol::Solution; T::Int = 200_000, b0_idx::Int = 1,
    l0_idx::Int = 1, g0_idx::Int = 1, s0_idx::Int = Int(cld(model.Ns, 2)),
    eps0_idx::Int = Int(cld(model.Ne, 2)),
    seed::Union{Nothing,Int} = nothing)

    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

    b = model.b
    l = model.l
    g = model.g
    eps = model.eps

    bprime = b ./ reshape(g, 1, model.Ng)
    bprime_idx = nearest_index(b, bprime)
    kbprime_idx = nearest_index(b, model.kappa .* bprime)

    P_g_cdf = cumsum(model.P_g, dims = 2)
    P_s_cdf = cumsum(model.P_s, dims = 2)
    pi_eps_cdf = cumsum(model.pi_eps)

    @assert 1 <= b0_idx <= model.Nb
    @assert 1 <= l0_idx <= model.Nl
    @assert 1 <= g0_idx <= model.Ng
    @assert 1 <= s0_idx <= model.Ns
    @assert 1 <= eps0_idx <= model.Ne

    b_idx = b0_idx
    l_idx = l0_idx
    g_idx = g0_idx
    s_idx = s0_idx
    eps_idx = eps0_idx

    in_default = false

    b_path = Vector{Float64}(undef, T)
    l_path = Vector{Float64}(undef, T)
    y_path = Vector{Float64}(undef, T)
    n_path = Vector{Float64}(undef, T)
    n_l_path = Vector{Float64}(undef, T)
    R_path = Vector{Float64}(undef, T)
    default_path = Vector{Bool}(undef, T)
    g_idx_path = Vector{Int}(undef, T)
    s_idx_path = Vector{Int}(undef, T)
    eps_idx_path = Vector{Int}(undef, T)

    for t in 1:T
        g_val = g[g_idx]
        eps_val = eps[eps_idx]
        y_base = g_val * exp(model.sigma_eps * eps_val)

        if !in_default && sol.d[b_idx, l_idx, g_idx, s_idx, eps_idx]
            in_default = true
        end

        default_path[t] = in_default
        b_path[t] = b[b_idx]
        l_path[t] = l[l_idx]
        g_idx_path[t] = g_idx
        s_idx_path[t] = s_idx
        eps_idx_path[t] = eps_idx

        if in_default
            y_path[t] = model.phi_g[g_idx] * y_base
            n_path[t] = 0.0
            R_path[t] = NaN
            l_idx_next_default = sol.l_policy_idx_d[b_idx, l_idx, g_idx, s_idx, eps_idx]
            n_l_path[t] = g_val * l[l_idx_next_default] / model.R_l
            b_idx_next_default = bprime_idx[b_idx, g_idx]
            kb_idx_next_default = kbprime_idx[b_idx, g_idx]
        else
            y_path[t] = y_base
            b_idx_next = sol.b_policy_idx[b_idx, l_idx, g_idx, s_idx, eps_idx]
            l_idx_next = sol.l_policy_idx[b_idx, l_idx, g_idx, s_idx, eps_idx]
            n_path[t] = sol.n[b_idx_next, l_idx_next, g_idx, s_idx]
            n_l_path[t] = g_val * l[l_idx_next] / model.R_l
            R_path[t] = sol.R[b_idx_next, l_idx_next, g_idx, s_idx]
        end

        g_idx_next = _draw_index(view(P_g_cdf, g_idx, :), rand(rng), model.Ng)
        s_idx_next = _draw_index(view(P_s_cdf, s_idx, :), rand(rng), model.Ns)
        eps_idx_next = _draw_index(pi_eps_cdf, rand(rng), model.Ne)

        if !default_path[t]
            b_idx = b_idx_next
            l_idx = l_idx_next
            in_default = false
        else
            if rand(rng) < model.theta
                in_default = true
                b_idx = b_idx_next_default
                l_idx = l_idx_next_default
            else
                v_repay = sol.vnd[kb_idx_next_default, l_idx_next_default, g_idx_next, s_idx_next, eps_idx_next]
                v_def = sol.vd[b_idx_next_default, l_idx_next_default, g_idx_next, s_idx_next, eps_idx_next]
                if v_repay >= v_def
                    in_default = false
                    b_idx = kb_idx_next_default
                    l_idx = l_idx_next_default
                else
                    in_default = true
                    b_idx = b_idx_next_default
                    l_idx = l_idx_next_default
                end
            end
        end

        g_idx = g_idx_next
        s_idx = s_idx_next
        eps_idx = eps_idx_next
    end

    return (b = b_path, l = l_path, y = y_path, n = n_path, n_l = n_l_path, R = R_path,
        default = default_path, g_idx = g_idx_path, s_idx = s_idx_path, eps_idx = eps_idx_path)
end

function summarize_simulation(sim, model::Model; sol::Union{Nothing, Solution} = nothing, burnin::Int = 10_000, nbins::Int = 40)
    T = length(sim.b)
    @assert 0 <= burnin < T
    idx = (burnin + 1):T

    default_sample = sim.default[idx]
    b_sample = sim.b[idx]
    l_sample = sim.l[idx]
    y_sample = sim.y[idx]
    n_sample = sim.n[idx]
    n_l_sample = sim.n_l[idx]
    R_sample = sim.R[idx]
    g_idx_sample = sim.g_idx[idx]
    s_idx_sample = sim.s_idx[idx]
    eps_idx_sample = sim.eps_idx[idx]

    prev_default = burnin == 0 ? false : sim.default[burnin]
    starts = _episode_starts(default_sample; prev_default = prev_default)
    episodes = sum(starts)
    default_rate = episodes / length(default_sample)

    mean_b_to_gdp = mean(b_sample ./ y_sample)
    mean_l_to_gdp = mean(l_sample ./ y_sample)
    mean_n_to_gdp = mean(n_sample ./ y_sample)
    mean_nl_to_gdp = mean(n_l_sample ./ y_sample)

    nondefault = .!default_sample
    spread = R_sample .- (1 + model.rstar)
    spread_mask = nondefault .& .!isnan.(spread) .& .!isinf.(spread)
    spreads = spread[spread_mask]
    mean_spread = isempty(spreads) ? NaN : mean(spreads)

    if sol === nothing
        qb_to_y = fill(NaN, length(idx))
        f_to_y = fill(NaN, length(idx))
        tb_to_y = fill(NaN, length(idx))
    else
        b_idx_sample = nearest_index(model.b, b_sample)
        l_idx_sample = nearest_index(model.l, l_sample)
        bprime_idx = similar(b_idx_sample)

        @inbounds for t in eachindex(b_idx_sample)
            bprime_idx[t] = sol.b_policy_idx[
                b_idx_sample[t],
                l_idx_sample[t],
                g_idx_sample[t],
                s_idx_sample[t],
                eps_idx_sample[t],
            ]
        end

        f_sample = model.g[g_idx_sample] .* model.b[bprime_idx]
        qb_sample = b_sample ./ R_sample
        qb_to_y = qb_sample ./ y_sample
        f_to_y = f_sample ./ y_sample
        tb_to_y = (b_sample .+ l_sample .- n_sample .- n_l_sample) ./ y_sample
    end

    low_state = g_idx_sample .== argmin(model.g)
    high_state = g_idx_sample .== argmax(model.g)

    default_rate_low = _safe_rate(sum(starts .& low_state), sum(low_state))
    default_rate_high = _safe_rate(sum(starts .& high_state), sum(high_state))

    b_to_y = b_sample ./ y_sample
    l_to_y = l_sample ./ y_sample
    n_to_y = n_sample ./ y_sample
    nl_to_y = n_l_sample ./ y_sample

    function _moment_block(mask, this_default_rate)
        spread_valid = mask .& .!isnan.(spread) .& .!isinf.(spread)
        return (
            avg_spread = _safe_mean(spread, spread_valid),
            avg_qb_to_y = _safe_mean(qb_to_y, mask),
            avg_f_to_y = _safe_mean(f_to_y, mask),
            avg_b_to_y = _safe_mean(b_to_y, mask),
            avg_l_to_y = _safe_mean(l_to_y, mask),
            avg_n_to_y = _safe_mean(n_to_y, mask),
            avg_nl_to_y = _safe_mean(nl_to_y, mask),
            avg_tb_to_y = _safe_mean(tb_to_y, mask),
            default_rate = this_default_rate,
        )
    end

    all_moments = _moment_block(nondefault, default_rate)
    low_moments = _moment_block(nondefault .& low_state, default_rate_low)
    high_moments = _moment_block(nondefault .& high_state, default_rate_high)

    edges = collect(range(model.b[1], model.b[end], length = nbins + 1))
    counts = zeros(Int, nbins)
    for x in b_sample
        bin = searchsortedlast(edges, x)
        if bin == length(edges)
            bin -= 1
        end
        bin = clamp(bin, 1, nbins)
        counts[bin] += 1
    end
    total_counts = sum(counts)
    density = total_counts > 0 ? counts / total_counts : zeros(Float64, nbins)

    share_b_min = mean(b_sample .== model.b[1])
    share_b_max = mean(b_sample .== model.b[end])

    return (default_rate = default_rate,
        mean_b_to_gdp = mean_b_to_gdp,
        mean_l_to_gdp = mean_l_to_gdp,
        mean_n_to_gdp = mean_n_to_gdp,
        mean_nl_to_gdp = mean_nl_to_gdp,
        mean_credit_spread = mean_spread,
        credit_spread_obs = sum(spread_mask),
        hist_edges = edges,
        hist_counts = counts,
        hist_density = density,
        share_b_min = share_b_min,
        share_b_max = share_b_max,
        sample_size = length(default_sample),
        all = all_moments,
        low = low_moments,
        high = high_moments)
end
