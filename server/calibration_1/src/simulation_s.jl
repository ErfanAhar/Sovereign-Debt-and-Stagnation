using Random
using Statistics

@inline function _draw_index(cdf::AbstractVector, u::Real, n::Int)
    idx = searchsortedfirst(cdf, u)
    return idx > n ? n : idx
end

function simulate(model::Model, sol::Solution; T::Int = 200_000, b0_idx::Int = 1,
    g0_idx::Int = 1, s0_idx::Int = Int(cld(model.Ns, 2)),
    eps0_idx::Int = Int(cld(model.Ne, 2)),
    seed::Union{Nothing,Int} = nothing)

    rng = seed === nothing ? Random.default_rng() : MersenneTwister(seed)

    b = model.b
    g = model.g
    eps = model.eps

    bprime = b ./ reshape(g, 1, model.Ng)
    bprime_idx = nearest_index(b, bprime)
    kbprime_idx = nearest_index(b, model.kappa .* bprime)

    P_g_cdf = cumsum(model.P_g, dims = 2)
    P_s_cdf = cumsum(model.P_s, dims = 2)
    pi_eps_cdf = cumsum(model.pi_eps)

    @assert 1 <= b0_idx <= model.Nb
    @assert 1 <= g0_idx <= model.Ng
    @assert 1 <= s0_idx <= model.Ns
    @assert 1 <= eps0_idx <= model.Ne

    b_idx = b0_idx
    g_idx = g0_idx
    s_idx = s0_idx
    eps_idx = eps0_idx

    in_default = false

    b_path = Vector{Float64}(undef, T)
    y_path = Vector{Float64}(undef, T)
    n_path = Vector{Float64}(undef, T)
    R_path = Vector{Float64}(undef, T)
    default_path = Vector{Bool}(undef, T)
    g_idx_path = Vector{Int}(undef, T)
    s_idx_path = Vector{Int}(undef, T)
    eps_idx_path = Vector{Int}(undef, T)

    for t in 1:T
        g_val = g[g_idx]
        eps_val = eps[eps_idx]
        y_base = g_val * exp(model.sigma_eps * eps_val)

        if !in_default && sol.d[b_idx, g_idx, s_idx, eps_idx]
            in_default = true
        end

        default_path[t] = in_default
        b_path[t] = b[b_idx]
        g_idx_path[t] = g_idx
        s_idx_path[t] = s_idx
        eps_idx_path[t] = eps_idx

        if in_default
            y_path[t] = model.phi_g[g_idx] * y_base
            n_path[t] = 0.0
            R_path[t] = NaN
            b_idx_next_default = bprime_idx[b_idx, g_idx]
            kb_idx_next_default = kbprime_idx[b_idx, g_idx]
        else
            y_path[t] = y_base
            b_idx_next = sol.b_policy_idx[b_idx, g_idx, s_idx, eps_idx]
            n_path[t] = sol.n[b_idx_next, g_idx, s_idx]
            R_path[t] = sol.R[b_idx_next, g_idx, s_idx]
        end

        g_idx_next = _draw_index(view(P_g_cdf, g_idx, :), rand(rng), model.Ng)
        s_idx_next = _draw_index(view(P_s_cdf, s_idx, :), rand(rng), model.Ns)
        eps_idx_next = _draw_index(pi_eps_cdf, rand(rng), model.Ne)

        if !default_path[t]
            b_idx = b_idx_next
            in_default = false
        else
            if rand(rng) < model.theta
                in_default = true
                b_idx = b_idx_next_default
            else
                v_repay = sol.vnd[kb_idx_next_default, g_idx_next, s_idx_next, eps_idx_next]
                v_def = sol.vd[b_idx_next_default, g_idx_next, s_idx_next, eps_idx_next]
                if v_repay >= v_def
                    in_default = false
                    b_idx = kb_idx_next_default
                else
                    in_default = true
                    b_idx = b_idx_next_default
                end
            end
        end

        g_idx = g_idx_next
        s_idx = s_idx_next
        eps_idx = eps_idx_next
    end

    return (b = b_path, y = y_path, n = n_path, R = R_path, default = default_path,
        g_idx = g_idx_path, s_idx = s_idx_path, eps_idx = eps_idx_path)
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

function summarize_simulation(sim, model::Model; sol::Solution, burnin::Int = 10_000, nbins::Int = 40)
    T = length(sim.b) # total simulated periods
    @assert 0 <= burnin < T # burn-in must be a valid prefix length
    idx = (burnin + 1):T # kept sample after burn-in

    default_sample = sim.default[idx] # default indicator in estimation window
    b_sample = sim.b[idx] # current debt state path
    y_sample = sim.y[idx] # output path
    n_sample = sim.n[idx] # issuance path
    R_sample = sim.R[idx] # gross borrowing rate path
    g_idx_sample = sim.g_idx[idx] # growth-state index path
    s_idx_sample = sim.s_idx[idx] # sunspot-state index path
    eps_idx_sample = sim.eps_idx[idx] # shock-state index path

    prev_default = burnin == 0 ? false : sim.default[burnin] # default status right before kept sample
    starts = _episode_starts(default_sample; prev_default = prev_default) # flags for starts of default episodes
    episodes = sum(starts) # number of default episodes
    default_rate = episodes / length(default_sample) # episodes per period in kept sample

    mean_b_to_gdp = mean(b_sample ./ y_sample) # unconditional avg(b/y)
    mean_n_to_gdp = mean(n_sample ./ y_sample) # unconditional avg(n/y)

    nondefault = .!default_sample # non-default observations
    spread = R_sample .- (1 + model.rstar) # spread over risk-free gross rate
    spread_mask = nondefault .& .!isnan.(spread) .& .!isinf.(spread) # valid spread observations
    spreads = spread[spread_mask] # filtered spread series
    mean_spread = isempty(spreads) ? NaN : mean(spreads) # avg spread conditional on valid non-default

    b_idx_sample = nearest_index(model.b, b_sample) # map simulated b to grid indices
    bprime_idx = similar(b_idx_sample) # container for chosen next-period debt indices
    @inbounds for t in eachindex(b_idx_sample)
        bprime_idx[t] = sol.b_policy_idx[b_idx_sample[t], g_idx_sample[t], s_idx_sample[t], eps_idx_sample[t]] # optimal b' index at each simulated state
    end # fill optimal next-period debt choice

    f_sample = model.g[g_idx_sample] .* model.b[bprime_idx] # face value implied by chosen b'
    qb_sample = f_sample ./ R_sample # market value qb = f / R

    imf_net = similar(y_sample) # net IMF transfer term (1 - R_l) * l(g,eps)
    @inbounds for t in eachindex(imf_net)
        imf_net[t] = (1 - model.R_l) * model.l_policy[g_idx_sample[t], eps_idx_sample[t]] # state-contingent net IMF resources
    end # fill IMF net resources path
    tb_sample = b_sample .- n_sample .- imf_net # trade balance from flow budget identity

    low_state = g_idx_sample .== argmin(model.g) # low-growth state indicator
    high_state = g_idx_sample .== argmax(model.g) # high-growth state indicator

    default_rate_low = _safe_rate(sum(starts .& low_state), sum(low_state)) # low-growth default-episode rate
    default_rate_high = _safe_rate(sum(starts .& high_state), sum(high_state)) # high-growth default-episode rate
    all_mask = nondefault # mask for all non-default observations
    low_mask = nondefault .& low_state # mask for low-growth non-default observations
    high_mask = nondefault .& high_state # mask for high-growth non-default observations
    qb_to_y = qb_sample ./ y_sample # qb/y series
    f_to_y = f_sample ./ y_sample # f/y series
    n_to_y = n_sample ./ y_sample # n/y series
    b_to_y = b_sample ./ y_sample # b/y series
    tb_to_y = tb_sample ./ y_sample # tb/y series
    function _moment_block(mask, this_default_rate) # helper used for different conditioning sets of moments (all/low/high)
        spread_valid = mask .& .!isnan.(spread) .& .!isinf.(spread) # valid spread mask for current block
        return ( # grouped moments for one conditioning set
            avg_spread = _safe_mean(spread, spread_valid), # avg spread in current block
            avg_qb_to_y = _safe_mean(qb_to_y, mask), # avg(qb/y) in current block
            avg_f_to_y = _safe_mean(f_to_y, mask), # avg(f/y) in current block
            avg_n_to_y = _safe_mean(n_to_y, mask), # avg(n/y) in current block
            avg_b_to_y = _safe_mean(b_to_y, mask), # avg(b/y) in current block
            avg_tb_to_y = _safe_mean(tb_to_y, mask), # avg(tb/y) in current block
            default_rate = this_default_rate, # default-episode rate in current block
        ) # return moment tuple for current block
    end # _moment_block
    all_moments = _moment_block(all_mask, default_rate) # moments in all non-default observations
    low_moments = _moment_block(low_mask, default_rate_low) # moments in low-growth non-default observations
    high_moments = _moment_block(high_mask, default_rate_high) # moments in high-growth non-default observations

    edges = collect(range(model.b[1], model.b[end], length = nbins + 1)) # histogram bin edges for b
    counts = zeros(Int, nbins) # histogram counts for b
    for x in b_sample
        bin = searchsortedlast(edges, x) # locate bin index for current b
        if bin == length(edges)
            bin -= 1 # place right-endpoint observations in last valid bin
        end # adjust endpoint case
        bin = clamp(bin, 1, nbins) # guard against boundary issues
        counts[bin] += 1 # increment bin count
    end # finished histogram count loop
    total_counts = sum(counts) # total histogram mass
    density = total_counts > 0 ? counts / total_counts : zeros(Float64, nbins) # normalized histogram

    share_b_min = mean(b_sample .== model.b[1]) # fraction at lower debt bound
    share_b_max = mean(b_sample .== model.b[end]) # fraction at upper debt bound

    return (default_rate = default_rate, # overall default-episode rate
        mean_b_to_gdp = mean_b_to_gdp, # unconditional avg(b/y)
        mean_n_to_gdp = mean_n_to_gdp, # unconditional avg(n/y)
        mean_credit_spread = mean_spread, # avg spread in valid non-default sample
        credit_spread_obs = sum(spread_mask), # number of valid spread observations
        hist_edges = edges, # histogram edges
        hist_counts = counts, # histogram counts
        hist_density = density, # histogram density
        share_b_min = share_b_min, # mass at lower debt boundary
        share_b_max = share_b_max, # mass at upper debt boundary
        sample_size = length(default_sample), # post-burn sample size
        all = all_moments, # all-state non-default moments
        low = low_moments, # low-growth non-default moments
        high = high_moments) # high-growth non-default moments
end # summarize_simulation
