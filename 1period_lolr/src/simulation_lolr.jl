using Random
using Statistics

@inline function _draw_index(cdf::AbstractVector, u::Real, n::Int)
    idx = searchsortedfirst(cdf, u)
    return idx > n ? n : idx
end

function simulate(model::Model, sol::Solution; T::Int = 200_000, b0_idx::Int = 1,
    l0_idx::Int = 1, g0_idx::Int = 1, eps0_idx::Int = Int(cld(model.Ne, 2)),
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
    pi_eps_cdf = cumsum(model.pi_eps)

    @assert 1 <= b0_idx <= model.Nb
    @assert 1 <= l0_idx <= model.Nl
    @assert 1 <= g0_idx <= model.Ng
    @assert 1 <= eps0_idx <= model.Ne

    b_idx = b0_idx
    l_idx = l0_idx
    g_idx = g0_idx
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
    eps_idx_path = Vector{Int}(undef, T)

    for t in 1:T
        g_val = g[g_idx]
        eps_val = eps[eps_idx]
        y_base = g_val * exp(model.sigma_eps * eps_val)

        if !in_default && sol.d[b_idx, l_idx, g_idx, eps_idx]
            in_default = true
        end

        default_path[t] = in_default
        b_path[t] = b[b_idx]
        l_path[t] = l[l_idx]
        g_idx_path[t] = g_idx
        eps_idx_path[t] = eps_idx

        if in_default
            y_path[t] = model.phi_g[g_idx] * y_base
            n_path[t] = 0.0
            R_path[t] = NaN
            l_idx_next_default = sol.l_policy_idx_d[b_idx, l_idx, g_idx, eps_idx]
            n_l_path[t] = g_val * l[l_idx_next_default] / model.R_l_d
            b_idx_next_default = bprime_idx[b_idx, g_idx]
            kb_idx_next_default = kbprime_idx[b_idx, g_idx]
        else
            y_path[t] = y_base
            b_idx_next = sol.b_policy_idx[b_idx, l_idx, g_idx, eps_idx]
            l_idx_next = sol.l_policy_idx[b_idx, l_idx, g_idx, eps_idx]
            n_path[t] = sol.n[b_idx_next, l_idx_next, g_idx]
            n_l_path[t] = g_val * l[l_idx_next] / model.R_l_nd
            R_path[t] = sol.R[b_idx_next, l_idx_next, g_idx]
        end

        g_idx_next = _draw_index(view(P_g_cdf, g_idx, :), rand(rng), model.Ng)
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
                v_repay = sol.vnd[kb_idx_next_default, l_idx_next_default, g_idx_next, eps_idx_next]
                v_def = sol.vd[b_idx_next_default, l_idx_next_default, g_idx_next, eps_idx_next]
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
        eps_idx = eps_idx_next
    end

    return (b = b_path, l = l_path, y = y_path, n = n_path, n_l = n_l_path, R = R_path,
        default = default_path, g_idx = g_idx_path, eps_idx = eps_idx_path)
end

function summarize_simulation(sim, model::Model; burnin::Int = 10_000, nbins::Int = 40)
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

    episodes = 0
    prev_def = burnin == 0 ? false : sim.default[burnin]
    for d in default_sample
        if d && !prev_def
            episodes += 1
        end
        prev_def = d
    end
    default_rate = episodes / length(default_sample)

    mean_b_to_gdp = mean(b_sample ./ y_sample)
    mean_l_to_gdp = mean(l_sample ./ y_sample)
    mean_n_to_gdp = mean(n_sample ./ y_sample)
    mean_nl_to_gdp = mean(n_l_sample ./ y_sample)

    nondefault = .!default_sample
    spreads = R_sample[nondefault] .- (1 + model.rstar)
    mean_spread = isempty(spreads) ? NaN : mean(spreads)

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
    density = counts / sum(counts)

    share_b_min = mean(b_sample .== model.b[1])
    share_b_max = mean(b_sample .== model.b[end])

    return (default_rate = default_rate,
        mean_b_to_gdp = mean_b_to_gdp,
        mean_l_to_gdp = mean_l_to_gdp,
        mean_n_to_gdp = mean_n_to_gdp,
        mean_nl_to_gdp = mean_nl_to_gdp,
        mean_credit_spread = mean_spread,
        credit_spread_obs = sum(nondefault),
        hist_edges = edges,
        hist_counts = counts,
        hist_density = density,
        share_b_min = share_b_min,
        share_b_max = share_b_max,
        sample_size = length(default_sample))
end
