function _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx, l, l0_idx)
    Nb, Nl, Ng, Ne = size(vd)
    err = Inf
    iters = 0
    l_policy_idx_d = zeros(Int, Nb, Nl, Ng, Ne)
    l_prime = l
    l_allowed = l_prime .<= model.lbar_d
    if !any(l_allowed)
        l_allowed[l0_idx] = true
    end
    w_reentry = Array{Float64}(undef, Nb, Ng, Ne)
    cont = Array{Float64}(undef, Nb, Ng, Ne)
    contE_tmp = Array{Float64}(undef, Nb, Ng)
    contE_tmp_vec = Vector{Float64}(undef, Nb * Ng)

    for _ in 1:model.max_iter_vd
        iters += 1
        vd_new = similar(vd)
        for gi in 1:Ng
            bidx = bprime_idx[:, gi]
            kbidx = kbprime_idx[:, gi]
            contE_mat = zeros(Nb, Nl)
            for lprimej in 1:Nl
                @views vnd_kb = vnd[kbidx, lprimej, :, :]
                @views vd_b = vd[bidx, lprimej, :, :]
                @. w_reentry = max(vnd_kb, vd_b)
                @. cont = model.theta * vd_b + (1 - model.theta) * w_reentry
                expected_next_iid!(contE_tmp, cont, model.P_g, model.pi_eps, contE_tmp_vec)
                @views contE_mat[:, lprimej] .= contE_tmp[:, gi]
            end

            n_l_vec = (model.g[gi] .* l_prime) ./ model.R_l_d
            g_contE = g_beta[gi] .* contE_mat
            Threads.@threads for tidx in 1:(Ne * Nl)
                ei = (tidx - 1) ÷ Nl + 1
                lj = tidx - (ei - 1) * Nl # current l index
                y_def = model.phi_g[gi] * y[gi, ei]
                util_vec = similar(n_l_vec)
                @. util_vec = u(y_def + n_l_vec - l[lj], model.gamma)
                @inbounds for bi in 1:Nb
                    best_val = -Inf
                    best_l = l0_idx
                    @inbounds for lprimej in 1:Nl
                        if l_allowed[lprimej]
                            val = g_contE[bi, lprimej] + util_vec[lprimej]
                            if val > best_val
                                best_val = val
                                best_l = lprimej
                            end
                        end
                    end
                    vd_new[bi, lj, gi, ei] = best_val
                    l_policy_idx_d[bi, lj, gi, ei] = best_l
                end
            end
        end
        err = maximum(abs.(vd_new .- vd))
        vd .= vd_new
        if err < model.tol_vd
            break
        end
    end
    return err, iters, l_policy_idx_d
end

function _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx, l_policy_idx_d)
    Nb, Nl, Ng, Ne = size(X)
    err = Inf
    iters = 0

    for _ in 1:model.max_iter_x
        iters += 1
        Q = (1 .- d) .+ d .* X
        X_new = similar(X)
        for gi in 1:Ng
            bidx = bprime_idx[:, gi]
            kbidx = kbprime_idx[:, gi]
            for ei in 1:Ne
                lidx = l_policy_idx_d[:, :, gi, ei]
                termE = zeros(Nb, Nl)
                for gpi in 1:Ng
                    Pg = model.P_g[gi, gpi]
                    for epi in 1:Ne
                        weight = Pg * model.pi_eps[epi]
                        for bi in 1:Nb
                            bpi = bidx[bi]
                            kbpi = kbidx[bi]
                            for li in 1:Nl
                                lpi = lidx[bi, li]
                                X_next = X[bpi, lpi, gpi, epi]
                                e_next = e[bpi, lpi, gpi, epi]
                                Q_kb = Q[kbpi, lpi, gpi, epi]
                                term = model.theta * X_next + (1 - model.theta) * ((1 - e_next) * X_next + e_next * (model.kappa * Q_kb))
                                termE[bi, li] += weight * term
                            end
                        end
                    end
                end
                X_new[:, :, gi, ei] .= (1 / (1 + model.rstar)) .* termE
            end
        end
        err = maximum(abs.(X_new .- X))
        X .= X_new
        if err < model.tol_x
            break
        end
    end
    Q = (1 .- d) .+ d .* X
    return Q, err, iters
end

function _compute_schedule(b, l, g, model, Q, d)
    Nb = length(b)
    Nl = length(l)
    Ng = length(g)
    QE = expected_next_iid(Q, model.P_g, model.pi_eps)
    R = (1 + model.rstar) ./ max.(QE, 1e-8)
    gb = reshape(b, Nb, 1, 1) .* reshape(g, 1, 1, Ng)
    n = gb ./ R
    pdefault = expected_next_iid(Float64.(d), model.P_g, model.pi_eps)
    return QE, R, n, pdefault
end

function _update_vnd!(vnd, vd, model, b, l, y, g_beta, n, pdefault, b0_idx, l0_idx)
    Nb, Nl, Ng, Ne = size(vnd)
    w = max.(vnd, vd)
    wE = expected_next_iid(w, model.P_g, model.pi_eps)
    vnd_new = similar(vnd)
    l_allowed = l .<= model.lbar_nd
    if !any(l_allowed)
        l_allowed[l0_idx] = true
    end
    γ = model.gamma

    for gi in 1:Ng
        n_choice = n[:, :, gi]
        p_choice = pdefault[:, :, gi]
        wE_mat = wE[:, :, gi]
        allowed = (p_choice .<= model.pub) .& (n_choice .<= model.nbar_g[gi]) .& reshape(l_allowed, 1, Nl)
        if !any(allowed)
            allowed[b0_idx, l0_idx] = true
        end
        n_l_vec = (model.g[gi] .* l) ./ model.R_l_nd
        n_choice_plus = n_choice .+ reshape(n_l_vec, 1, Nl)
        wE_scaled = g_beta[gi] .* wE_mat

        Threads.@threads for tidx in 1:(Ne * Nl)
            ei = (tidx - 1) ÷ Nl + 1
            li = tidx - (ei - 1) * Nl
            y_val = y[gi, ei]
            l_cur = l[li]
            for bi in 1:Nb
                c_base = y_val - b[bi] - l_cur
                best_val = -Inf
                @inbounds for bj in 1:Nb
                    @inbounds for lprimej in 1:Nl
                        if allowed[bj, lprimej]
                            c = c_base + n_choice_plus[bj, lprimej]
                            val = u(c, γ) + wE_scaled[bj, lprimej]
                            if val > best_val
                                best_val = val
                            end
                        end
                    end
                end
                vnd_new[bi, li, gi, ei] = best_val
            end
        end
    end
    err = maximum(abs.(vnd_new .- vnd))
    vnd .= vnd_new
    return err
end

function _compute_policy_idx(model, b, l, y, g_beta, n, pdefault, vnd, vd, b0_idx, l0_idx)
    Nb, Nl, Ng, Ne = size(vnd)
    w = max.(vnd, vd)
    wE = expected_next_iid(w, model.P_g, model.pi_eps)
    b_policy_idx = zeros(Int, Nb, Nl, Ng, Ne)
    l_policy_idx = zeros(Int, Nb, Nl, Ng, Ne)
    l_allowed = l .<= model.lbar_nd
    if !any(l_allowed)
        l_allowed[l0_idx] = true
    end
    γ = model.gamma

    for gi in 1:Ng
        n_choice = n[:, :, gi]
        p_choice = pdefault[:, :, gi]
        wE_mat = wE[:, :, gi]
        allowed = (p_choice .<= model.pub) .& (n_choice .<= model.nbar_g[gi]) .& reshape(l_allowed, 1, Nl)
        if !any(allowed)
            allowed[b0_idx, l0_idx] = true
        end
        n_l_vec = (model.g[gi] .* l) ./ model.R_l_nd
        n_choice_plus = n_choice .+ reshape(n_l_vec, 1, Nl)
        wE_scaled = g_beta[gi] .* wE_mat
        Threads.@threads for tidx in 1:(Ne * Nl)
            ei = (tidx - 1) ÷ Nl + 1
            li = tidx - (ei - 1) * Nl
            y_val = y[gi, ei]
            l_cur = l[li]
            for bi in 1:Nb
                c_base = y_val - b[bi] - l_cur
                best_val = -Inf
                best_b = b0_idx
                best_l = l0_idx
                @inbounds for bj in 1:Nb
                    @inbounds for lprimej in 1:Nl
                        if allowed[bj, lprimej]
                            c = c_base + n_choice_plus[bj, lprimej]
                            val = u(c, γ) + wE_scaled[bj, lprimej]
                            if val > best_val
                                best_val = val
                                best_b = bj
                                best_l = lprimej
                            end
                        end
                    end
                end
                b_policy_idx[bi, li, gi, ei] = best_b
                l_policy_idx[bi, li, gi, ei] = best_l
            end
        end
    end
    return b_policy_idx, l_policy_idx
end

function solve_model(model::Model; sol::Union{Nothing, Solution} = nothing)
    b = model.b
    l = model.l
    g = model.g
    eps = model.eps

    Nb = model.Nb
    Nl = model.Nl
    Ng = model.Ng
    Ne = model.Ne
    @assert length(b) == Nb && length(l) == Nl && length(g) == Ng && length(eps) == Ne

    y = reshape(g, Ng, 1) .* exp.(model.sigma_eps .* reshape(eps, 1, Ne))
    g_beta = model.beta .* (g .^ (1 - model.gamma))

    bprime = b ./ reshape(g, 1, Ng)
    bprime_idx = nearest_index(b, bprime)
    kbprime_idx = nearest_index(b, model.kappa .* bprime)
    b0_idx = findmin(abs.(b))[2]
    l0_idx = findmin(abs.(l))[2]

    if sol === nothing || isempty(sol.vnd)
        base = @. u(y, model.gamma) / (1 - model.beta)
        vnd = repeat(reshape(base, 1, 1, Ng, Ne), Nb, Nl, 1, 1)
    else
        @assert size(sol.vnd) == (Nb, Nl, Ng, Ne)
        vnd = copy(sol.vnd)
    end

    if sol === nothing || isempty(sol.vd)
        vd = copy(vnd .- 1.0)
    else
        @assert size(sol.vd) == (Nb, Nl, Ng, Ne)
        vd = copy(sol.vd)
    end
    X = ones(Float64, Nb, Nl, Ng, Ne)

    d = falses(Nb, Nl, Ng, Ne)
    e = trues(Nb, Nl, Ng, Ne)

    outer_errs = Float64[]
    vd_errs = Float64[]
    x_errs = Float64[]
    vd_iters = Int[]
    x_iters = Int[]
    outer_times = Float64[]

    l_policy_idx_d = zeros(Int, Nb, Nl, Ng, Ne)

    for _ in 1:model.max_iter
        t0 = time()
        # Step 2: solve for default value and default LOLR policy
        vd_err, vd_iter, l_policy_idx_d = _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx, l, l0_idx)

        # Step 3: default and re-entry decisions
        d = vnd .< vd
        e = .!d

        # Step 4: solve prices Q and X
        Q, x_err, x_iter = _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx, l_policy_idx_d)

        # Step 5: schedule, issuance, and default probabilities
        _, R, n, pdefault = _compute_schedule(b, l, g, model, Q, d)

        # Step 6: update value in repayment
        err = _update_vnd!(vnd, vd, model, b, l, y, g_beta, n, pdefault, b0_idx, l0_idx)
        push!(vd_errs, vd_err)
        push!(x_errs, x_err)
        push!(vd_iters, vd_iter)
        push!(x_iters, x_iter)
        push!(outer_errs, err)
        push!(outer_times, time() - t0)
        print(".")
        flush(stdout)
        if err < model.tol
            break
        end
    end

    # Final policy evaluation
    Q = (1 .- d) .+ d .* X
    _, R, n, pdefault = _compute_schedule(b, l, g, model, Q, d)
    b_policy_idx, l_policy_idx = _compute_policy_idx(model, b, l, y, g_beta, n, pdefault, vnd, vd, b0_idx, l0_idx)

    outer_iters = length(outer_errs)
    return Solution(vnd, vd, d, e, Q, X, R, n, pdefault, b_policy_idx, l_policy_idx, l_policy_idx_d,
        outer_iters, vd_iters, x_iters, outer_errs, vd_errs, x_errs, outer_times)
end
