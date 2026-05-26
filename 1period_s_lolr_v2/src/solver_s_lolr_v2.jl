# Sunspot schedule selection

function _increasing_part_mask(n::AbstractVector{<:Real}; window::Int = 2, tol::Float64 = 1e-10)
    N = length(n)
    mask = falses(N)
    if N == 0
        return mask
    end

    if N <= 2 * window
        mask .= true
        return mask
    end

    # Match solver_s: mark the "increasing part" using rolling averages
    # of issuance on each side of a point.
    for i in (window + 1):(N - window)
        left_avg = sum(@view n[(i - window):i]) / (window + 1)
        right_avg = sum(@view n[i:(i + window)]) / (window + 1)
        if left_avg < right_avg - tol
            mask[i] = true
        end
    end

    mask[1:window] .= mask[window + 1]
    mask[(N - (window - 1)):N] .= mask[N - window]
    return mask
end

function _bad_branch_mask(n::AbstractVector{<:Real}, pdefault::AbstractVector{<:Real}, pub::Float64;
    window::Int = 2, tol::Float64 = 1e-10)
    N = length(n)
    @assert length(pdefault) == N

    # Match solver_s:
    # lsch: below the default-probability cutoff
    # isch: locally increasing part of the issuance correspondence
    # msch: multiplicity filter that removes points dominated by a later
    #       feasible point with lower issuance
    lsch = pdefault .< (0.999 * pub)
    isch = _increasing_part_mask(n; window = window, tol = tol)
    lisch = lsch .& isch
    msch = trues(N)

    for i in 1:N
        min_future_n = Inf
        for j in i:N
            if lisch[j] && n[j] < min_future_n
                min_future_n = n[j]
            end
        end
        if min_future_n < n[i] - tol
            msch[i] = false
        end
    end

    return msch
end

function _select_schedule_mask(n::Array{Float64,4}, R::Array{Float64,4}, pdefault::Array{Float64,4}, model)
    Nb, Nl, Ng, Ns = size(n)
    mask = falses(Nb, Nl, Ng, Ns)
    for gi in 1:Ng
        for si in 1:Ns
            if si == 1
                # Bad sunspot: for each fixed l', apply the solver_s branch
                # selector over the private-debt schedule b' only.
                @inbounds for lprimej in 1:Nl
                    mask[:, lprimej, gi, si] .= _bad_branch_mask(
                        view(n, :, lprimej, gi, si),
                        view(pdefault, :, lprimej, gi, si),
                        model.pub,
                    )
                end
            else
                # Good sunspot: everything is accessible before applying the
                # repayment constraints in the choice problem.
                mask[:, :, gi, si] .= true
            end
        end
    end
    return mask
end

# Default value

function _default_continuation_by_lprime!(contE_g, vd, vnd, model, gi, bidx, kbidx, w_reentry, cont, contE_tmp, contE_tmp_vec)
    Nl = size(vd, 2)
    Ns = size(vd, 4)
    @inbounds for lprimej in 1:Nl
        @views vnd_kb = vnd[kbidx, lprimej, :, :, :]
        @views vd_b = vd[bidx, lprimej, :, :, :]
        @. w_reentry = max(vnd_kb, vd_b)
        @. cont = model.theta * vd_b + (1 - model.theta) * w_reentry
        expected_next!(contE_tmp, cont, model.K, model.pi_eps, contE_tmp_vec)
        @views contE_g[:, lprimej, :] .= contE_tmp[:, gi, :]
    end
    return contE_g
end

function _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx)
    Nb, Nl, Ng, Ns, Ne = size(vd)
    l = model.l
    l_prime = model.l
    γ = model.gamma
    l0_idx = argmin(abs.(l))
    l_policy_idx_d = fill(l0_idx, Nb, Nl, Ng, Ns, Ne)
    err = Inf
    w_reentry = Array{Float64}(undef, Nb, Ng, Ns, Ne)
    cont = Array{Float64}(undef, Nb, Ng, Ns, Ne)
    contE_tmp = Array{Float64}(undef, Nb, Ng, Ns)
    contE_tmp_vec = Vector{Float64}(undef, Nb * Ng * Ns)

    for _ in 1:model.max_iter_vd
        vd_new = similar(vd)
        for gi in 1:Ng
            bidx = bprime_idx[:, gi]
            kbidx = kbprime_idx[:, gi]
            contE_g = Array{Float64}(undef, Nb, Nl, Ns)
            _default_continuation_by_lprime!(contE_g, vd, vnd, model, gi, bidx, kbidx, w_reentry, cont, contE_tmp, contE_tmp_vec)
            n_l_vec = (model.g[gi] .* l_prime) ./ model.R_l

            for si in 1:Ns
                l_last = max(searchsortedlast(l_prime, model.lbar_dgs[2, gi, si]), l0_idx)
                contE_mat = @view contE_g[:, :, si]
                g_contE = g_beta[gi] .* contE_mat[:, 1:l_last]
                for ei in 1:Ne, li in 1:Nl
                    y_def = model.phi_g[gi] * y[gi, ei]
                    util_row = reshape(u.(y_def .+ n_l_vec[1:l_last] .- l[li], γ), 1, l_last)
                    val_mat = g_contE .+ util_row
                    best_val, best_idx = findmax(val_mat, dims = 2)
                    vd_new[:, li, gi, si, ei] .= vec(best_val)
                    l_policy_idx_d[:, li, gi, si, ei] .= getindex.(best_idx, 2)
                end
            end
        end
        err = max_abs_diff_safe(vd_new, vd)
        vd .= vd_new
        err < model.tol_vd && break
    end

    return err, l_policy_idx_d
end

# Price recursion

function _gather_policy_slice!(out, v, bidx, lidx, gp, sp, ep)
    Nb, Nl = size(out)
    @inbounds for li in 1:Nl, bi in 1:Nb
        out[bi, li] = v[bidx[bi], lidx[bi, li], gp, sp, ep]
    end
    return out
end

function _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx, l_policy_idx_d)
    Nb, Nl, Ng, Ns, Ne = size(X)
    err = Inf
    discount = 1 / (1 + model.rstar)
    x_next = Array{Float64}(undef, Nb, Nl)
    q_next = Array{Float64}(undef, Nb, Nl)
    e_next = Matrix{Bool}(undef, Nb, Nl)
    term = Array{Float64}(undef, Nb, Nl)
    termE = Array{Float64}(undef, Nb, Nl)

    for _ in 1:model.max_iter_x
        Q = (1 .- d) .+ d .* X
        X_new = similar(X)
        for gi in 1:Ng
            bidx = bprime_idx[:, gi]
            kbidx = kbprime_idx[:, gi]
            for si in 1:Ns, ei in 1:Ne
                lidx = @view l_policy_idx_d[:, :, gi, si, ei]
                fill!(termE, 0.0)
                for gp in 1:Ng, sp in 1:Ns, ep in 1:Ne
                    weight = model.P_g[gi, gp] * model.P_s[si, sp] * model.pi_eps[ep]
                    _gather_policy_slice!(x_next, X, bidx, lidx, gp, sp, ep)
                    _gather_policy_slice!(e_next, e, bidx, lidx, gp, sp, ep)
                    _gather_policy_slice!(q_next, Q, kbidx, lidx, gp, sp, ep)
                    @. term = model.theta * x_next + (1 - model.theta) * ifelse(e_next, model.kappa * q_next, x_next)
                    @. termE += weight * term
                end
                @views X_new[:, :, gi, si, ei] .= discount .* termE
            end
        end
        err = maximum(abs.(X_new .- X))
        X .= X_new
        err < model.tol_x && break
    end

    return (1 .- d) .+ d .* X, err
end

# Repayment schedule

function _compute_schedule!(QE, R, n, pdefault, model::Model, Q, d, gb, QE_tmp, pdefault_tmp)
    expected_next!(QE, Q, model.K, model.pi_eps, QE_tmp)
    @. R = (1 + model.rstar) / max(QE, 1e-8)
    @. n = gb / R
    expected_next!(pdefault, d, model.K, model.pi_eps, pdefault_tmp)
    return nothing
end

function _repayment_choice_mask(model::Model, n_choice, p_choice, schedule_mask_gs, gi, l_last, b0_idx, l0_idx)
    allowed = falses(model.Nb, l_last)
    @inbounds for lprimej in 1:l_last
        allowed[:, lprimej] .= view(schedule_mask_gs, :, lprimej) .&
            (view(p_choice, :, lprimej) .<= model.pub) .&
            (view(n_choice, :, lprimej) .<= model.nbar_g[gi])
    end
    if !any(allowed)
        allowed[b0_idx, l0_idx] = true
    end
    return allowed
end

# Repayment value and policies

function _update_vnd!(vnd, vd, model, y, g_beta, n, pdefault, schedule_mask, b0_idx, l0_idx; damp::Float64 = 0.5)
    Nb, Nl, Ng, Ns, Ne = size(vnd)
    b = model.b
    l = model.l
    l_prime = model.l
    w = similar(vnd)
    @. w = max(vnd, vd)
    wE = Array{Float64}(undef, Nb, Nl, Ng, Ns)
    wE_tmp = Vector{Float64}(undef, Nb * Nl * Ng * Ns)
    expected_next!(wE, w, model.K, model.pi_eps, wE_tmp)
    vnd_new = similar(vnd)
    γ = model.gamma

    for gi in 1:Ng, si in 1:Ns
        l_last = max(searchsortedlast(l_prime, model.lbar_dgs[1, gi, si]), l0_idx)
        @views n_choice = n[:, :, gi, si]
        @views p_choice = pdefault[:, :, gi, si]
        @views wE_gs = wE[:, :, gi, si]
        @views schedule_mask_gs = schedule_mask[:, :, gi, si]
        allowed = _repayment_choice_mask(model, n_choice, p_choice, schedule_mask_gs, gi, l_last, b0_idx, l0_idx)
        n_l_vec = (model.g[gi] .* l_prime[1:l_last]) ./ model.R_l

        Threads.@threads for tidx in 1:(Ne * Nl)
            ei = (tidx - 1) ÷ Nl + 1
            li = tidx - (ei - 1) * Nl
            y_cur = y[gi, ei]
            l_cur = l[li]
            best = fill(-Inf, Nb)
            val_vec = Vector{Float64}(undef, Nb)

            for lprimej in 1:l_last
                @views resource_vec = n_choice[:, lprimej]
                @views continuation_vec = wE_gs[:, lprimej]
                @views mask_vec = allowed[:, lprimej]
                n_l_scalar = n_l_vec[lprimej]
                g_beta_g = g_beta[gi]

                for bi in 1:Nb
                    base = y_cur - l_cur - b[bi]
                    @. val_vec = ifelse(mask_vec, u(base + resource_vec + n_l_scalar, γ) + g_beta_g * continuation_vec, -Inf)
                    best[bi] = max(best[bi], maximum(val_vec))
                end
            end

            bad = .!isfinite.(best)
            if any(bad)
                # If repayment is impossible, keep the value finite but strictly
                # below default so the state defaults without breaking iteration.
                vd_cur = view(vd, :, li, gi, si, ei)
                best[bad] .= vd_cur[bad] .- 1.0
            end
            vnd_new[:, li, gi, si, ei] = best
        end
    end

    err = max_abs_diff_safe(vnd_new, vnd)
    damp_update_safe!(vnd, vnd, vnd_new, damp)
    return err
end

function _compute_policy_idx(model::Model, vnd, vd, y, g_beta, n, pdefault, schedule_mask, b0_idx, l0_idx)
    Nb, Nl, Ng, Ns, Ne = size(vnd)
    b = model.b
    l = model.l
    l_prime = model.l
    w = similar(vnd)
    @. w = max(vnd, vd)
    wE = Array{Float64}(undef, Nb, Nl, Ng, Ns)
    wE_tmp = Vector{Float64}(undef, Nb * Nl * Ng * Ns)
    expected_next!(wE, w, model.K, model.pi_eps, wE_tmp)
    b_policy_idx = fill(b0_idx, Nb, Nl, Ng, Ns, Ne)
    l_policy_idx = fill(l0_idx, Nb, Nl, Ng, Ns, Ne)
    γ = model.gamma

    for gi in 1:Ng, si in 1:Ns
        l_last = max(searchsortedlast(l_prime, model.lbar_dgs[1, gi, si]), l0_idx)
        @views n_choice = n[:, :, gi, si]
        @views p_choice = pdefault[:, :, gi, si]
        @views wE_gs = wE[:, :, gi, si]
        @views schedule_mask_gs = schedule_mask[:, :, gi, si]
        allowed = _repayment_choice_mask(model, n_choice, p_choice, schedule_mask_gs, gi, l_last, b0_idx, l0_idx)
        n_l_vec = (model.g[gi] .* l_prime[1:l_last]) ./ model.R_l

        Threads.@threads for tidx in 1:(Ne * Nl)
            ei = (tidx - 1) ÷ Nl + 1
            li = tidx - (ei - 1) * Nl
            y_cur = y[gi, ei]
            l_cur = l[li]
            best_val = fill(-Inf, Nb)
            best_b = fill(b0_idx, Nb)
            best_l = fill(l0_idx, Nb)
            val_vec = Vector{Float64}(undef, Nb)

            for lprimej in 1:l_last
                @views resource_vec = n_choice[:, lprimej]
                @views continuation_vec = wE_gs[:, lprimej]
                @views mask_vec = allowed[:, lprimej]
                n_l_scalar = n_l_vec[lprimej]
                g_beta_g = g_beta[gi]

                for bi in 1:Nb
                    base = y_cur - l_cur - b[bi]
                    @. val_vec = ifelse(mask_vec, u(base + resource_vec + n_l_scalar, γ) + g_beta_g * continuation_vec, -Inf)
                    row_best, row_b = findmax(val_vec)
                    if row_best > best_val[bi]
                        best_val[bi] = row_best
                        best_b[bi] = row_b
                        best_l[bi] = lprimej
                    end
                end
            end

            b_policy_idx[:, li, gi, si, ei] = best_b
            l_policy_idx[:, li, gi, si, ei] = best_l
        end
    end

    return b_policy_idx, l_policy_idx
end

# Main solver

function solve_model(model::Model; sol::Union{Nothing, Solution} = nothing, verbose::Bool = true)
    sol === nothing && (sol = initial_solution(model))

    l = model.l
    y = reshape(model.g, model.Ng, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, model.Ne))
    g_beta = model.beta .* (model.g .^ (1 - model.gamma))
    gb = reshape(model.b, model.Nb, 1, 1, 1) .* reshape(model.g, 1, 1, model.Ng, 1)
    bprime_idx = nearest_index(model.b, model.b ./ reshape(model.g, 1, model.Ng))
    kbprime_idx = nearest_index(model.b, model.kappa .* model.b ./ reshape(model.g, 1, model.Ng))
    b0_idx = findmin(abs.(model.b))[2]
    l0_idx = findmin(abs.(l))[2]

    vnd = copy(sol.vnd)
    vd = copy(sol.vd)
    X = copy(sol.X)
    d = copy(sol.d)
    e = copy(sol.e)
    QE = Array{Float64}(undef, model.Nb, model.Nl, model.Ng, model.Ns)
    R = copy(sol.R)
    n = copy(sol.n)
    pdefault = copy(sol.pdefault)
    QE_tmp = Vector{Float64}(undef, model.Nb * model.Nl * model.Ng * model.Ns)
    pdefault_tmp = Vector{Float64}(undef, model.Nb * model.Nl * model.Ng * model.Ns)

    outer_errs = Float64[]
    vd_errs = Float64[]
    x_errs = Float64[]
    l_policy_idx_d = copy(sol.l_policy_idx_d)

    for it in 1:model.max_iter
        damp = max(0.05, 0.9 - 0.0 * ((it - 1) / 10))
        # Step 2: solve for default value
        vd_err, l_policy_idx_d = _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx)

        # Step 3: default and re-entry decisions
        d = vnd .< vd
        e = .!d

        # Step 4: solve prices Q and X
        Q, x_err = _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx, l_policy_idx_d)

        # Step 5: schedule, issuance, and default probabilities
        _compute_schedule!(QE, R, n, pdefault, model, Q, d, gb, QE_tmp, pdefault_tmp)
        schedule_mask = _select_schedule_mask(n, R, pdefault, model)

        # Step 6: update value in repayment
        err = _update_vnd!(vnd, vd, model, y, g_beta, n, pdefault, schedule_mask, b0_idx, l0_idx; damp = damp)

        push!(outer_errs, err)
        push!(vd_errs, vd_err)
        push!(x_errs, x_err)

        if verbose
            println("iter=", it, ", vnd_err=", err, ", vd_err=", vd_err, ", x_err=", x_err, ", damp=", damp)
        end
        err < model.tol && break
    end

    d = vnd .< vd
    e = .!d
    Q = (1 .- d) .+ d .* X
    _compute_schedule!(QE, R, n, pdefault, model, Q, d, gb, QE_tmp, pdefault_tmp)
    schedule_mask = _select_schedule_mask(n, R, pdefault, model)
    b_policy_idx, l_policy_idx = _compute_policy_idx(model, vnd, vd, y, g_beta, n, pdefault, schedule_mask, b0_idx, l0_idx)

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
        outer_errs,
        vd_errs,
        x_errs,
    )
end
