function _enforce_monotone(n::AbstractVector{<:Real}, idxs::AbstractVector{Int};
    tol::Float64 = 1e-10)
    mask = falses(length(n))
    last_n = -Inf
    for i in idxs
        if n[i] >= last_n - tol
            mask[i] = true
            last_n = n[i]
        end
    end
    return mask
end

function _branch_masks(n::AbstractVector{<:Real}; tol::Float64 = 1e-10)
    Nb = length(n)
    dn = diff(n)
    if all(dn .>= -tol) || all(dn .<= tol)
        return trues(Nb), trues(Nb)
    end

    imax = findmax(n)[2]
    low_idxs = collect(1:imax)
    high_idxs = collect(imax:Nb)
    low_mask = _enforce_monotone(n, low_idxs; tol = tol)
    high_mask = _enforce_monotone(n, reverse(high_idxs); tol = tol)
    return low_mask, high_mask
end

function _select_schedule_mask(n::Array{Float64,3})
    Nb, Ng, Ns = size(n)
    mask = falses(Nb, Ng, Ns)
    for gi in 1:Ng
        for si in 1:Ns
            low_mask, high_mask = _branch_masks(view(n, :, gi, si))
            # s=1 is bad (high-rate schedule), s=2 is good (low-rate schedule)
            mask[:, gi, si] .= (si == 1) ? high_mask : low_mask
        end
    end
    return mask
end

function _solve_vd!(vd, vnd, model, y, g_beta, imf_net, bprime_idx, kbprime_idx)
    Nb, Ng, Ns, Ne = size(vd)
    err = Inf
    iters = 0
    for _ in 1:model.max_iter_vd
        iters += 1
        vd_new = similar(vd)
        for gi in 1:Ng
            vnd_kbprime = vnd[kbprime_idx[:, gi], :, :, :]
            vd_bprime = vd[bprime_idx[:, gi], :, :, :]
            w_reentry = max.(vnd_kbprime, vd_bprime)
            cont = model.theta .* vd_bprime .+ (1 - model.theta) .* w_reentry
            contE = expected_next_iid(cont, model.P_g, model.P_s, model.pi_eps)
            udef = u(model.phi_g[gi] .* y[gi, :] .+ imf_net[gi, :], model.gamma)
            for si in 1:Ns
                contE_gs = reshape(contE[:, gi, si], Nb, 1)
                vd_new[:, gi, si, :] = reshape(udef, 1, Ne) .+ g_beta[gi] .* contE_gs
            end
        end
        err = maximum(abs.(vd_new .- vd))
        vd .= vd_new
        if err < model.tol_vd
            break
        end
    end
    return err, iters
end

function _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx)
    Nb, Ng, Ns, Ne = size(X)
    err = Inf
    iters = 0
    for _ in 1:model.max_iter_x
        iters += 1
        Q = (1 .- d) .+ d .* X
        X_new = similar(X)
        for gi in 1:Ng
            X_bprime = X[bprime_idx[:, gi], :, :, :]
            e_bprime = e[bprime_idx[:, gi], :, :, :]
            Q_kbprime = Q[kbprime_idx[:, gi], :, :, :]
            term = model.theta .* X_bprime .+ (1 - model.theta) .* ((1 .- e_bprime) .* X_bprime .+ e_bprime .* (model.kappa .* Q_kbprime))
            termE = expected_next_iid(term, model.P_g, model.P_s, model.pi_eps)
            for si in 1:Ns
                termE_gs = reshape(termE[:, gi, si], Nb, 1)
                X_new[:, gi, si, :] .= (1 / (1 + model.rstar)) .* termE_gs
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

function _compute_schedule(b, g, model, Q, d)
    Nb = length(b)
    Ng = length(g)
    QE = expected_next_iid(Q, model.P_g, model.P_s, model.pi_eps)
    R = (1 + model.rstar) ./ max.(QE, 1e-8)
    gb = reshape(b, Nb, 1, 1) .* reshape(g, 1, Ng, 1)
    n = gb ./ R
    pdefault = expected_next_iid(Float64.(d), model.P_g, model.P_s, model.pi_eps)
    return QE, R, n, pdefault
end

function _update_vnd!(vnd, vd, model, b, y, g_beta, imf_net, n, pdefault, schedule_mask, b0_idx)
    Nb, Ng, Ns, Ne = size(vnd)
    w = max.(vnd, vd)
    wE = expected_next_iid(w, model.P_g, model.P_s, model.pi_eps)
    vnd_new = similar(vnd)
    bcol = reshape(b, Nb, 1)
    for gi in 1:Ng
        for si in 1:Ns
            n_choice = n[:, gi, si]
            p_choice = pdefault[:, gi, si]
            allowed = schedule_mask[:, gi, si] .& (p_choice .<= model.pub) .& (n_choice .<= model.nbar_g[gi])
            if !any(allowed)
                allowed[b0_idx] = true
            end

            nrow = reshape(n_choice, 1, Nb)
            wE_row = reshape(wE[:, gi, si], 1, Nb)
            mask = reshape(allowed, 1, Nb)
            Threads.@threads for ei in 1:Ne
                c = y[gi, ei] + imf_net[gi, ei] .+ nrow .- bcol
                util = u(c, model.gamma)
                val = util .+ g_beta[gi] .* wE_row
                val_masked = ifelse.(mask, val, -Inf)
                vnd_new[:, gi, si, ei] = dropdims(maximum(val_masked, dims = 2), dims = 2)
            end
        end
    end
    err = maximum(abs.(vnd_new .- vnd))
    vnd .= vnd_new
    return err
end

function _compute_policy_idx(model, b, y, g_beta, imf_net, n, pdefault, schedule_mask, vnd, vd, b0_idx)
    Nb, Ng, Ns, Ne = size(vnd)
    w = max.(vnd, vd)
    wE = expected_next_iid(w, model.P_g, model.P_s, model.pi_eps)
    b_policy_idx = zeros(Int, Nb, Ng, Ns, Ne)
    bcol = reshape(b, Nb, 1)
    for gi in 1:Ng
        for si in 1:Ns
            n_choice = n[:, gi, si]
            p_choice = pdefault[:, gi, si]
            allowed = schedule_mask[:, gi, si] .& (p_choice .<= model.pub) .& (n_choice .<= model.nbar_g[gi])
            if !any(allowed)
                allowed[b0_idx] = true
            end

            nrow = reshape(n_choice, 1, Nb)
            wE_row = reshape(wE[:, gi, si], 1, Nb)
            mask = reshape(allowed, 1, Nb)
            Threads.@threads for ei in 1:Ne
                c = (y[gi, ei] + imf_net[gi, ei]) .+ nrow .- bcol
                util = u(c, model.gamma)
                val = util .+ g_beta[gi] .* wE_row
                val_masked = ifelse.(mask, val, -Inf)

                idx = findmax(val_masked, dims = 2)[2]
                b_policy_idx[:, gi, si, ei] = map(i -> i[2], idx)
            end
        end
    end
    return b_policy_idx
end

function solve_model(model::Model; sol::Union{Nothing, Solution} = nothing)
    b = model.b
    g = model.g
    eps = model.eps

    Nb = model.Nb
    Ng = model.Ng
    Ns = model.Ns
    Ne = model.Ne
    @assert length(b) == Nb && length(g) == Ng && length(eps) == Ne

    imf_net = (1 - model.R_l) .* model.l_policy
    y = reshape(g, Ng, 1) .* exp.(model.sigma_eps .* reshape(eps, 1, Ne))
    g_beta = model.beta .* (g .^ (1 - model.gamma))

    bprime = b ./ reshape(g, 1, Ng)
    bprime_idx = nearest_index(b, bprime)
    kbprime_idx = nearest_index(b, model.kappa .* bprime)
    b0_idx = findmin(abs.(b))[2]

    if sol === nothing || isempty(sol.vnd)
        vnd_base = reshape(u(y .+ imf_net, model.gamma) ./ (1 - model.beta), 1, Ng, 1, Ne)
        vnd = repeat(vnd_base, Nb, 1, Ns, 1)
    else
        @assert size(sol.vnd) == (Nb, Ng, Ns, Ne)
        vnd = copy(sol.vnd)
    end

    if sol === nothing || isempty(sol.vd)
        vd = copy(vnd .- 1.0)
    else
        @assert size(sol.vd) == (Nb, Ng, Ns, Ne)
        vd = copy(sol.vd)
    end
    X = ones(Float64, Nb, Ng, Ns, Ne)

    d = falses(Nb, Ng, Ns, Ne)
    e = trues(Nb, Ng, Ns, Ne)

    outer_errs = Float64[]
    vd_errs = Float64[]
    x_errs = Float64[]
    vd_iters = Int[]
    x_iters = Int[]
    outer_times = Float64[]

    for _ in 1:model.max_iter
        t0 = time()
        # Step 2: solve for default value
        vd_err, vd_iter = _solve_vd!(vd, vnd, model, y, g_beta, imf_net, bprime_idx, kbprime_idx)

        # Step 3: default and re-entry decisions
        d = vnd .< vd
        e = .!d

        # Step 4: solve prices Q and X
        Q, x_err, x_iter = _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx)

        # Step 5: schedule, issuance, and default probabilities
        _, R, n, pdefault = _compute_schedule(b, g, model, Q, d)
        schedule_mask = _select_schedule_mask(n)

        # Step 6: update value in repayment
        err = _update_vnd!(vnd, vd, model, b, y, g_beta, imf_net, n, pdefault, schedule_mask, b0_idx)
        push!(vd_errs, vd_err)
        push!(x_errs, x_err)
        push!(vd_iters, vd_iter)
        push!(x_iters, x_iter)
        push!(outer_errs, err)
        push!(outer_times, time() - t0)
        if err < model.tol
            break
        end
    end

    # Final policy evaluation
    Q = (1 .- d) .+ d .* X
    _, R, n, pdefault = _compute_schedule(b, g, model, Q, d)
    schedule_mask = _select_schedule_mask(n)
    b_policy_idx = _compute_policy_idx(model, b, y, g_beta, imf_net, n, pdefault, schedule_mask, vnd, vd, b0_idx)

    outer_iters = length(outer_errs)
    return Solution(vnd, vd, d, e, Q, X, R, n, schedule_mask, pdefault, b_policy_idx,
        outer_iters, vd_iters, x_iters, outer_errs, vd_errs, x_errs, outer_times)
end
