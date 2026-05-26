# Schedule selection in the bad sunspot follows the existing sunspot model.

# Identify locally increasing points on the issuance correspondence.
function _increasing_part_mask(n::AbstractVector{<:Real}; window::Int = 2, tol::Float64 = 1e-10)
    Nb = length(n)
    mask = falses(Nb)
    Nb == 0 && return mask
    if Nb <= 2 * window
        mask .= true
        return mask
    end

    for i in (window + 1):(Nb - window)
        left_avg = sum(@view n[(i - window):i]) / (window + 1)
        right_avg = sum(@view n[i:(i + window)]) / (window + 1)
        mask[i] = left_avg < right_avg - tol
    end
    mask[1:window] .= mask[window + 1]
    mask[(Nb - (window - 1)):Nb] .= mask[Nb - window]
    return mask
end

# Select the high-rate branch that remains available in the bad sunspot.
function _bad_branch_mask(n::AbstractVector{<:Real}, pdefault::AbstractVector{<:Real}, pub::Real)
    feasible = pdefault .< (0.999 * pub)
    increasing = _increasing_part_mask(n)
    candidate = feasible .& increasing
    mask = trues(length(n))

    for i in eachindex(n)
        min_future_n = Inf
        for j in i:length(n)
            if candidate[j]
                min_future_n = min(min_future_n, n[j])
            end
        end
        mask[i] = !(min_future_n < n[i] - 1e-10)
    end
    return mask
end

# Restrict the schedule in bad sunspots and leave the good-sunspot schedule open.
function _select_schedule_mask(n::Array{Float64,3}, pdefault::Array{Float64,3}, model::Model)
    mask = trues(size(n))
    for gi in 1:model.Ng
        mask[:, gi, 1] .= _bad_branch_mask(view(n, :, gi, 1), view(pdefault, :, gi, 1), model.pub)
    end
    return mask
end

# Iterate on the value of default, equation (12), taking repayment value as given.
function _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx)
    err = Inf
    for _ in 1:model.max_iter_vd
        vd_new = similar(vd)
        for gi in 1:model.Ng
            vd_next = vd[bprime_idx[:, gi], :, :, :]
            repay_after_haircut = vnd[kbprime_idx[:, gi], :, :, :]
            continuation = model.theta .* vd_next .+
                (1 - model.theta) .* max.(repay_after_haircut, vd_next)
            continuation_E = expected_next(continuation, model)
            utility_default = u(model.phi_g[gi] .* y[gi, :], model.gamma)
            for si in 1:model.Ns
                vd_new[:, gi, si, :] .= reshape(utility_default, 1, model.Ne) .+
                    g_beta[gi] .* reshape(continuation_E[:, gi, si], model.Nb, 1)
            end
        end
        err = max_abs_diff_safe(vd_new, vd)
        vd .= vd_new
        err < model.tol_vd && break
    end
    return err
end

# Iterate on recovery value X and construct bond price Q, equations (17)-(18).
function _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx)
    err = Inf
    for _ in 1:model.max_iter_x
        Q = (1 .- d) .+ d .* X
        X_new = similar(X)
        for gi in 1:model.Ng
            X_next = X[bprime_idx[:, gi], :, :, :]
            e_next = e[bprime_idx[:, gi], :, :, :]
            Q_haircut = Q[kbprime_idx[:, gi], :, :, :]
            recovery = (1 .- e_next) .* X_next .+ e_next .* model.kappa .* Q_haircut
            continuation = model.theta .* X_next .+ (1 - model.theta) .* recovery
            continuation_E = expected_next(continuation, model)
            for si in 1:model.Ns
                X_new[:, gi, si, :] .= reshape(continuation_E[:, gi, si], model.Nb, 1) ./
                    (1 + model.rstar)
            end
        end
        err = max_abs_diff_safe(X_new, X)
        X .= X_new
        err < model.tol_x && break
    end
    return (1 .- d) .+ d .* X, err
end

# Construct the lender schedule with ECB intervention, equations (13)-(16).
# H is the payoff from holding; L is the payoff from queueing at the facility.
function _compute_schedule(model::Model, Q, d)
    QE = expected_next(Q, model)
    H = QE ./ (1 + model.rstar)
    L = copy(H)

    # A bondholder queues at the facility only for positive new sovereign debt.
    for gi in 1:model.Ng, si in 1:model.Ns
        prob = ecb_sale_probability(model, gi, si)
        price = ecb_purchase_price(model, gi, si)
        for bi in 1:model.Nb
            if model.b[bi] > 0 && prob > 0
                L[bi, gi, si] = prob * price + (1 - prob) * H[bi, gi, si]
            end
        end
    end

    lenders_value = max.(H, L)
    R = 1.0 ./ max.(lenders_value, 1e-8)
    gb = reshape(model.b, model.Nb, 1, 1) .* reshape(model.g, 1, model.Ng, 1)
    n = gb ./ R
    pdefault = expected_next(Float64.(d), model)
    ecb_used = L .> H .+ 1e-12
    return H, L, R, n, pdefault, ecb_used
end

# Enforce risk, issuance-cap, and sunspot-branch restrictions on debt choices.
function _allowed_choices(model::Model, n_choice, p_choice, schedule_choice, gi::Int, b0_idx::Int)
    allowed = schedule_choice .& (p_choice .<= model.pub) .& (n_choice .<= model.nbar_g[gi])
    any(allowed) || (allowed[b0_idx] = true)
    return allowed
end

# Update repayment value by maximizing over permitted choices, equations (21)-(22).
function _update_vnd!(vnd, vd, model, y, g_beta, n, pdefault, schedule_mask, b0_idx;
    damp::Float64 = 0.5)
    wE = expected_next(max.(vnd, vd), model)
    vnd_new = similar(vnd)
    b_current = reshape(model.b, model.Nb, 1)

    for gi in 1:model.Ng, si in 1:model.Ns
        n_choice = view(n, :, gi, si)
        allowed = _allowed_choices(model, n_choice, view(pdefault, :, gi, si),
            view(schedule_mask, :, gi, si), gi, b0_idx)
        continuation = reshape(view(wE, :, gi, si), 1, model.Nb)
        issuance = reshape(n_choice, 1, model.Nb)
        choice_mask = reshape(allowed, 1, model.Nb)

        Threads.@threads for ei in 1:model.Ne
            c = y[gi, ei] .+ issuance .- b_current
            values = u(c, model.gamma) .+ g_beta[gi] .* continuation
            values = ifelse.(choice_mask, values, -Inf)
            vnd_new[:, gi, si, ei] .= vec(maximum(values, dims = 2))
        end
    end

    err = max_abs_diff_safe(vnd_new, vnd)
    damp_update_safe!(vnd, vnd_new, damp)
    return err
end

# Recover the optimal next-period private-debt index at each repayment state.
function _compute_policy_idx(model, vnd, vd, y, g_beta, n, pdefault, schedule_mask, b0_idx)
    wE = expected_next(max.(vnd, vd), model)
    policy = fill(b0_idx, model.Nb, model.Ng, model.Ns, model.Ne)
    b_current = reshape(model.b, model.Nb, 1)

    for gi in 1:model.Ng, si in 1:model.Ns
        n_choice = view(n, :, gi, si)
        allowed = _allowed_choices(model, n_choice, view(pdefault, :, gi, si),
            view(schedule_mask, :, gi, si), gi, b0_idx)
        continuation = reshape(view(wE, :, gi, si), 1, model.Nb)
        issuance = reshape(n_choice, 1, model.Nb)
        choice_mask = reshape(allowed, 1, model.Nb)

        Threads.@threads for ei in 1:model.Ne
            c = y[gi, ei] .+ issuance .- b_current
            values = u(c, model.gamma) .+ g_beta[gi] .* continuation
            values = ifelse.(choice_mask, values, -Inf)
            idx = findmax(values, dims = 2)[2]
            policy[:, gi, si, ei] .= map(x -> x[2], idx)
        end
    end
    return policy
end

# Solve the ECB fixed point using the algorithm in Section 2.2 of the note.
function solve_model(model::Model; sol::Union{Nothing,Solution} = nothing, verbose::Bool = true)
    # Step 1: start from a guess for repayment value vnd.
    sol === nothing && (sol = initial_solution(model))

    # Step 0 is performed by init_model; prepare flow output and debt mappings here.
    y = reshape(model.g, model.Ng, 1) .* exp.(model.sigma_eps .* reshape(model.eps, 1, model.Ne))
    g_beta = model.beta .* model.g .^ (1 - model.gamma)
    bprime_idx = nearest_index(model.b, model.b ./ reshape(model.g, 1, model.Ng))
    kbprime_idx = nearest_index(model.b, model.kappa .* model.b ./ reshape(model.g, 1, model.Ng))
    b0_idx = argmin(abs.(model.b))

    vnd = copy(sol.vnd)
    vd = copy(sol.vd)
    X = copy(sol.X)
    d = copy(sol.d)
    e = copy(sol.e)
    outer_errs = Float64[]
    vd_errs = Float64[]
    x_errs = Float64[]

    for it in 1:model.max_iter
        # Step 2: given vnd, solve the fixed point for default value vd.
        vd_err = _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx)

        # Step 3: determine default and re-entry, equations (19)-(20).
        d = vnd .< vd
        e = .!d

        # Step 4: solve the bond-price and recovery recursions, equations (17)-(18).
        Q, x_err = _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx)

        # Step 5: build H, L, R, issuance, and expected default probabilities.
        H, L, R, n, pdefault, ecb_used = _compute_schedule(model, Q, d)

        # Step 6: remove debt choices unavailable under the realized sunspot.
        schedule_mask = _select_schedule_mask(n, pdefault, model)

        # Step 7: update repayment value using the available debt choices.
        err = _update_vnd!(vnd, vd, model, y, g_beta, n, pdefault, schedule_mask, b0_idx)

        push!(outer_errs, err)
        push!(vd_errs, vd_err)
        push!(x_errs, x_err)
        verbose && println("iter=", it, ", vnd_err=", err, ", vd_err=", vd_err, ", x_err=", x_err)
        err < model.tol && break
    end

    # Recompute equilibrium objects so they correspond to the final vnd iterate.
    _solve_vd!(vd, vnd, model, y, g_beta, bprime_idx, kbprime_idx)
    d = vnd .< vd
    e = .!d
    Q, _ = _solve_prices!(X, d, e, model, bprime_idx, kbprime_idx)
    H, L, R, n, pdefault, ecb_used = _compute_schedule(model, Q, d)
    schedule_mask = _select_schedule_mask(n, pdefault, model)
    b_policy_idx = _compute_policy_idx(model, vnd, vd, y, g_beta, n, pdefault,
        schedule_mask, b0_idx)

    return Solution(vnd, vd, d, e, Q, X, H, L, R, n, pdefault, schedule_mask,
        ecb_used, b_policy_idx, outer_errs, vd_errs, x_errs)
end
