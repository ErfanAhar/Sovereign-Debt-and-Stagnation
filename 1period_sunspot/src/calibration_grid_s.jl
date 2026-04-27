using JLD2
using Printf

# Clear filename so each saved file tells you the parameterization directly.
case_filename(phi_bad, phi_good, beta, ps_bad) = @sprintf(
    "sol_phi_bad=%0.4f_phi_good=%0.4f_beta=%0.4f_ps_bad=%0.2f.jld2",
    phi_bad, phi_good, beta, ps_bad,
)

# Distance on the calibration grid for nearest-neighbor warm starts.
case_distance(p1, p2) =
    20 * abs(p1.phi_bad - p2.phi_bad) +
    20 * abs(p1.phi_good - p2.phi_good) +
    (5 / 0.45) * abs(p1.beta - p2.beta) +
    (1 / 0.24) * abs(p1.ps_bad - p2.ps_bad)

# The raw model is mutable, so each loop only changes the calibration targets.
function set_params!(raw, phi_bad, phi_good, beta, ps_bad)
    raw.beta = beta
    raw.phi_g .= [phi_bad, phi_good]
    raw.P_s .= [ps_bad 1 - ps_bad; ps_bad 1 - ps_bad]
end

# For warm starts we only need the two value functions.
function guess_solution(solution_data, model)
    sz4 = (model.Nb, model.Ng, model.Ns, model.Ne)
    sz3 = (model.Nb, model.Ng, model.Ns)
    vnd = Float64.(solution_data.vnd)
    vd = Float64.(solution_data.vd)
    d = vnd .< vd
    Q = hasproperty(solution_data, :Q) ? Float64.(solution_data.Q) : ones(Float64, sz4)
    X = ones(Float64, sz4)
    X[d] .= Q[d]
    return Solution(
        vnd,
        vd,
        d,
        .!d,
        Q,
        X,
        ones(Float64, sz3),
        zeros(Float64, sz3),
        falses(sz3...),
        zeros(Float64, sz3),
        zeros(Int, sz4),
        0,
        Int[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
    )
end

# Save plain arrays and scalars so the file is small and easy to reload.
function save_case(path, model, sol, phi_bad, phi_good, beta, ps_bad, status, converged, guess_source, error_message)
    solution_data = sol === nothing ? nothing : (
        vnd = Float32.(sol.vnd),
        vd = Float32.(sol.vd),
        Q = Float32.(sol.Q),
        R = Float32.(sol.R),
        n = Float32.(sol.n),
        pdefault = Float32.(sol.pdefault),
        schedule_mask = BitArray(sol.schedule_mask),
        b_policy_idx = Int32.(sol.b_policy_idx),
        outer_iters = sol.outer_iters,
        final_outer_err = sol.outer_errs[end],
    )

    model_data = (
        Nb = model.Nb,
        Ng = model.Ng,
        Ns = model.Ns,
        Ne = model.Ne,
        b = Float32.(model.b),
        g = Float32.(model.g),
        eps = Float32.(model.eps),
        P_g = Float32.(model.P_g),
        pi_eps = Float32.(model.pi_eps),
        P_s = Float32.(model.P_s),
        K = Float32.(model.K),
        beta = Float32(model.beta),
        gamma = Float32(model.gamma),
        sigma_eps = Float32(model.sigma_eps),
        rstar = Float32(model.rstar),
        theta = Float32(model.theta),
        kappa = Float32(model.kappa),
        phi_g = Float32.(model.phi_g),
        nbar_g = Float32.(model.nbar_g),
        R_l = Float32(model.R_l),
        l_policy = Float32.(model.l_policy),
        pub = Float32(model.pub),
    )

    JLD2.jldopen(path, "w"; compress = true) do f
        f["phi_bad"] = phi_bad
        f["phi_good"] = phi_good
        f["beta"] = beta
        f["ps_bad"] = ps_bad
        f["status"] = status
        f["converged"] = converged
        f["guess_source"] = guess_source
        f["error_message"] = error_message
        f["model_data"] = model_data
        f["solution_data"] = solution_data
    end
end

# Reload one saved case later for plotting.
function load_case(path)
    file = JLD2.load(path)
    m = file["model_data"]
    model = Model(
        Nb = m.Nb,
        Ng = m.Ng,
        Ns = m.Ns,
        Ne = m.Ne,
        b = Float64.(m.b),
        g = Float64.(m.g),
        eps = Float64.(m.eps),
        P_g = Float64.(m.P_g),
        pi_eps = Float64.(m.pi_eps),
        P_s = Float64.(m.P_s),
        K = Float64.(m.K),
        beta = Float64(m.beta),
        gamma = Float64(m.gamma),
        sigma_eps = Float64(m.sigma_eps),
        rstar = Float64(m.rstar),
        theta = Float64(m.theta),
        kappa = Float64(m.kappa),
        phi_g = Float64.(m.phi_g),
        nbar_g = Float64.(m.nbar_g),
        R_l = Float64(m.R_l),
        l_policy = Float64.(m.l_policy),
        pub = Float64(m.pub),
    )

    s = file["solution_data"]
    sol = nothing
    if s !== nothing
        d = BitArray(Float64.(s.vnd) .< Float64.(s.vd))
        sz4 = (model.Nb, model.Ng, model.Ns, model.Ne)
        Q = Float64.(s.Q)
        X = ones(Float64, sz4)
        X[d] .= Q[d]
        sol = Solution(
            Float64.(s.vnd),
            Float64.(s.vd),
            d,
            .!d,
            Q,
            X,
            Float64.(s.R),
            Float64.(s.n),
            BitArray(s.schedule_mask),
            Float64.(s.pdefault),
            Int.(s.b_policy_idx),
            Int(s.outer_iters),
            Int[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        )
    end

    return (
        phi_bad = file["phi_bad"],
        phi_good = file["phi_good"],
        beta = file["beta"],
        ps_bad = file["ps_bad"],
        status = file["status"],
        converged = file["converged"],
        guess_source = file["guess_source"],
        error_message = file["error_message"],
        model = model,
        sol = sol,
    )
end

function write_manifest(manifest, results_dir)
    open(joinpath(results_dir, "manifest.csv"), "w") do io
        println(io, "case_idx,phi_bad,phi_good,beta,ps_bad,status,converged,outer_iters,final_outer_err,guess_source,path")
        for row in manifest
            println(
                io,
                string(
                    row.case_idx, ",",
                    @sprintf("%0.4f", row.phi_bad), ",",
                    @sprintf("%0.4f", row.phi_good), ",",
                    @sprintf("%0.4f", row.beta), ",",
                    @sprintf("%0.2f", row.ps_bad), ",",
                    row.status, ",",
                    row.converged, ",",
                    row.outer_iters, ",",
                    isnan(row.final_outer_err) ? "" : @sprintf("%.8e", row.final_outer_err), ",",
                    row.guess_source, ",",
                    row.path,
                ),
            )
        end
    end
end

function read_saved_cases(results_dir)
    saved_cases = NamedTuple[]
    for path in sort(readdir(results_dir; join = true))
        endswith(path, ".jld2") || continue
        file = JLD2.load(path)
        haskey(file, "phi_bad") || continue
        push!(saved_cases, (
            params = (
                phi_bad = file["phi_bad"],
                phi_good = file["phi_good"],
                beta = file["beta"],
                ps_bad = file["ps_bad"],
            ),
            path = path,
            converged = file["converged"],
            has_solution = file["solution_data"] !== nothing,
        ))
    end
    return saved_cases
end

function run_calibration_grid(base_raw, phi_bad_vals, phi_good_vals, beta_vals, ps_bad_vals, results_dir)
    mkpath(results_dir)
    total_cases = length(phi_bad_vals) * length(phi_good_vals) * length(beta_vals) * length(ps_bad_vals)
    saved_cases = read_saved_cases(results_dir)
    manifest = NamedTuple[]
    last_sol = nothing

    for (case_idx, (ps_bad, beta, phi_bad, phi_good)) in enumerate([
        (ps_bad, beta, phi_bad, phi_good)
        for ps_bad in ps_bad_vals, beta in beta_vals, phi_bad in phi_bad_vals, phi_good in phi_good_vals
    ])
        p = (phi_bad = phi_bad, phi_good = phi_good, beta = beta, ps_bad = ps_bad)
        path = joinpath(results_dir, case_filename(phi_bad, phi_good, beta, ps_bad))

        set_params!(base_raw, phi_bad, phi_good, beta, ps_bad)
        model = init_model(base_raw)

        if isfile(path)
            file = JLD2.load(path)
            if haskey(file, "status") && haskey(file, "converged") && haskey(file, "solution_data")
                solution_data = file["solution_data"]
                push!(manifest, (
                    case_idx = case_idx,
                    phi_bad = phi_bad,
                    phi_good = phi_good,
                    beta = beta,
                    ps_bad = ps_bad,
                    status = file["status"],
                    converged = file["converged"],
                    outer_iters = solution_data === nothing ? 0 : Int(solution_data.outer_iters),
                    final_outer_err = solution_data === nothing ? NaN : Float64(solution_data.final_outer_err),
                    guess_source = "cached_file",
                    path = path,
                ))
                if solution_data !== nothing
                    last_sol = guess_solution(solution_data, model)
                end
                println("[$case_idx/$total_cases] cached: ", basename(path))
                continue
            end
            println("[$case_idx/$total_cases] rewriting old cache: ", basename(path))
        end

        guess = nothing
        guess_source = "none"
        pool = [row for row in saved_cases if row.has_solution]
        if !isempty(pool)
            converged_pool = [row for row in pool if row.converged]
            if !isempty(converged_pool)
                pool = converged_pool
            end
            row = pool[argmin([case_distance(p, row.params) for row in pool])]
            guess = guess_solution(JLD2.load(row.path)["solution_data"], model)
            guess_source = row.converged ? "nearest_converged" : "nearest"
        elseif last_sol !== nothing
            guess = last_sol
            guess_source = "last_solution"
        end

        println("[$case_idx/$total_cases] solving: ", basename(path))
        try
            sol = guess === nothing ? solve_model(model; verbose = false) : solve_model(model; sol = guess, verbose = false)
            converged = sol.outer_errs[end] <= model.tol
            status = converged ? "solved" : "not_converged"
            save_case(path, model, sol, phi_bad, phi_good, beta, ps_bad, status, converged, guess_source, "")

            push!(manifest, (
                case_idx = case_idx,
                phi_bad = phi_bad,
                phi_good = phi_good,
                beta = beta,
                ps_bad = ps_bad,
                status = status,
                converged = converged,
                outer_iters = sol.outer_iters,
                final_outer_err = sol.outer_errs[end],
                guess_source = guess_source,
                path = path,
            ))
            push!(saved_cases, (params = p, path = path, converged = converged, has_solution = true))
            last_sol = sol
        catch err
            error_message = sprint(showerror, err)
            save_case(path, model, nothing, phi_bad, phi_good, beta, ps_bad, "error", false, guess_source, error_message)

            push!(manifest, (
                case_idx = case_idx,
                phi_bad = phi_bad,
                phi_good = phi_good,
                beta = beta,
                ps_bad = ps_bad,
                status = "error",
                converged = false,
                outer_iters = 0,
                final_outer_err = NaN,
                guess_source = guess_source,
                path = path,
            ))
            push!(saved_cases, (params = p, path = path, converged = false, has_solution = false))
            println("[$case_idx/$total_cases] error: ", basename(path))
        end
    end

    write_manifest(manifest, results_dir)
    return (
        manifest = manifest,
        results_dir = results_dir,
        solved_cases = count(row -> row.status == "solved", manifest),
        not_converged_cases = count(row -> row.status == "not_converged", manifest),
        error_cases = count(row -> row.status == "error", manifest),
    )
end
