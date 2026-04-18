include(joinpath(@__DIR__, "src", "ArellanoModel_s.jl"))
include(joinpath(@__DIR__, "src", "calibration_grid_s.jl"))

function calibration_cases(phi_bad_vals, phi_good_vals, beta_vals, ps_bad_vals)
    return vec([
        (
            phi_bad = phi_bad,
            phi_good = phi_good,
            beta = beta,
            ps_bad = ps_bad,
        )
        for ps_bad in ps_bad_vals, beta in beta_vals, phi_bad in phi_bad_vals, phi_good in phi_good_vals
    ])
end

phi_bad_vals = collect(0.80:0.01:0.99)
phi_good_vals = collect(0.80:0.01:0.99)
beta_vals = collect(range(0.50, 0.95, length = 5))
ps_bad_vals = [0.01, 0.25]

cases = calibration_cases(phi_bad_vals, phi_good_vals, beta_vals, ps_bad_vals)
results_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(@__DIR__, "result", "calibration_grid")

manifest = NamedTuple[]
for (case_idx, params) in enumerate(cases)
    path = joinpath(
        results_dir,
        case_filename(params.phi_bad, params.phi_good, params.beta, params.ps_bad),
    )

    if !isfile(path)
        push!(manifest, (
            case_idx = case_idx,
            phi_bad = params.phi_bad,
            phi_good = params.phi_good,
            beta = params.beta,
            ps_bad = params.ps_bad,
            status = "missing",
            converged = false,
            outer_iters = 0,
            final_outer_err = NaN,
            guess_source = "none",
            path = path,
        ))
        continue
    end

    file = JLD2.load(path)
    solution_data = haskey(file, "solution_data") ? file["solution_data"] : nothing
    push!(manifest, (
        case_idx = case_idx,
        phi_bad = params.phi_bad,
        phi_good = params.phi_good,
        beta = params.beta,
        ps_bad = params.ps_bad,
        status = haskey(file, "status") ? file["status"] : "unknown",
        converged = haskey(file, "converged") ? file["converged"] : false,
        outer_iters = solution_data === nothing ? 0 : Int(solution_data.outer_iters),
        final_outer_err = solution_data === nothing ? NaN : Float64(solution_data.final_outer_err),
        guess_source = haskey(file, "guess_source") ? file["guess_source"] : "unknown",
        path = path,
    ))
end

write_manifest(manifest, results_dir)

println("wrote manifest to ", joinpath(results_dir, "manifest.csv"))
println("solved=", count(row -> row.status == "solved", manifest))
println("not_converged=", count(row -> row.status == "not_converged", manifest))
println("error=", count(row -> row.status == "error", manifest))
println("missing=", count(row -> row.status == "missing", manifest))
