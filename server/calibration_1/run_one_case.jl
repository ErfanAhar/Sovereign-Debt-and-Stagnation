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
total_cases = length(cases)

case_idx = if !isempty(ARGS)
    parse(Int, ARGS[1])
elseif haskey(ENV, "SLURM_ARRAY_TASK_ID")
    parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
else
    error("Pass a case index as ARGS[1] or run inside a Slurm array job.")
end

1 <= case_idx <= total_cases || error("case_idx=$(case_idx) is out of range 1:$(total_cases)")

results_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(@__DIR__, "result", "calibration_grid")
mkpath(results_dir)

p = cases[case_idx]
path = joinpath(results_dir, case_filename(p.phi_bad, p.phi_good, p.beta, p.ps_bad))

if isfile(path)
    file = JLD2.load(path)
    if haskey(file, "status") && haskey(file, "converged") && haskey(file, "solution_data")
        println("[$case_idx/$total_cases] cached: ", basename(path))
        exit()
    end
end

base_raw = Model(
    Nb = 3000,
    Ne = 17,
    max_iter = 100,
    max_iter_vd = 200,
    max_iter_x = 2000,
    pub = 0.65,
)

set_params!(base_raw, p.phi_bad, p.phi_good, p.beta, p.ps_bad)
model = init_model(base_raw)

guess = nothing
guess_source = "none"
saved_cases = [row for row in read_saved_cases(results_dir) if row.has_solution && row.path != path]
if !isempty(saved_cases)
    converged_pool = [row for row in saved_cases if row.converged]
    pool = isempty(converged_pool) ? saved_cases : converged_pool
    row = pool[argmin([case_distance(p, row.params) for row in pool])]
    guess = guess_solution(JLD2.load(row.path)["solution_data"], model)
    guess_source = row.converged ? "nearest_converged" : "nearest"
end

println("[$case_idx/$total_cases] solving: ", basename(path))

try
    sol = guess === nothing ? solve_model(model; verbose = false) : solve_model(model; sol = guess, verbose = false)
    converged = sol.outer_errs[end] <= model.tol
    status = converged ? "solved" : "not_converged"
    save_case(path, model, sol, p.phi_bad, p.phi_good, p.beta, p.ps_bad, status, converged, guess_source, "")
    println(
        "[$case_idx/$total_cases] done status=", status,
        " converged=", converged,
        " guess_source=", guess_source,
    )
catch err
    error_message = sprint(showerror, err)
    save_case(path, model, nothing, p.phi_bad, p.phi_good, p.beta, p.ps_bad, "error", false, guess_source, error_message)
    rethrow(err)
end
