using LinearAlgebra
using QuantEcon

Base.@kwdef struct Model
    # Grid sizes
    Nb::Int = 1000            # debt grid size
    Ng::Int = 2               # growth-state grid size
    Ns::Int = 2               # sunspot-state grid size
    Ne::Int = 17              # shock grid size

    # Grids
    b::Vector{Float64} = collect(range(-0.05, 1.0, length = Nb))  # debt grid
    g::Vector{Float64} = [0.96, 1.04]                             # growth states
    eps::Vector{Float64} = collect(range(-0.1, 0.1, length = Ne)) # shock grid

    # Transitions
    P_g::Matrix{Float64} = [0.60 0.40; 0.25 0.75]    # growth transition
    pi_eps::Vector{Float64} = fill(1.0 / Ne, Ne)     # iid shock probabilities
    P_s::Matrix{Float64} = [0.25 0.75; 0.25 0.75]    # sunspot transition (iid default)
    K::Matrix{Float64} = zeros(0, 0)                # kron(P_s, P_g), filled in init_model

    # Preferences and technology
    beta::Float64 = 0.75                 # discount factor
    gamma::Float64 = 2.0                 # CRRA risk aversion
    sigma_eps::Float64 = 0.023           # shock scale in output
    rstar::Float64 = 0.01                # world interest rate
    theta::Float64 = 0.1                 # re-entry probability
    kappa::Float64 = 0.75                # debt haircut in default
    phi_g::Vector{Float64} = [0.92, 0.88] # output loss in default
    nbar_g::Vector{Float64} = [1.0, 1.0] # issuance cap

    # IMF policy
    R_l::Float64 = 1.02                            # IMF gross rate
    l_policy::Array{Float64,2} = fill(0.0, Ng, Ne)  # IMF loan rule

    # Schedule constraint
    pub::Float64 = 0.65                # default-prob cutoff

    # Solver options
    max_iter::Int = 200               # outer iterations
    max_iter_vd::Int = 200            # default value iterations
    max_iter_x::Int = 2000             # price iterations
    tol::Float64 = 1e-6               # outer tolerance
    tol_vd::Float64 = 1e-6            # default tolerance
    tol_x::Float64 = 1e-6             # price tolerance
end

function init_model(raw::Model)
    b = raw.b
    g = raw.g
    P_g = raw.P_g
    P_s = raw.P_s

    Nb = length(b)
    Ne = raw.Ne
    Ng = length(g)
    Ns = size(P_s, 1)

    @assert size(P_g) == (Ng, Ng)
    @assert size(P_s, 1) == size(P_s, 2)

    sigma_u = 1.0
    n_std = 2.5
    mc = QuantEcon.tauchen(Ne, 0.0, sigma_u, 0.0, n_std)
    eps = collect(mc.state_values)
    P_eps = mc.p
    P_eps_pow = P_eps^10000
    pi_eps = vec(P_eps_pow[1, :])
    pi_eps ./= sum(pi_eps)

    l_policy = raw.l_policy
    phi_g = raw.phi_g
    nbar_g = raw.nbar_g

    K = kron(P_s, P_g)

    return Model(; Nb = Nb, Ng = Ng, Ns = Ns, Ne = Ne, b = b, g = g, eps = eps, P_g = P_g,
        pi_eps = pi_eps, P_s = P_s, K = K, beta = raw.beta, gamma = raw.gamma,
        sigma_eps = raw.sigma_eps, rstar = raw.rstar, theta = raw.theta, kappa = raw.kappa,
        phi_g = phi_g, nbar_g = nbar_g, R_l = raw.R_l, l_policy = l_policy, pub = raw.pub,
        max_iter = raw.max_iter, max_iter_vd = raw.max_iter_vd, max_iter_x = raw.max_iter_x,
        tol = raw.tol, tol_vd = raw.tol_vd, tol_x = raw.tol_x)
end

struct Solution
    # Value functions
    vnd::Array{Float64,4}       # value in repayment
    vd::Array{Float64,4}        # value in default

    # Default and re-entry
    d::BitArray{4}              # default indicator
    e::BitArray{4}              # re-entry indicator

    # Prices and schedules
    Q::Array{Float64,4}         # bond price
    X::Array{Float64,4}         # continuation price
    R::Array{Float64,3}         # gross rate schedule
    n::Array{Float64,3}         # issuance schedule
    schedule_mask::BitArray{3} # selected schedule by sunspot

    # Default probabilities
    pdefault::Array{Float64,3}  # default probability

    # Policy indexes
    b_policy_idx::Array{Int,4}  # optimal b' index

    # Diagnostics
    outer_iters::Int            # number of outer iterations
    vd_iters::Vector{Int}       # default-value iterations per outer loop
    x_iters::Vector{Int}        # price iterations per outer loop
    outer_errs::Vector{Float64} # vnd error per outer loop
    vd_errs::Vector{Float64}    # vd error per outer loop
    x_errs::Vector{Float64}     # price error per outer loop
    outer_times::Vector{Float64} # seconds per outer loop
end
