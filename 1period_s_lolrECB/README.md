# One-period Sunspot Model With ECB Intervention

This folder implements Section 2, "LOLR Model: ECB version," of
`docs/LOLR_oneperioddebt_vsunspot.pdf`.

The sovereign chooses only private one-period debt. After the primary market
clears, private lenders compare holding the bond with queueing at the ECB
facility:

```julia
H = expected_Q / (1 + rstar)
L = intervention_prob * purchase_price + (1 - intervention_prob) * H
R = 1 / max(H, L)
```

The benchmark policy is active only during repayment in the low-growth,
bad-sunspot state. Its purchase price is
`1 / (1 + rstar + ecb_spread)`.

Run from this folder in Julia:

```julia
include("src/simple_s_lolrECB.jl")
model = init_model(Model(Nb = 1000, Ne = 17, ecb_spread = 0.04, intervention_prob = 0.50))
sol = solve_model(model)
sim = simulate(model, sol; T = 200_000, seed = 1234)
stats = summarize_simulation(sim, model; burnin = 10_000)
```

`main_s_lolrECB.ipynb` provides the same workflow and lists the policy grid
suggested in the note.
