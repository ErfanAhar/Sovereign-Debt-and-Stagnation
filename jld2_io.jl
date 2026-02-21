using Dates
using JLD2

@inline _has_key(d, key::String) = haskey(d, key) || haskey(d, Symbol(key))
@inline _get_key(d, key::String) = haskey(d, key) ? d[key] : d[Symbol(key)]

"""
    save_jld2_file(path; model=nothing, sol=nothing, data=nothing, metadata=NamedTuple())

Standardized saver for JLD2 files used across notebooks.
Stores `model`/`sol` when provided, optional `data`, and optional scalar metadata.
Always writes `saved_at`.
"""
function save_jld2_file(path::AbstractString; model = nothing, sol = nothing, data = nothing,
    metadata::NamedTuple = NamedTuple())
    if model === nothing && sol === nothing && data === nothing
        error("save_jld2_file requires at least one of `model`, `sol`, or `data`.")
    end

    mkpath(dirname(path))
    JLD2.jldopen(path, "w") do f
        f["saved_at"] = string(now())
        if model !== nothing
            f["model"] = model
        end
        if sol !== nothing
            f["sol"] = sol
        end
        if data !== nothing
            f["data"] = data
        end
        for (k, v) in pairs(metadata)
            f[string(k)] = v
        end
    end
    return path
end

"""
    load_jld2_file(path) -> (raw, model, sol, data)

Standardized loader for JLD2 files used across notebooks.
Supports direct `model`/`sol` keys and legacy payloads nested under `data`.
`model` and/or `sol` are `nothing` if not found.
"""
function load_jld2_file(path::AbstractString)
    raw = JLD2.load(path)

    model = _has_key(raw, "model") ? _get_key(raw, "model") : nothing
    sol = _has_key(raw, "sol") ? _get_key(raw, "sol") : nothing
    data = _has_key(raw, "data") ? _get_key(raw, "data") : nothing

    if (model === nothing || sol === nothing) && data !== nothing
        if data isa NamedTuple
            model = model === nothing && hasproperty(data, :model) ? data.model : model
            sol = sol === nothing && hasproperty(data, :sol) ? data.sol : sol
        elseif data isa AbstractDict
            model = (model === nothing && _has_key(data, "model")) ? _get_key(data, "model") : model
            sol = (sol === nothing && _has_key(data, "sol")) ? _get_key(data, "sol") : sol
        end
    end

    return (raw = raw, model = model, sol = sol, data = data)
end
