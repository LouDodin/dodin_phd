module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SV_r"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    K = 6.675449070379925e7,
    δ = 0.02,
)

# ── Fitted parameters : names, bounds (natural space), log-fit flag ──────────
#
# FITTED_PARAMS  : Vector of NamedTuples
#   name   – symbol used in logs / output
#   lower  – lower bound  (natural space)
#   upper  – upper bound  (natural space)


const FITTED_PARAMS = [
    (name = :ϕ,   lower = 1e-14, upper = 1e-7, description = "interaction rate"),
    (name = :β,   lower = 1, upper = 500, description = "burst size"),
    (name = :r,   lower = 1e-5, upper = 1e1, description = "host growth rate"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p_vec, t)
    ϕ, β, r_func = p_vec[1], p_vec[2], p_vec[3]
    r = r_func(t)
    S, Vi = Y[1], Y[2]
    K, δ = FIXED_PARAMS.K, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - S/K) - ϕ*S*Vi
    dY[2] = β*ϕ*S*Vi - δ*Vi
end

end