module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, INITIAL_CONDITION, PHI_EQUIV

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "S_R"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    r = 0.574619342477644,
    K = 6.675449070379925e7,
    β = 144.0,
    δ = 0.02,
)

# ── Fitted parameters : names, bounds (natural space), log-fit flag ──────────
#
# FITTED_PARAMS  : Vector of NamedTuples
#   name   – symbol used in logs / output
#   lower  – lower bound  (natural space)
#   upper  – upper bound  (natural space)


const FITTED_PARAMS = [
    (name = :ϕ, lower = 1e-14, upper = 1e-7,  description = "interaction rate"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    ϕ = p[1]
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - H/K) - ϕ*S*Vi
    dY[2] = r*R*(1 - H/K)
    dY[3] = β*ϕ*S*Vi - δ*Vi
end

# ── Initial condition builder ─────────────────────────────────────────────────
# Called at the start of each cycle.
#   H0, V0    – observed mean initial abundances
#   prop_S    – fraction of hosts that are susceptible (carried over from
#               previous cycle end; = 1.0 at cycle 1)
# Returns u0 :: Vector{Float64}
function INITIAL_CONDITION(H0, V0, prop_S)
    return [prop_S * H0, (1 - prop_S) * H0, V0]
end

# ── Phi_equiv ─────────────────────────────────────────────────
function PHI_EQUIV(Y_t, p)
    ϕ = p[1]
    S  = Y_t[1]
    H  = Y_t[1] + Y_t[2]
    return ϕ * S / H
end

end