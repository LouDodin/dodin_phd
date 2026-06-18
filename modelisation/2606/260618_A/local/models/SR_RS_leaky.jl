module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, INITIAL_CONDITION

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SR_RS_leaky"

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
    (name = :α, lower = 1e-10, upper = 1,  description = "mutation rate S→R"),
    (name = :γ, lower = 1e-10, upper = 1,  description = "reversion rate R→S"),
    (name = :φ, lower = 1e-14, upper = 1e-7,  description = "interaction rate"),
    (name = :λ, lower = 1e-14, upper = 1.0,  description = "leakage parameter"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    α, γ, φ, λ = p[1], p[2], p[3], p[4]
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - H/K) - φ*S*Vi - α*S + γ*R
    dY[2] = r*R*(1 - H/K) + α*S - γ*R - λ*φ*R*Vi
    dY[3] = β*φ*(S + λ*R)*Vi - δ*Vi
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

end