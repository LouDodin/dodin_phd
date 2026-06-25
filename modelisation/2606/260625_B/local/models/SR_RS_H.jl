module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, INITIAL_CONDITION, PHI_EQUIV

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SR_RS_H"

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
    (name = :eps_H,   lower = 1e-10, upper = 1, description = "total rate (α+β)"),
    (name = :K_H,   lower = 1e-10, upper = 1e10, description = "Kp"),
    (name = :ϕ,   lower = 1e-14, upper = 1e-7, description = "interaction rate"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p_vec, t)
    eps_H, K_H, ϕ = p_vec[1], p_vec[2], p_vec[3]
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    alpha = eps_H*H/(H+K_H)
    gamma = eps_H*K_H/(H+K_H)
    dY[1] = r*S*(1 - H/K) - ϕ*S*Vi - alpha*S + gamma*R
    dY[2] = r*R*(1 - H/K) + alpha*S - gamma*R
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
    ϕ = p[3]
    S  = Y_t[1]
    H  = Y_t[1] + Y_t[2]
    return ϕ * S / H
end

end