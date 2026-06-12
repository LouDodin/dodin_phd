## ===== Model : SRVi =====
# dS/dt  = r·S·(1 - H/K) - φ·S·Vi - α·S
# dR/dt  = r·R·(1 - H/K) + α·S
# dVi/dt = β·φ·S·Vi - δ·Vi
#
# State : Y = [S, R, Vi]
# H = S + R
#
# φ_equiv(t) = φ · S(t) / H(t)
# → cohérence avec φ'(t) mesurée expérimentalement

module ModelDef

using Colors

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, PHI_equiv,
       INITIAL_CONDITION, LOWER_BOUNDS, UPPER_BOUNDS, STATE_LABELS, LOG_FIT

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SVi"

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
#
# LOG_FIT : if true, optimisation runs in log-space for all parameters
#           (recommended when bounds span several orders of magnitude)

const LOG_FIT = true

const FITTED_PARAMS = [
    (name = :φ, lower = 1e-14, upper = 1e-7,  description = "interaction rate"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    φ = p[1]
    S, Vi = Y[1], Y[2]
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - S/K) - φ*S*Vi
    dY[2] = β*φ*S*Vi - δ*Vi
end

# ── φ_equiv ─────────────────────────────────────────────────────────────
# Receives the ODE solution at a single time point (a vector) and the fitted
# parameter vector p.  Returns the apparent infection rate comparable to φ'(t).
function PHI_equiv(Y_t, p)
    φ = p[1]
    return φ
end

# ── Initial condition builder ─────────────────────────────────────────────────
# Called at the start of each cycle.
#   H0, V0    – observed mean initial abundances
#   prop_S    – fraction of hosts that are susceptible (carried over from
#               previous cycle end; = 1.0 at cycle 1)
# Returns u0 :: Vector{Float64}
function INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)
    return [H0, V0]
end

# ── State labels (for plots / logs) ──────────────────────────────────────────
# Recognised states : :H_total, :S, :R  (host panel)
#                     :V_total, :Vi, :Vd (virus panel)
# indices : tuple of state-vector positions that contribute to this curve.
const STATE_LABELS = [
    (state = :H_total, indices = (1,), label = "Model H",  color = RGB(255/255, 127/255, 14/255), lw = 4, ls = :solid),
    (state = :V_total, indices = (2,),   label = "Model V",  color = RGB(255/255, 127/255, 14/255), lw = 4, ls = :solid),
]

end # module