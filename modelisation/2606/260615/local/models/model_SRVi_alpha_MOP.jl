module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, PHI_equiv,
       INITIAL_CONDITION, LOWER_BOUNDS, UPPER_BOUNDS, STATE_LABELS, LOG_FIT

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SRVi"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    r = 0.574619342477644,
    K = 6.675449070379925e7,
    β = 144.0,
    δ = 0.02,
    a = 1e-10,
    b = 0.0,
    φ = 5.458443669387811e-9
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

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r, K, β, δ, a, b, φ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ, FIXED_PARAMS.a, FIXED_PARAMS.b, FIXED_PARAMS.φ
    dY[1] = r*S*(1 - H/K) - φ*S*Vi - (a*Vi/H +b)*S
    dY[2] = r*R*(1 - H/K) + (a*Vi/H +b)
    dY[3] = β*φ*S*Vi - δ*Vi
end

# ── φ_equiv ─────────────────────────────────────────────────────────────
# Receives the ODE solution at a single time point (a vector) and the fitted
# parameter vector p.  Returns the apparent infection rate comparable to φ'(t).
function PHI_equiv(Y_t, p)
    S  = Y_t[1]
    H  = Y_t[1] + Y_t[2]
    return FIXED_PARAMS.φ * S / H
end

# ── Initial condition builder ─────────────────────────────────────────────────
# Called at the start of each cycle.
#   H0, V0    – observed mean initial abundances
#   prop_S    – fraction of hosts that are susceptible (carried over from
#               previous cycle end; = 1.0 at cycle 1)
# Returns u0 :: Vector{Float64}
function INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)
    return [prop_S * H0, (1 - prop_S) * H0, V0]
end

# ── State labels (for plots / logs) ──────────────────────────────────────────
# Recognised states : :H_total, :S, :R  (host panel)
#                     :V_total, :Vi, :Vd (virus panel)
# indices : tuple of state-vector positions that contribute to this curve.
const STATE_LABELS = [
    (state = :H_total, indices = (1, 2), label = "Model H",  color = :black, lw = 4, ls = :solid),
    (state = :S,       indices = (1,),   label = "Model S",  color = :red,   lw = 2, ls = :dash),
    (state = :R,       indices = (2,),   label = "Model R",  color = :green, lw = 2, ls = :dash),
    (state = :V_total, indices = (3,),   label = "Model V",  color = :black, lw = 4, ls = :solid),
    (state = :Vi,      indices = (3,),   label = "Model Vi", color = :red, lw = 2, ls = :dash),
]

end # module