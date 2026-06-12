## ===== Model : SViVd =====
# dS/dt  = r·S·(1 - S/K) - φ·S·Vi
# dVi/dt = ε·β·φ·S·Vi - δ·Vi
# dVd/dt = (1-ε)·β·φ·S·Vi - δ·Vd
#
# State : Y = [S, Vi, Vd]
#
# φ_equiv(t) = φ · Vi(t) / V(t)

module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, PHI_equiv,
       INITIAL_CONDITION, LOWER_BOUNDS, UPPER_BOUNDS, STATE_LABELS, LOG_FIT

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SViVd"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    r = 0.5592225270686286,
    K = 7.29695252684594e7,
    β = 144.0,
    δ = 0.02,
)

# ── Fitted parameters ────────────────────────────────────────────────────────
const LOG_FIT = true

const FITTED_PARAMS = [
    (name = :φ,       lower = 1e-14, upper = 1e-5,  description = "interaction rate"),
    (name = :epsilon, lower = 1e-6,  upper = 1.0,   description = "fraction of infective virions"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    φ, ε = p[1], p[2]
    S, Vi, Vd = Y[1], Y[2], Y[3]
    r, K, β, δ = FIXED_PARAMS.r, FIXED_PARAMS.K, FIXED_PARAMS.β, FIXED_PARAMS.δ
    dY[1] = r*S*(1 - S/K) - φ*S*Vi
    dY[2] = ε*β*φ*S*Vi - δ*Vi
    dY[3] = (1 - ε)*β*φ*S*Vi - δ*Vd
end

# ── φ_equiv ─────────────────────────────────────────────────────────────
function PHI_equiv(Y_t, p)
    φ = p[1]
    Vi = Y_t[2]
    V = Y_t[2] + Y_t[3]
    return φ * Vi / V
end

# ── Initial condition builder ─────────────────────────────────────────────────
# prop_S  is ignored (no resistant class), included for API compatibility.
# prop_Vi splits V0 between Vi and Vd (carried over from previous cycle end).
function INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)
    return [H0, prop_Vi * V0, (1 - prop_Vi) * V0]
end

# ── State labels ──────────────────────────────────────────────────────────────
# Recognised states : :H_total, :S, :R  (host panel)
#                     :V_total, :Vi, :Vd (virus panel)
# indices : tuple of state-vector positions that contribute to this curve.
const STATE_LABELS = [
    (state = :H_total, indices = (1,),    label = "Model H",  color = :black,  lw = 4, ls = :solid),
    (state = :S,       indices = (1,),    label = "Model S",  color = :red,  lw = 2, ls = :dash),
    (state = :V_total, indices = (2, 3),  label = "Model V",  color = :black,  lw = 4, ls = :solid),
    (state = :Vi,      indices = (2,),    label = "Model Vi", color = :red,   lw = 2, ls = :dash),
    (state = :Vd,      indices = (3,),    label = "Model Vd", color = :green, lw = 2, ls = :dash),
]

end # module