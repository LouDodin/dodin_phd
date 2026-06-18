module ModelDef

export MODEL_NAME, FIXED_PARAMS, FITTED_PARAMS, ODE_MODEL!, PHI_equiv,
       INITIAL_CONDITION, LOWER_BOUNDS, UPPER_BOUNDS, STATE_LABELS, LOG_FIT

# ── Identity ────────────────────────────────────────────────────────────────
const MODEL_NAME = "SRVi_fit_a_b_phi"

# ── Fixed parameters ────────────────────────────────────────────────────────
const FIXED_PARAMS = (
    r = 0.574619342477644,
    K = 6.675449070379925e7,
    β = 144.0,
    δ = 0.02,
)

# ── Fitted parameters ────────────────────────────────────────────────────────
const LOG_FIT = true

struct FittedParam
    name::String
    lower::Float64
    upper::Float64
    description::String
end

const FITTED_PARAMS = [
    FittedParam("a",   1e-14, 1e-6,  "Host-density-dependent resistance acquisition rate"),
    FittedParam("b",   1e-6,  1.0,   "Basal resistance acquisition rate (day⁻¹)"),
    FittedParam("phi", 1e-12, 1e-6,  "Adsorption rate φ (mL/part/day)"),
]

# ── ODE ─────────────────────────────────────────────────────────────────────
function ODE_MODEL!(dY, Y, p, t)
    S, R, Vi = Y[1], Y[2], Y[3]
    H = S + R
    r  = FIXED_PARAMS.r
    K  = FIXED_PARAMS.K
    β  = FIXED_PARAMS.β
    δ  = FIXED_PARAMS.δ
    a, b, φ = p[1], p[2], p[3]

    dY[1] = r*S*(1 - H/K) - φ*S*Vi - (a*Vi/H + b)*S
    dY[2] = r*R*(1 - H/K) + (a*Vi/H + b)*S
    dY[3] = β*φ*S*Vi - δ*Vi
end

# ── φ_equiv ──────────────────────────────────────────────────────────────────
function PHI_equiv(Y_t, p)
    S = Y_t[1]
    H = Y_t[1] + Y_t[2]
    φ = p[3]
    return φ * S / H
end

# ── Initial condition builder ─────────────────────────────────────────────────
function INITIAL_CONDITION(H0, V0, prop_S, prop_Vi)
    return [prop_S * H0, (1 - prop_S) * H0, V0]
end

# ── State labels ──────────────────────────────────────────────────────────────
const STATE_LABELS = [
    (state = :H_total, indices = (1, 2), label = "Model H",  color = :black, lw = 4, ls = :solid),
    (state = :S,       indices = (1,),   label = "Model S",  color = :red,   lw = 2, ls = :dash),
    (state = :R,       indices = (2,),   label = "Model R",  color = :green, lw = 2, ls = :dash),
    (state = :V_total, indices = (3,),   label = "Model V",  color = :black, lw = 4, ls = :solid),
    (state = :Vi,      indices = (3,),   label = "Model Vi", color = :red,   lw = 2, ls = :dash),
]

end