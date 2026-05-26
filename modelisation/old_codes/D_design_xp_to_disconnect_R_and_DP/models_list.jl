struct ModelSpec
    name::String
    params::Vector{Symbol}
    dynamics!::Function
end

SIVi_model = ModelSpec(
    "SIVi",
    [],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η,
        εdp, σdp,
        μ_r, k_r, ν = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi
        dI   = φi*S*Vi - η*I
        dR   = zero(R)
        dVi  = β*η*I - φi*H*Vi - δ*Vi
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

SIViVdp_model = ModelSpec(
    "SIViVdp",
    [],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η,
        εdp, σdp,
        μ_r, k_r, ν = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi
        dI   = φi*S*Vi - η*I
        dR   = zero(R)
        dVi  = (1-εdp)*β*η*I - φi*H*Vi - σdp*Vi
        dVdp = εdp*β*η*I - δ*Vdp + σdp*Vi
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = dVdp / max(Vdp, 1e-12)
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

SIRVi_SR_model = ModelSpec(
    "SIRVi_SR",
    [:α],
    function (dY, Y, p, t)

        μ, k, φi, ω, δ, η,
        εdp, σdp,
        μ_r, k_r, ν, α = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi - α*S
        dI   = φi*S*Vi - η*I
        dR   = μ_r*R*(1-H/k_r) + α*S
        dVi  = β*η*I - φi*H*Vi - δ*Vi
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = dR / max(R, 1e-12)
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

SIRVi_IR_model = ModelSpec(
    "SIRVi_IR",
    [:α],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η,
        εdp, σdp,
        μ_r, k_r, ν, α = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi
        dI   = φi*S*Vi - η*I - α*I
        dR   = μ_r*R*(1-H/k_r) + α*I
        dVi  = β*η*I - φi*H*Vi - δ*Vi
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = dR / max(R, 1e-12)
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

SIViEv_model = ModelSpec(
    "SIViEv",
    [],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η,
        εdp, σdp,
        μ_r, k_r, ν = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi - ν*S
        dI   = φi*S*Vi - η*I
        dR   = zero(R)
        dVi  = β*η*I - φi*(H+Ev)*Vi - δ*Vi
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = ν*S

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = dEv / max(Ev, 1e-12)
    end
)

const MODELS = Dict(
    "SIVi"     => SIVi_model,
    "SIViVdp"  => SIViVdp_model,
    "SIRVi_SR"     => SIRVi_SR_model,
    "SIRVi_IR"     => SIRVi_IR_model,
    "SIViEv"     => SIViEv_model
)