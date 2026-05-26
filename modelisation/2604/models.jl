struct ModelSpec
    name::String
    fit_params::Vector{Symbol}     # parameters to fit
    full_params::Vector{Symbol}    # all parameters
    dynamics!::Function
end

SIVi_model = ModelSpec(
    "SIVi",
    [:φi, :β, :δ, :η],
    [:μ, :k, :φi, :β, :δ, :η],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        Vi  = exp(Y[4])

        H = S + I
        V = Vi

        dS   = μ*S*(1-H/k) - φi*S*Vi
        dI   = φi*S*Vi - η*I
        dVi  = β*η*I - φi*H*Vi - δ*Vi

        dY[1] = dS/S
        dY[2] = dI/I
        dY[4] = dVi/Vi
    end
)

SIVi_2_model = ModelSpec(
    "SIVi_2",
    [:φi],
    [:μ, :k, :φi, :β, :δ, :η],
    function (dY, Y, p, t)
        μ, k, φi, β, δ, η = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        Vi  = exp(Y[4])

        H = S + I
        V = Vi

        dS   = μ*S*(1-H/k) - φi*S*Vi
        dI   = φi*S*Vi - η*I
        dVi  = β*η*I - φi*H*Vi - δ*Vi

        dY[1] = dS / S
        dY[2] = dI / I
        dY[4] = dVi / Vi
    end
)

const MODELS = Dict(
    "SIVi" => SIVi_model,
    "SIVi_2" => SIVi_2_model
)