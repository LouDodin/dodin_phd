function model_SV(dY, Y, p, t)

    # parameters
    ϕ, β, δ = p
    
    # variables
    S, V = Y

    # ode
    dY[1] = r*S*(1-S/K) - ϕ*S*V - m*S
    dY[2] = β*ϕ*S*V - δ*V

end
