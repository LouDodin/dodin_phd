function model_SV_no_delta(dY, Y, p, t)

    # parameters
    ϕ, β = p
    
    # variables
    S, V = Y

    # ode
    dY[1] = r*S*(1-S/K) - ϕ*S*V - m*S
    dY[2] = β*ϕ*S*V

end
