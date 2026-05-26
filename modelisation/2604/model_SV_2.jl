function model(dY, Y, phi_func, t)

    # parameters
    ϕ = phi_func(t)
    
    # variables
    S, V = Y

    # ode
    dY[1] = r*S*(1 - S/K) - ϕ*S*V
    dY[2] = β*ϕ*S*V - δ*V

end