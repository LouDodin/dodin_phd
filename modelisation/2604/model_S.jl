function model_S(dY, Y, p, t)

    # parameters
    r, K, m = p
    
    # variables
    S = Y[1]

    # ode
    dY[1] = r*S*(1-S/K) - m*S

end
