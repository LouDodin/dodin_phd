function model(dY, Y, p, t)

    # parameters
    r, K = p
    
    # variables
    S = Y[1]

    # ode
    dY[1] = r*S*(1-S/K)

end
