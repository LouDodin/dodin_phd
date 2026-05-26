function SIV_model(dY, Y, p, t)

    # this should work but it is better to be explicit and really see the model in a more userfriendly way (see bellow)
    # also it is weird to have a decay rate for the virus and nothing for the host

    #φi, β, δ, η = p
    #dY[1]  = μ*Y[1]*(1-(Y[1]+Y[2])/k) - φi*Y[1]*Y[3]
    #dY[2]  = φi*Y[1]*Y[3] - η*Y[2]
    #dY[3] = β*η*Y[2] - φi*Y[1]*Y[3] - δ*Y[3]

    # parameters
    r, K, m, ϕ, β, η, δ = p
    
    # variables
    S, E1, E2, E3, E4, I, V = Y

    # ode
    dY[1] = r*S*(1-(S+E1+E2+E3+E4+I)/K) - ϕ*S*V -m*S
    dY[2] = ϕ*S*V - 5*η*E1 -m*E1
    dY[3] = 5*η*E1 - 5*η*E2 -m*E2
    dY[4] = 5*η*E2 - 5*η*E3 -m*E3
    dY[5] = 5*η*E3 - 5*η*E4 -m*E4
    dY[6] = 5*η*E4 - 5*η*I -m*I
    dY[7] = β*η*I - ϕ*S*V - δ*V

end
