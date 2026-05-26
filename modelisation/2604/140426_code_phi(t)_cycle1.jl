using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Statistics
using Plots
using Measures

t_global_start = now()

## ===== INPUT DATA =====
t_H, H = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
t_V, V = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()

replicates = ("A", "B", "C")

for rep in replicates
    for data in ("host", "virus")
        df = CSV.read(joinpath(@__DIR__, "input/xp_input_20°/$(data)Data_coevoCondition_Temperature20_Replicate$(rep)_cycle1.csv"), DataFrame)
        t = collect(skipmissing(df[:, 1])) ./ 24
        x = collect(skipmissing(df[:, 2]))
        if data == "host"
            push!(t_H, t); push!(H, x)
        else
            push!(t_V, t); push!(V, x)
        end
    end
end


## ===== CONSTANTS =====
const μ = 0.5881765172
const k = 6.0e7


## ===== MODEL =====
SIVi_model = function (dY, Y, p, t)
    φi, β, δ, η = p
    dY[1]  = μ*Y[1]*(1-(Y[1]+Y[2])/k) - φi*Y[1]*Y[3]
    dY[2]  = φi*Y[1]*Y[3] - η*Y[2]
    dY[3] = β*η*Y[2] - φi*Y[1]*Y[3] - δ*Y[3]
end


## ===== INITIAL CONDITIONS =====
Y0 = [(H[1][1]+H[2][1]+H[3][1])/3, 0, (V[1][1]+V[2][1]+V[3][1])/3]
tspan = (t_H[1][1], t_H[1][end])


## ===== PARAMETERS =====
#          phi   beta  delta  eta
θ0 = log.([1e-7, 100.0, 1e-3, 5.0])
lb = log.([1e-12, 10.0, 1e-5, 1.0])
ub = log.([1e-6, 300.0, 0.1,  8.0])


## ===== FIT =====
function objective(θ)
    prob = ODEProblem(SIVi_model, Y0, tspan, exp.(θ))
    sol = solve(prob, Tsit5(), abstol=1E-8, reltol=1E-8)
    err = 0.0
    for i in eachindex(t_H) # Each replicate
        YH = sol(t_H[i])
        YV = sol(t_V[i])
        println(YH)
        println(YV)
        err += sum((log.(YH[1, :] .+ YH[2, :]) .- log.(H[i])).^2)
        err += sum((log.(YV[3, :]) .- log.(V[i])).^2)
    end
    return err
end

res = optimize(objective, lb, ub, θ0, Fminbox(BFGS()))
θopt = exp.(Optim.minimizer(res))

println("Optimal parameters:")
println("φ  = ", θopt[1])
println("β  = ", θopt[2])
println("δ  = ", θopt[3])
println("η  = ", θopt[4])


## ===== SIMULATION =====
prob = ODEProblem(SIVi_model, Y0, tspan, θopt)
sol  = solve(prob, Rodas5(), saveat=range(tspan[1], tspan[2], length=500))
S_model  = sol[1, :]
I_model  = sol[2, :]
Vi_model = sol[3, :]
H_model = S_model .+ I_model


## ===== PLOTS =====
host_palette  = cgrad([:darkgreen, :chartreuse])
virus_palette = cgrad([:darkred, :orangered])

col_S = host_palette[0.0]
col_I = host_palette[0.5]

p = plot(layout = (1, 2), size = (1200, 500), yscale = :log10, legend=:bottomleft, margins=10mm)

# SUBPLOT 1 : H + Vi + data
for i in eachindex(t_H)
    scatter!(p[1], t_H[i], H[i], color = :green, alpha = 0.5, label = i == 1 ? "Host data" : false)
end
for i in eachindex(t_V)
    scatter!(p[1], t_V[i], V[i], color = :red, alpha = 0.5, label = i == 1 ? "Virus data" : false)
end

plot!(p[1], sol.t, H_model, lw = 3, color = :darkgreen, label = "H model")
plot!(p[1], sol.t, Vi_model, lw = 3, color = :darkred, label = "V model")

title!(p[1], "H and V")
xlabel!(p[1], "Time (days)")
ylabel!(p[1], "Concentration (parts/mL)")


# SUBPLOT 2 : S + I + Vi
plot!(p[2], sol.t, S_model, lw = 3, color = col_S, label = "S")
plot!(p[2], sol.t, I_model, lw = 3, color = col_I, label = "I")
plot!(p[2], sol.t, Vi_model, lw = 3, color = :darkred, label = "Vi")

title!(p[2], "S, I and V")
xlabel!(p[2], "Time (days)")
ylabel!(p[2], "Concentration (parts/mL)")

display(p)

println("Done in $(round((now() - t_global_start).value / 1000, digits=1))s")