## ===== Packages =====
using CSV
using DataFrames
using DifferentialEquations
using DataInterpolations
using Plots
using Measures
using Sundials
using LaTeXStrings

n_interv = 3


## ===== Input =====
df_HA = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_VA = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_HB = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_VB = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_HC = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_VC = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_HD = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_VD = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_HE = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_VE = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)


## ===== Time alignment =====
t_HA = df_HA[:,1]./24; HA = df_HA[:,2]
t_VA = df_VA[:,1]./24; VA = df_VA[:,2]

t_HB = df_HB[:,1]./24; HB = df_HB[:,2]
t_HB .-= (t_HB[1] - t_HA[end])
t_VB = df_VB[:,1]./24; VB = df_VB[:,2]
t_VB .-= (t_VB[1] - t_HA[end])

t_HC = df_HC[:,1]./24; HC = df_HC[:,2]
t_HC .-= (t_HC[1] - t_HB[end])
t_VC = df_VC[:,1]./24; VC = df_VC[:,2]
t_VC .-= (t_VC[1] - t_HB[end])

t_HD = df_HD[:,1]./24; HD = df_HD[:,2]
t_HD .-= (t_HD[1] - t_HC[end])
t_VD = df_VD[:,1]./24; VD = df_VD[:,2]
t_VD .-= (t_VD[1] - t_HC[end])

t_HE = df_HE[:,1]./24; HE = df_HE[:,2]
t_HE .-= (t_HE[1] - t_HD[end])
t_VE = df_VE[:,1]./24; VE = df_VE[:,2]
t_VE .-= (t_VE[1] - t_HD[end])


## ===== Cycles =====
cycles = [
    (tH=t_HA, H=HA, tV=t_VA, V=VA, u0=[HA[1], VA[1]]),
    (tH=t_HB, H=HB, tV=t_VB, V=VB, u0=[HB[1], VB[1]]),
    (tH=t_HC, H=HC, tV=t_VC, V=VC, u0=[HC[1], VC[1]]),
    (tH=t_HD, H=HD, tV=t_VD, V=VD, u0=[HD[1], VD[1]]),
    (tH=t_HE, H=HE, tV=t_VE, V=VE, u0=[HE[1], VE[1]]),
]


## ===== GLOBAL VECTORS =====
t_H = vcat([c.tH for c in cycles]...)
H   = vcat([c.H for c in cycles]...)

t_V = vcat([c.tV for c in cycles]...)
V   = vcat([c.V for c in cycles]...)


t_start = minimum([c.tH[1] for c in cycles])
t_end   = maximum([c.tH[end] for c in cycles])


## ===== Constants =====
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02


## ===== Polynomial (φ) =====
function parse_polynome(filepath::String)
    lines = readlines(filepath)

    function extract_float(s)
        m = match(r"([+-])?\s*([0-9]+\.[0-9]+(?:[eE][+-]?[0-9]+)?)", s)
        sign = (m.captures[1] == "-") ? "-" : ""
        return parse(Float64, sign * m.captures[2])
    end

    intervals = []
    i = 1
    while i <= length(lines)
        m = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]", lines[i])
        if m !== nothing
            t0 = parse(Float64, m.captures[1])
            a  = extract_float(split(lines[i+1],"=")[2])
            b  = extract_float(split(lines[i+2],"*")[1])
            c  = extract_float(split(lines[i+3],"*")[1])
            d  = extract_float(split(lines[i+4],"*")[1])
            push!(intervals, (t0,a,b,c,d))
            i += 5
        else
            i += 1
        end
    end

    function φ(t)
        for i in 1:length(intervals)-1
            t0 = intervals[i][1]
            t1 = intervals[i+1][1]
            if t ≥ t0 && t < t1
                dt = t - t0
                a,b,c,d = intervals[i][2:end]
                return exp(a + b*dt + c*dt^2 + d*dt^3)
            end
        end
        a,b,c,d = intervals[end][2:end]
        dt = t - intervals[end][1]
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return φ, intervals
end


poly_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_interv)_polynome.txt")
φ_global, intervals = parse_polynome(poly_file)

# check
t_check = collect(range(t_H[1], t_H[end], length=1000))
plot(t_check, φ_global.(t_check), yscale=:log10)


## ===== MODEL (IMPORTANT FIX) =====
function model!(dY, Y, p, t)
    S, V = Y
    φ = p(t)

    dY[1] = r*S*(1 - S/K) - φ*S*V
    dY[2] = β*φ*S*V - δ*V
end


## ===== SOLVE PER CYCLE =====
sols = []

for cyc in cycles

    t0 = cyc.tH[1]
    t1 = max(cyc.tH[end], cyc.tV[end])
    u0 = cyc.u0

    t_save = collect(range(t_H[1], t_H[end], length=1000))

    prob = ODEProblem(model!, u0, (t0, t1), φ_global)
    sol  = solve(prob, Tsit5(), saveat=t_save)

    push!(sols, sol)
end

L = 11

# plot phi, S, V
data_color = RGB(31/255, 119/255, 180/255)
model_color = RGB(255/255, 127/255, 14/255)
phi_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGBA(0.5, 0.5, 0.5, 0.4)

cycle_starts = [cyc.tH[1] for cyc in cycles[2:end]]

ytick_vals1   = [10.0^i for i in 3:1:8]
ytick_labels1 = [L"10^{%$i}" for i in 3:1:8]
ytick_vals2   = [10.0^i for i in 6:1:10]
ytick_labels2 = [L"10^{%$i}" for i in 6:1:10]
ytick_vals3   = [10.0^i for i in -L:1:-7]
ytick_labels3 = [L"10^{%$i}" for i in -L:1:-7]

p1 = plot(yscale=:log10, xlabel="t", ylabel="S(t)", legend=:bottomright, ylims=(1e3, 1e8), yticks=(ytick_vals1, ytick_labels1))
p2 = plot(yscale=:log10, xlabel="t", ylabel="V(t)", legend=:bottomleft, ylims=(1e6, 1e10), yticks=(ytick_vals2, ytick_labels2))
p3 = plot(yscale=:log10, xlabel="t", ylabel="phi(t)", legend=:bottomleft, ylims=(1e-11, 1e-7), yticks=(ytick_vals3, ytick_labels3))
vline!(p1, cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")
vline!(p2, cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")
vline!(p3, cycle_starts, lw=1, linestyle=:dash, color=cycle_color, label="Dilution")

scatter!(p1, t_H, H, label=" H data", color=data_color, markershape=:circle, markersize=6, markerstrokewidth=0)
scatter!(p2, t_V, V, label=" V data", color=data_color, markershape=:square, markersize=6, markerstrokewidth=0)
t_knots = [it[1] for it in intervals]
scatter!(p3, t_knots, φ_global.(t_knots), label=" Knots", color=model_color, markersize=6, markerstrokewidth=0)

for (i, sol) in enumerate(sols)
    println([φ_global(t) for t in sol.t])
    plot!(p1, sol.t, sol[1,:], lw=4, label= i==1 ? " H model" : nothing, color=model_color)
    plot!(p2, sol.t, sol[2,:], lw=4, label= i==1 ? " V model" : nothing, color=model_color)
    plot!(p3, sol.t, [φ_global(t) for t in sol.t], lw=4, label= i==1 ? " phi model" : nothing, color=model_color)
end
plot(plot_title="$(n_interv) intervalle(s) par cycle", p1, p2, p3, layout=(3, 1), size=(1200,1000), margins=10mm)

# Vecteurs globaux concaténés
t_all = Float64[]
H_all = Float64[]
V_all = Float64[]
phi_all = Float64[]

for sol in sols
    # extraction
    t = sol.t
    H = sol[1, :]
    V = sol[2, :]
    phi = [φ_global(tt) for tt in t]

    # concaténation
    append!(t_all, t)
    append!(H_all, H)
    append!(V_all, V)
    append!(phi_all, phi)
end

# Export en fichier texte
open(joinpath(@__DIR__, "240426_output/vectors.txt"), "w") do io
    for i in eachindex(t_all)
        println(io, "$(t_all[i]) $(H_all[i]) $(V_all[i]) $(phi_all[i])")
    end
end



## ===== PLOTS =====
p4 = plot(yscale=:log10, xscale=:log10, xlabel="S(t)", ylabel="phi(t)")
p5 = plot(yscale=:log10, xscale=:log10, xlabel="V(t)", ylabel="phi(t)")

ytick_vals   = [10.0^i for i in -L:1:-7]
ytick_labels = [L"10^{%$i}" for i in -L:1:-7]
xtick_vals1   = [10.0^i for i in 3:1:8]
xtick_labels1 = [L"10^{%$i}" for i in 3:1:8]
xtick_vals2   = [10.0^i for i in 6:10]
xtick_labels2 = [L"10^{%$i}" for i in 6:10]

for (i, sol) in enumerate(sols)
    plot!(p4, sol[1,:], [φ_global(t) for t in sol.t], lw=2, yscale=:log10, label="Cycle $i", xlims=(1e3, 1e8), xticks=(xtick_vals1, xtick_labels1))
    plot!(p5, sol[2,:], [φ_global(t) for t in sol.t], lw=2, yscale=:log10, label="Cycle $i", xlims=(1e6, 1e10), xticks=(xtick_vals2, xtick_labels2))
    scatter!(p4, [sol[1,1]], [φ_global(sol.t[1])], color=:black, label="", legend=:bottomleft)
    scatter!(p5, [sol[2,1]], [φ_global(sol.t[1])], color=:black, label="", legend=:topright)
end

plot(p4, p5, plot_title="$(n_interv) intervalles par cycle", layout=(1,2), size=(1200,450), margins=10mm, ylims=(1e-14, 1e-7), yticks=(ytick_vals, ytick_labels))