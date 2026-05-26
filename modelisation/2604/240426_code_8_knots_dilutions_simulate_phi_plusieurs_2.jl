using Plots
using Measures
using DifferentialEquations
using CSV
using DataFrames
using DataInterpolations
using Sundials
using LaTeXStrings

n_interv_list = [1, 2, 3]


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


## ===== POLYNOME =====
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

    return φ
end

## ===== MODEL =====
function model!(dY, Y, p, t)
    S, V = Y
    φ = p(t)

    dY[1] = r*S*(1 - S/K) - φ*S*V
    dY[2] = β*φ*S*V - δ*V
end

## ===== CONSTANTS =====
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02

colors = [:blue, :red, :green, :purple, :orange]

## ===== PLOTS =====
xtick_vals   = [10.0^i for i in -17:1:-7]
xtick_labels = [L"10^{%$i}" for i in -17:1:-7]
ytick_vals   = [10.0^i for i in -2:1:7]
ytick_labels = [L"10^{%$i}" for i in -2:1:7]

p1 = plot(
    xlabel="t",
    ylabel="φ·S·V",
    yscale=:log10,
    title="φ·S·V(t)",
    legend=:topright,
    ylims=(1e-2, 1e7),
    yticks=(ytick_vals, ytick_labels)
)

p2 = plot(
    xlabel="φ(t)",
    ylabel="φ·S·V",
    xscale=:log10,
    yscale=:log10,
    title="φ·S·V vs φ",
    legend=:bottomright,
    xlims=(1e-17, 1e-7),
    xticks=(xtick_vals, xtick_labels),
    ylims=(1e-2, 1e7),
    yticks=(ytick_vals, ytick_labels)
)

## ===== LOOP =====
for (k, n_interv) in enumerate(n_interv_list)

    poly_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_interv)_polynome.txt")
    φ = parse_polynome(poly_file)

    all_t = Float64[]
    all_phi = Float64[]
    all_flux = Float64[]

    ## ===== SOLVE cycles =====
    for cyc in cycles

        t0 = cyc.tH[1]
        t1 = max(cyc.tH[end], cyc.tV[end])

        prob = ODEProblem(model!, cyc.u0, (t0, t1), φ)
        sol  = solve(prob, Tsit5(), saveat=range(t0, t1, length=400))

        for (i, t) in enumerate(sol.t)
            S = sol[1,i]
            V = sol[2,i]
            phi_t = φ(t)

            push!(all_t, t)
            push!(all_phi, phi_t)
            push!(all_flux, phi_t * S * V)
        end
    end

    ## ===== SORT (important pour plot propre) =====
    idx = sortperm(all_t)
    all_t = all_t[idx]
    all_phi = all_phi[idx]
    all_flux = all_flux[idx]

    ## ===== PLOT =====
    plot!(p1, all_t, all_flux, lw=3, color=colors[k], label="n = $n_interv")
    plot!(p2, all_phi, all_flux, lw=3, color=colors[k], label="n = $n_interv")
end

## ===== FINAL FIGURE =====
fig = plot(p1, p2, layout=(1,2), size=(1500,600), margins=10mm, xtickfontsize=14, ytickfontsize=16, guidefontsize=16, legendfontsize=16, legend=:bottomright)

savefig(fig, joinpath(@__DIR__, "240426_output/phiSV_comparison.png"))
display(fig)