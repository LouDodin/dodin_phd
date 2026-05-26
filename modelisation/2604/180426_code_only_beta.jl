## ===== Import packages =====
using Dates
using CSV
using DataFrames
using Optim
using DifferentialEquations
using LogExpFunctions
using Statistics
using Plots
using Measures
using BlackBoxOptim
using LaTeXStrings


## ===== INPUT DATA =====
df_H1 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_V1 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle1.csv"), DataFrame)
df_H2 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_V2 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle2.csv"), DataFrame)
df_H3 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_V3 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle3.csv"), DataFrame)
df_H4 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_V4 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle4.csv"), DataFrame)
df_H5 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)
df_V5 = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_ReplicateA_Cycle5.csv"), DataFrame)

# --- β1 ---
β1 = 38.09257390233523
t_H1 = df_H1[1:6,1]./24; H1 = df_H1[1:6,2]
t_V1 = df_V1[1:6,1]./24; V1 = df_V1[1:6,2]

# --- β2 ---
β2 = 10
t_H2 = df_H1[6:end,1]./24; H2 = df_H1[6:end,2]
t_V2 = df_V1[6:end,1]./24; V2 = df_V1[6:end,2]

t_H3 = df_H2[1:4,1]./24; H3 = df_H2[1:4,2]
t_V3 = df_V2[1:4,1]./24; V3 = df_V2[1:4,2]

ratio = (t_H3[1]-t_H2[end])
t_H3 = t_H3 .- ratio
t_V3 = t_V3 .- ratio

# --- β3 ---
β3 = 4.28254103565e-10
t_H4 = df_H2[4:end,1]./24; H4 = df_H2[4:end,2]
t_V4 = df_V2[4:end,1]./24; V4 = df_V2[4:end,2]

t_H4 = t_H4 .- ratio
t_V4 = t_V4 .- ratio

t_H5 = df_H3[1:4,1]./24; H5 = df_H3[1:4,2]
t_V5 = df_V3[1:4,1]./24; V5 = df_V3[1:4,2]

ratio = (t_H5[1]-t_H4[end])
t_H5 = t_H5 .- ratio
t_V5 = t_V5 .- ratio

# --- β4 ---
β4 = 1.1296329738793916e-9
t_H6 = df_H3[4:end,1]./24; H6 = df_H3[4:end,2]
t_V6 = df_V3[4:end,1]./24; V6 = df_V3[4:end,2]

t_H6 = t_H6 .- ratio
t_V6 = t_V6 .- ratio

t_H7 = df_H4[1:5,1]./24; H7 = df_H4[1:5,2]
t_V7 = df_V4[1:5,1]./24; V7 = df_V4[1:5,2]

ratio = (t_H7[1]-t_H6[end])
t_H7 = t_H7 .- ratio
t_V7 = t_V7 .- ratio

# --- β5 ---
β5 = 2.01289564822654e-10
t_H8 = df_H4[5:end,1]./24; H8 = df_H4[5:end,2]
t_V8 = df_V4[5:end,1]./24; V8 = df_V4[5:end,2]

t_H8 = t_H8 .- ratio
t_V8 = t_V8 .- ratio

t_H9 = df_H5[1:3,1]./24; H9 = df_H5[1:3,2]
t_V9 = df_V5[1:3,1]./24; V9 = df_V5[1:3,2]

ratio = (t_H9[1]-t_H8[end])
t_H9 = t_H9 .- ratio
t_V9 = t_V9 .- ratio

# --- β6 ---
β6 = 1.2514025718429622e-9
t_H10 = df_H5[3:end,1]./24; H10 = df_H5[3:end,2]
t_V10 = df_V5[3:end,1]./24; V10 = df_V5[3:end,2]

t_H10 = t_H10 .- ratio
t_V10 = t_V10 .- ratio



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

t_HA = df_HA[1:end,1]./24; HA = df_HA[1:end,2]
t_VA = df_VA[1:end,1]./24; VA = df_VA[1:end,2]

t_HB = df_HB[1:end,1]./24; HB = df_HB[1:end,2]
ratio = (t_HB[1]-t_HA[end])
t_HB= t_HB .- ratio
t_VB = df_VB[1:end,1]./24; VB = df_VB[1:end,2]
t_VB= t_VB .- ratio

t_HC = df_HC[1:end,1]./24; HC = df_HC[1:end,2]
ratio = (t_HC[1]-t_HB[end])
t_HC= t_HC .- ratio
t_VC = df_VC[1:end,1]./24; VC = df_VC[1:end,2]
t_VC= t_VC .- ratio

t_HD = df_HD[1:end,1]./24; HD = df_HD[1:end,2]
ratio = (t_HD[1]-t_HC[end])
t_HD= t_HD .- ratio
t_VD = df_VD[1:end,1]./24; VD = df_VD[1:end,2]
t_VD= t_VD .- ratio

t_HE = df_HE[1:end,1]./24; HE = df_HE[1:end,2]
ratio = (t_HE[1]-t_HD[end])
t_HE= t_HE .- ratio
t_VE = df_VE[1:end,1]./24; VE = df_VE[1:end,2]
t_VE= t_VE .- ratio

t_H_all = vcat(t_HA, t_HB, t_HC, t_HD, t_HE)
H_all   = vcat(HA, HB, HC, HD, HE)

t_V_all = vcat(t_VA, t_VB, t_VC, t_VD, t_VE)
V_all   = vcat(VA, VB, VC, VD, VE)


## ===== CONSTANTS =====
r = 0.5746194091297323
K = 6.675446257207877e7
m = 1.4266424138490254e-11




## ===== β1 =====
include("model_SV_no_delta.jl")

ϕ = 1.4914407318477944e-8

p = [ϕ, β1]
u0 = [H1[1],V1[1]]
tspan = (t_H1[1],t_H1[6])
prob = ODEProblem(model_SV_no_delta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol1 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )



## ===== β2 =====
include("model_SV_no_delta_only_beta.jl")
p = [β2]

u0 = [H2[1],V2[1]]
tspan = (t_H2[1],t_H2[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol2 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


u0 = [H3[1],V3[1]]
tspan = (t_H3[1],t_H3[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol3 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


## ===== β3 =====
p = [β3]

u0 = [H4[1],V4[1]]
tspan = (t_H4[1],t_H4[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol4 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


u0 = [H5[1],V5[1]]
tspan = (t_H5[1],t_H5[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol5 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


## ===== β4 =====
p = [β4]

u0 = [H6[1],V6[1]]
tspan = (t_H6[1],t_H6[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol6 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


u0 = [H7[1],V7[1]]
tspan = (t_H7[1],t_H7[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol7 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


## ===== β5 =====
p = [β5]

u0 = [H8[1],V8[1]]
tspan = (t_H8[1],t_H8[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol8 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


u0 = [H9[1],V9[1]]
tspan = (t_H9[1],t_H9[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol9 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )


## ===== β6 =====
p = [β6]

u0 = [H10[1],V10[1]]
tspan = (t_H10[1],t_H10[end])
prob = ODEProblem(model_SV_no_delta_only_beta, u0, tspan, p)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

sol10 = solve(
            prob,
            Tsit5(),
            u0 = u0,
            p = p,
            reltol=1e-6,
            abstol=1e-6,
            isoutofdomain=isoutofdomain
        )





# ===== PLOT =====
ytick_vals1   = [10.0^i for i in 0:2:9]
ytick_labels1 = [L"10^{%$i}" for i in 0:2:9]
ytick_vals2   = [10.0^i for i in 2:2:10]
ytick_labels2 = [L"10^{%$i}" for i in 2:2:10]
ytick_vals3   = [10, 40, 500]
ytick_labels3 = ["10", "40", "500"]

pl_sim = plot(
    layout=(3,1),   # <-- 3 lignes maintenant
    size=(1500,1400),

    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,

    grid=true,

    yscale=:log10,
    xlims=(t_H_all[1], t_H_all[end]+1),

    ytickfontsize = 26,
    legendfontsize=20,
    guidefontsize=24,
    xtickfontsize=24,
    titlefontsize=24,
)

# ===== COLORS =====
data_color = RGB(31/255, 119/255, 180/255)
model_color = RGB(255/255, 127/255, 14/255)

phi_color = RGB(255/255, 127/255, 14/255)
cycle_color = RGBA(0.5, 0.5, 0.5, 0.4)
phi_change_color = RGB(255/255, 127/255, 14/255)

# ===== DATA =====
plot!(pl_sim[1], t_H_all, H_all,
    label=" Host data", lw=4,
    markershape=:circle, markersize=8,
    markerstrokewidth=0,
    color=data_color, linealpha=0.7,
    ylabel="Host abundance\n(cell/mL)",
    legend=:bottomright,
    ylims=(1e0, 1e9),
    yticks=(ytick_vals1, ytick_labels1)
)

plot!(pl_sim[2], t_V_all, V_all,
    label=" Virus data", lw=4,
    markershape=:square, markersize=8,
    markerstrokewidth=0,
    color=data_color, linealpha=0.7,
    ylabel="Virus abundance\n(particles/mL)",
    legend=:bottomright,
    ylims=(1e2, 1e10),
    yticks=(ytick_vals2, ytick_labels2)
)

# ===== SIMULATIONS =====
sols = [sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9, sol10]

for (i, sol) in enumerate(sols[1:2])
    lbl1 = i == 1 ? " Host model" : nothing
    lbl2 = i == 1 ? " Virus model" : nothing

    plot!(pl_sim[1], sol.t, sol[1,:], lw=6, color=model_color, label=lbl1)
    plot!(pl_sim[2], sol.t, sol[2,:], lw=6, color=model_color, label=lbl2)
end

# ===== CYCLE LINES =====
cycle_changes = [t_HB[1], t_HC[1], t_HD[1], t_HE[1]]

for (k, t_change) in enumerate(cycle_changes)
    for i in 1:3
        vline!(pl_sim[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label=k == 1 ? " Dilution" : nothing
        )
    end
end

# ===== PHI CHANGE LINES =====
phi_changes = [t_H2[1]]

for (k, t_change) in enumerate(phi_changes[1:1])
    for i in 1:3
        vline!(pl_sim[i], [t_change],
            color=phi_change_color,
            linestyle=:dot,
            lw=2,
            label=k == 1 ? " β shift" : nothing
        )
    end
end

# ===== PHI SEGMENTS =====
phi_segments = [
    (t_start=t_H1[1],  t_end=t_H1[end],  phi=β1),

    (t_start=t_H2[1],  t_end=t_H2[end],  phi=β2),
    (t_start=t_H3[1],  t_end=t_H3[end],  phi=β2),

    (t_start=t_H4[1],  t_end=t_H4[end],  phi=β3),
    (t_start=t_H5[1],  t_end=t_H5[end],  phi=β3),

    (t_start=t_H6[1],  t_end=t_H6[end],  phi=β4),
    (t_start=t_H7[1],  t_end=t_H7[end],  phi=β4),

    (t_start=t_H8[1],  t_end=t_H8[end],  phi=β5),
    (t_start=t_H9[1],  t_end=t_H9[end],  phi=β5),

    (t_start=t_H10[1], t_end=t_H10[end], phi=β6),
]

# ===== PHI PLOT (BOTTOM PANEL) =====
for seg in phi_segments[1:2]
    plot!(pl_sim[3],
        [seg.t_start, seg.t_end],
        [seg.phi, seg.phi],
        lw=6,
        color=phi_color,
        label=nothing
    )
end



plot!(pl_sim[3],
    yscale=:log10,
    ylims=(10, 510),
    yticks=(ytick_vals3, ytick_labels3),
    ylabel="β\n(–)",
    xlabel="Time (days)",
    legend=:bottomright
)