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


## ===== PARAMÈTRES CODE 2 (ϕ et β variables) =====

# --- ϕ1 ---
ϕ1 = 1.4914407318477944e-8
β1 = 38.09257390233523
t_H1 = df_H1[1:6,1]./24; H1 = df_H1[1:6,2]
t_V1 = df_V1[1:6,1]./24; V1 = df_V1[1:6,2]

# --- ϕ2 ---
ϕ2 = 7.275605942568924e-11
β2 = 155.56397628529712
t_H2 = df_H1[6:end,1]./24; H2 = df_H1[6:end,2]
t_V2 = df_V1[6:end,1]./24; V2 = df_V1[6:end,2]

t_H3 = df_H2[1:4,1]./24; H3 = df_H2[1:4,2]
t_V3 = df_V2[1:4,1]./24; V3 = df_V2[1:4,2]

ratio = (t_H3[1]-t_H2[end])
t_H3 = t_H3 .- ratio
t_V3 = t_V3 .- ratio

# --- ϕ3 ---
ϕ3 = 7.275605942568924e-11
β3 = 155.56397628529712
t_H4 = df_H2[4:end,1]./24; H4 = df_H2[4:end,2]
t_V4 = df_V2[4:end,1]./24; V4 = df_V2[4:end,2]

t_H4 = t_H4 .- ratio
t_V4 = t_V4 .- ratio

t_H5 = df_H3[1:5,1]./24; H5 = df_H3[1:5,2]
t_V5 = df_V3[1:5,1]./24; V5 = df_V3[1:5,2]

ratio = (t_H5[1]-t_H4[end])
t_H5 = t_H5 .- ratio
t_V5 = t_V5 .- ratio

# --- ϕ4 ---
ϕ4 = 2.289407615156652e-9
β4 = 19.68312324282512

t_H6 = df_H3[5:end,1]./24; H6 = df_H3[5:end,2]
t_V6 = df_V3[5:end,1]./24; V6 = df_V3[5:end,2]

t_H6 = t_H6 .- ratio
t_V6 = t_V6 .- ratio

t_H7 = df_H4[1:5,1]./24; H7 = df_H4[1:5,2]
t_V7 = df_V4[1:5,1]./24; V7 = df_V4[1:5,2]

ratio = (t_H7[1]-t_H6[end])
t_H7 = t_H7 .- ratio
t_V7 = t_V7 .- ratio

# --- ϕ5 ---
ϕ5 = 1.4739293408884745e-11
β5 = 499.9999999780063

t_H8 = df_H4[5:end,1]./24; H8 = df_H4[5:end,2]
t_V8 = df_V4[5:end,1]./24; V8 = df_V4[5:end,2]

t_H8 = t_H8 .- ratio
t_V8 = t_V8 .- ratio

t_H9 = df_H5[1:3,1]./24; H9 = df_H5[1:3,2]
t_V9 = df_V5[1:3,1]./24; V9 = df_V5[1:3,2]

ratio = (t_H9[1]-t_H8[end])
t_H9 = t_H9 .- ratio
t_V9 = t_V9 .- ratio

# --- ϕ6 ---
ϕ6 = 6.040101161232231e-10
β6 = 80.56266321712702
t_H10 = df_H5[3:end,1]./24; H10 = df_H5[3:end,2]
t_V10 = df_V5[3:end,1]./24; V10 = df_V5[3:end,2]

t_H10 = t_H10 .- ratio
t_V10 = t_V10 .- ratio


## ===== PARAMÈTRES CODE 1 (ϕ variable, β fixe) =====

ϕ1_m1 = 1.4914407318477944e-8
β = 38.09257390233523

ϕ2_m1 = 9.209511129767475e-11
ϕ3_m1 = 4.28254103565e-10
ϕ4_m1 = 1.1296329738793916e-9
ϕ5_m1 = 2.01289564822654e-10
ϕ6_m1 = 1.2514025718429622e-9

# Recalcul des time vectors du code 1 (même logique que code 1)
t_H1_m1 = df_H1[1:6,1]./24;   H1_m1 = df_H1[1:6,2]
t_V1_m1 = df_V1[1:6,1]./24;   V1_m1 = df_V1[1:6,2]

t_H2_m1 = df_H1[6:end,1]./24; H2_m1 = df_H1[6:end,2]
t_V2_m1 = df_V1[6:end,1]./24; V2_m1 = df_V1[6:end,2]

t_H3_m1 = df_H2[1:4,1]./24;   H3_m1 = df_H2[1:4,2]
t_V3_m1 = df_V2[1:4,1]./24;   V3_m1 = df_V2[1:4,2]
ratio_m1 = t_H3_m1[1] - t_H2_m1[end]
t_H3_m1 = t_H3_m1 .- ratio_m1
t_V3_m1 = t_V3_m1 .- ratio_m1

t_H4_m1 = df_H2[4:end,1]./24; H4_m1 = df_H2[4:end,2]
t_V4_m1 = df_V2[4:end,1]./24; V4_m1 = df_V2[4:end,2]
t_H4_m1 = t_H4_m1 .- ratio_m1
t_V4_m1 = t_V4_m1 .- ratio_m1

t_H5_m1 = df_H3[1:4,1]./24;   H5_m1 = df_H3[1:4,2]
t_V5_m1 = df_V3[1:4,1]./24;   V5_m1 = df_V3[1:4,2]
ratio_m1 = t_H5_m1[1] - t_H4_m1[end]
t_H5_m1 = t_H5_m1 .- ratio_m1
t_V5_m1 = t_V5_m1 .- ratio_m1

t_H6_m1 = df_H3[4:end,1]./24; H6_m1 = df_H3[4:end,2]
t_V6_m1 = df_V3[4:end,1]./24; V6_m1 = df_V3[4:end,2]
t_H6_m1 = t_H6_m1 .- ratio_m1
t_V6_m1 = t_V6_m1 .- ratio_m1

t_H7_m1 = df_H4[1:5,1]./24;   H7_m1 = df_H4[1:5,2]
t_V7_m1 = df_V4[1:5,1]./24;   V7_m1 = df_V4[1:5,2]
ratio_m1 = t_H7_m1[1] - t_H6_m1[end]
t_H7_m1 = t_H7_m1 .- ratio_m1
t_V7_m1 = t_V7_m1 .- ratio_m1

t_H8_m1 = df_H4[5:end,1]./24; H8_m1 = df_H4[5:end,2]
t_V8_m1 = df_V4[5:end,1]./24; V8_m1 = df_V4[5:end,2]
t_H8_m1 = t_H8_m1 .- ratio_m1
t_V8_m1 = t_V8_m1 .- ratio_m1

t_H9_m1 = df_H5[1:3,1]./24;   H9_m1 = df_H5[1:3,2]
t_V9_m1 = df_V5[1:3,1]./24;   V9_m1 = df_V5[1:3,2]
ratio_m1 = t_H9_m1[1] - t_H8_m1[end]
t_H9_m1 = t_H9_m1 .- ratio_m1
t_V9_m1 = t_V9_m1 .- ratio_m1

t_H10_m1 = df_H5[3:end,1]./24; H10_m1 = df_H5[3:end,2]
t_V10_m1 = df_V5[3:end,1]./24; V10_m1 = df_V5[3:end,2]
t_H10_m1 = t_H10_m1 .- ratio_m1
t_V10_m1 = t_V10_m1 .- ratio_m1


## ===== ALL DATA (pour xlims) =====
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


## ===== ODE SOLVING — CODE 2 (ϕ et β variables) =====
include("model_SV_no_delta.jl")
include("model_SV_no_delta_phi_beta.jl")

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

function solve_seg(model, u0, tspan, p)
    prob = ODEProblem(model, u0, tspan, p)
    return solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, isoutofdomain=isoutofdomain)
end

# Code 2 solutions
sol1  = solve_seg(model_SV_no_delta,          [H1[1],V1[1]],   (t_H1[1],  t_H1[6]),   [ϕ1, β1])
sol2  = solve_seg(model_SV_no_delta_phi_beta, [H2[1],V2[1]],   (t_H2[1],  t_H2[end]), [ϕ2, β2])
sol3  = solve_seg(model_SV_no_delta_phi_beta, [H3[1],V3[1]],   (t_H3[1],  t_H3[end]), [ϕ2, β2])
sol4  = solve_seg(model_SV_no_delta_phi_beta, [H4[1],V4[1]],   (t_H4[1],  t_H4[end]), [ϕ3, β3])
sol5  = solve_seg(model_SV_no_delta_phi_beta, [H5[1],V5[1]],   (t_H5[1],  t_H5[end]), [ϕ3, β3])
sol6  = solve_seg(model_SV_no_delta_phi_beta, [H6[1],V6[1]],   (t_H6[1],  t_H6[end]), [ϕ4, β4])
sol7  = solve_seg(model_SV_no_delta_phi_beta, [H7[1],V7[1]],   (t_H7[1],  t_H7[end]), [ϕ4, β4])
sol8  = solve_seg(model_SV_no_delta_phi_beta, [H8[1],V8[1]],   (t_H8[1],  t_H8[end]), [ϕ5, β5])
sol9  = solve_seg(model_SV_no_delta_phi_beta, [H9[1],V9[1]],   (t_H9[1],  t_H9[end]), [ϕ5, β5])
sol10 = solve_seg(model_SV_no_delta_phi_beta, [H10[1],V10[1]], (t_H10[1], t_H10[end]),[ϕ6, β6])


## ===== ODE SOLVING — CODE 1 (ϕ variable, β fixe) =====
include("model_SV_no_delta_only_phi.jl")

sol1_m1  = solve_seg(model_SV_no_delta,         [H1_m1[1],V1_m1[1]],   (t_H1_m1[1],  t_H1_m1[6]),    [ϕ1_m1, β_fixed])
sol2_m1  = solve_seg(model_SV_no_delta_only_phi, [H2_m1[1],V2_m1[1]],  (t_H2_m1[1],  t_H2_m1[end]),  [ϕ2_m1])
sol3_m1  = solve_seg(model_SV_no_delta_only_phi, [H3_m1[1],V3_m1[1]],  (t_H3_m1[1],  t_H3_m1[end]),  [ϕ2_m1])
sol4_m1  = solve_seg(model_SV_no_delta_only_phi, [H4_m1[1],V4_m1[1]],  (t_H4_m1[1],  t_H4_m1[end]),  [ϕ3_m1])
sol5_m1  = solve_seg(model_SV_no_delta_only_phi, [H5_m1[1],V5_m1[1]],  (t_H5_m1[1],  t_H5_m1[end]),  [ϕ3_m1])
sol6_m1  = solve_seg(model_SV_no_delta_only_phi, [H6_m1[1],V6_m1[1]],  (t_H6_m1[1],  t_H6_m1[end]),  [ϕ4_m1])
sol7_m1  = solve_seg(model_SV_no_delta_only_phi, [H7_m1[1],V7_m1[1]],  (t_H7_m1[1],  t_H7_m1[end]),  [ϕ4_m1])
sol8_m1  = solve_seg(model_SV_no_delta_only_phi, [H8_m1[1],V8_m1[1]],  (t_H8_m1[1],  t_H8_m1[end]),  [ϕ5_m1])
sol9_m1  = solve_seg(model_SV_no_delta_only_phi, [H9_m1[1],V9_m1[1]],  (t_H9_m1[1],  t_H9_m1[end]),  [ϕ5_m1])
sol10_m1 = solve_seg(model_SV_no_delta_only_phi, [H10_m1[1],V10_m1[1]],(t_H10_m1[1], t_H10_m1[end]), [ϕ6_m1])


# ===== PLOT =====
ytick_vals1   = [10.0^i for i in 0:2:9]
ytick_labels1 = [L"10^{%$i}" for i in 0:2:9]
ytick_vals2   = [10.0^i for i in 2:2:10]
ytick_labels2 = [L"10^{%$i}" for i in 2:2:10]
ytick_vals3   = [10.0^i for i in -11:-7]
ytick_labels3 = [L"10^{%$i}" for i in -11:-7]
ytick_vals4   = [10, 20, 40, 80, 150, 500]
ytick_labels4 = ["10", "20", "40", "80", "150", "500"]

pl_sim = plot(
    layout=(2, 2),
    size=(2800,1100),
    left_margin=20mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=20mm,
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
data_color      = RGB(31/255,  119/255, 180/255)
model_color     = RGB(255/255, 127/255, 14/255)   # code 2 — ϕ et β variables
model_color_m1  = RGB(44/255,  160/255, 44/255)   # code 1 — ϕ variable, β fixe
phi_color       = RGB(255/255, 127/255, 14/255)
phi_change_color= RGB(255/255, 127/255, 14/255)

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

plot!(pl_sim[3], t_V_all, V_all,
    label=" Virus data", lw=4,
    markershape=:square, markersize=8,
    markerstrokewidth=0,
    color=data_color, linealpha=0.7,
    ylabel="Virus abundance\n(particles/mL)",
    legend=:bottomright,
    ylims=(1e2, 1e10),
    yticks=(ytick_vals2, ytick_labels2),
    xlabel="Time (days)"
)

# ===== SIMULATIONS CODE 2 (ϕ et β variables) =====
sols2 = [sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9, sol10]
for (i, sol) in enumerate(sols2)
    lbl1 = i == 1 ? " Host model (ϕ & β var.)" : nothing
    lbl2 = i == 1 ? " Virus model (ϕ & β var.)" : nothing
    plot!(pl_sim[1], sol.t, sol[1,:], lw=6, color=model_color,    label=lbl1)
    plot!(pl_sim[3], sol.t, sol[2,:], lw=6, color=model_color,    label=lbl2)
end

# ===== SIMULATIONS CODE 1 (ϕ variable, β fixe) =====
sols1 = [sol1_m1, sol2_m1, sol3_m1, sol4_m1, sol5_m1, sol6_m1, sol7_m1, sol8_m1, sol9_m1, sol10_m1]
for (i, sol) in enumerate(sols1)
    lbl1 = i == 1 ? " Host model (ϕ var., β fixed)" : nothing
    lbl2 = i == 1 ? " Virus model (ϕ var., β fixed)" : nothing
    plot!(pl_sim[1], sol.t, sol[1,:], lw=4, color=model_color_m1, linestyle=:dash, label=lbl1)
    plot!(pl_sim[3], sol.t, sol[2,:], lw=4, color=model_color_m1, linestyle=:dash, label=lbl2)
end

# ===== CYCLE LINES =====
cycle_changes = [t_HB[1], t_HC[1], t_HD[1], t_HE[1]]
for (k, t_change) in enumerate(cycle_changes)
    for i in 1:4
        vline!(pl_sim[i], [t_change],
            color=data_color, linestyle=:dot, lw=2,
            label=k == 1 ? " Dilution" : nothing
        )
    end
end

# ===== PHI CHANGE LINES =====
phi_changes = [t_H2[1], t_H4[1], t_H6[1], t_H8[1], t_H10[1]]
for (k, t_change) in enumerate(phi_changes)
    for i in 1:4
        vline!(pl_sim[i], [t_change],
            color=phi_change_color, linestyle=:dot, lw=2,
            label=k == 1 ? " ϕ and β shift" : nothing
        )
    end
end

# ===== PHI SEGMENTS =====
phi_segments = [
    (t_start=t_H1[1],  t_end=t_H1[end],  phi=ϕ1, beta=β1),
    (t_start=t_H2[1],  t_end=t_H2[end],  phi=ϕ2, beta=β2),
    (t_start=t_H3[1],  t_end=t_H3[end],  phi=ϕ2, beta=β2),
    (t_start=t_H4[1],  t_end=t_H4[end],  phi=ϕ3, beta=β3),
    (t_start=t_H5[1],  t_end=t_H5[end],  phi=ϕ3, beta=β3),
    (t_start=t_H6[1],  t_end=t_H6[end],  phi=ϕ4, beta=β4),
    (t_start=t_H7[1],  t_end=t_H7[end],  phi=ϕ4, beta=β4),
    (t_start=t_H8[1],  t_end=t_H8[end],  phi=ϕ5, beta=β5),
    (t_start=t_H9[1],  t_end=t_H9[end],  phi=ϕ5, beta=β5),
    (t_start=t_H10[1], t_end=t_H10[end], phi=ϕ6, beta=β6),
]

# ===== PHI SEGMENTS — CODE 1 =====
phi_segments_m1 = [
    (t_start=t_H1_m1[1],  t_end=t_H1_m1[end],  phi=ϕ1_m1),
    (t_start=t_H2_m1[1],  t_end=t_H2_m1[end],  phi=ϕ2_m1),
    (t_start=t_H3_m1[1],  t_end=t_H3_m1[end],  phi=ϕ2_m1),
    (t_start=t_H4_m1[1],  t_end=t_H4_m1[end],  phi=ϕ3_m1),
    (t_start=t_H5_m1[1],  t_end=t_H5_m1[end],  phi=ϕ3_m1),
    (t_start=t_H6_m1[1],  t_end=t_H6_m1[end],  phi=ϕ4_m1),
    (t_start=t_H7_m1[1],  t_end=t_H7_m1[end],  phi=ϕ4_m1),
    (t_start=t_H8_m1[1],  t_end=t_H8_m1[end],  phi=ϕ5_m1),
    (t_start=t_H9_m1[1],  t_end=t_H9_m1[end],  phi=ϕ5_m1),
    (t_start=t_H10_m1[1], t_end=t_H10_m1[end], phi=ϕ6_m1),
]

# ===== PHI PANEL (pl_sim[2]) — les deux modèles =====
for seg in phi_segments
    plot!(pl_sim[2], [seg.t_start, seg.t_end], [seg.phi, seg.phi],
        lw=6, color=phi_color, label=nothing)
end
for seg in phi_segments_m1
    plot!(pl_sim[2], [seg.t_start, seg.t_end], [seg.phi, seg.phi],
        lw=4, color=model_color_m1, linestyle=:dash, label=nothing)
end

# transitions verticales ϕ
for k in 1:length(phi_segments)-1
    t_conn = phi_segments[k].t_end
    plot!(pl_sim[2], [t_conn, t_conn], [phi_segments[k].phi, phi_segments[k+1].phi],
        lw=2, linestyle=:dot, color=phi_change_color, label="")
end
for k in 1:length(phi_segments_m1)-1
    t_conn = phi_segments_m1[k].t_end
    plot!(pl_sim[2], [t_conn, t_conn], [phi_segments_m1[k].phi, phi_segments_m1[k+1].phi],
        lw=2, linestyle=:dot, color=model_color_m1, label="")
end

plot!(pl_sim[2],
    yscale=:log10,
    ylims=(1e-11, 1e-7),
    yticks=(ytick_vals3, ytick_labels3),
    ylabel="ϕ\n(mL/cell/day)",
    legend=:topright
)

# ===== BETA PANEL (pl_sim[4]) — code 2 uniquement =====
for seg in phi_segments
    plot!(pl_sim[4], [seg.t_start, seg.t_end], [seg.beta, seg.beta],
        lw=6, color=phi_color, label=nothing)
end
for k in 1:length(phi_segments)-1
    t_conn = phi_segments[k].t_end
    plot!(pl_sim[4], [t_conn, t_conn], [phi_segments[k].beta, phi_segments[k+1].beta],
        lw=2, linestyle=:dot, color=phi_change_color, label="")
end

plot!(pl_sim[4],
    yscale=:log10,
    ylims=(10, 510),
    yticks=(ytick_vals4, ytick_labels4),
    ylabel="β\n(-)",
    xlabel="Time (days)",
    legend=:bottomright
)