using Colors
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

## ===== COLORS =====
color_A = RGB(0.6, 0.8, 1.0)                 # light blue
color_B = RGB(31/255, 119/255, 180/255)     # original blue (KEEP THIS ONE)
color_C = RGB(0.0, 0.3, 0.7)                # dark blue

data_color = RGB(31/255, 119/255, 180/255)            # cycle change lines

vd_color = RGB(1.0, 0.2, 0.2)

## ===== FUNCTION: LOAD + SHIFT CYCLES =====
function load_rep(replicate, cycles)

    tH_all = Float64[]
    tV_all = Float64[]
    H_all  = Float64[]
    V_all  = Float64[]

    cycle_changes = Float64[]   # <-- store boundaries here
    prev_end = 0.0

    for c in cycles
        df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_all/hostData_coevoCondition_Temperature15_Replicate$(replicate)_Cycle$c.csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_all/virusData_coevoCondition_Temperature15_Replicate$(replicate)_Cycle$c.csv"), DataFrame)

        tH = df_H[:,1] ./ 24
        tV = df_V[:,1] ./ 24

        H = df_H[:,2]
        V = df_V[:,2]

        # shift time so cycles connect
        #shift = tH[1] - prev_end
        #tH .-= shift

        #shift2 = tV[1] - prev_end
        #tV .-= shift2

        append!(tH_all, tH)
        append!(tV_all, tV)
        append!(H_all, H)
        append!(V_all, V)

        prev_end = tH[end]

        # store cycle boundary (end of this cycle)
        push!(cycle_changes, prev_end)
    end
    cycle_changes = cycle_changes[1:end-1]

    return tH_all, H_all, tV_all, V_all, cycle_changes
end






## ===== LOAD EACH REPLICATE =====
t_HA, HA, t_VA, VA, cA = load_rep("A", 1:1)
t_HB, HB, t_VB, VB, cB = load_rep("B", 1:1)
t_HC, HC, t_VC, VC, cC = load_rep("C", 1:1)

cycle_changes = unique(sort(vcat(cA, cB, cC)))




## Videodrop 
# 26°C
VD_t = [8.77, 16.77, 24.77, 34.77, 43.79]
VD_A = [1.92E+09, 6.98E+08, 1.40E+08, 1.80E+08, 6.66E+08]
VD_B = [1.98E+09, 6.54E+08, 1e15, 1.58E+08, 3.09E+08]
VD_C = [2.89E+09, 5.20E+08, 8.10E+07, 2.63E+08, 1.95E+08]

# 20°C
VD_t = [24.77, 34.77, 43.79, 56.79, 66.79]
VD_A = [4.11E+08, 2.15E+08, 3.88E+08, 3.39E+08, 8.71E+08]
VD_B = [4.57E+08, 2.37E+08, 2.82E+08, 3.84E+08, 7.00E+08]
VD_C = [3.61E+08, 1.19E+08, 2.60E+08, 1e15, 9.08E+08]

# 15°C
VD_t = [66.81]
VD_A = [4.79E+07]
VD_B = [5.64E+07]
VD_C = [7.37E+07]

## ===== PLOT =====
ytick_vals1   = [10.0^i for i in 0:2:9]
ytick_labels1 = [L"10^{%$i}" for i in 0:2:9]

ytick_vals2   = [10.0^i for i in 2:2:10]
ytick_labels2 = [L"10^{%$i}" for i in 2:2:10]


pl_sim = plot(
    layout=(2,1),
    size=(1800,1000),

    left_margin=15mm,
    right_margin=10mm,
    top_margin=10mm,
    bottom_margin=10mm,

    grid=true,
    yscale=:log10,

    xlims=(0, 1603.0/24 + 1),

    ytickfontsize=26,
    legendfontsize=20,
    guidefontsize=24,
    xtickfontsize=24,
    titlefontsize=24,
)


## ===== HOST =====
plot!(pl_sim[1], t_HA, HA,
    label="Replicate A", lw=4,
    color=color_A, markershape=:circle, markersize=5,
    alpha=0.7,
    ylabel="Host abundance\n(cell/mL)",
    legend=nothing,
    ylims=(1e0, 1e9),
    yticks=(ytick_vals1, ytick_labels1)
)

plot!(pl_sim[1], t_HB, HB,
    label="Replicate B",
    lw=4,
    color=color_B, markershape=:circle, markersize=6,
    alpha=0.7
)

plot!(pl_sim[1], t_HC, HC,
    label="Replicate C",
    lw=4,
    color=color_C, markershape=:circle, markersize=5,
    alpha=0.7
)


## ===== VIRUS =====
plot!(pl_sim[2], t_VA, VA,
    label="Replicate A",
    lw=3,
    color=color_A, markershape=:square, markersize=5,
    alpha=0.7,
    ylabel="Virus abundance\n(part/mL)",
    legend=nothing,
    ylims=(1e2, 7e10),
    yticks=(ytick_vals2, ytick_labels2)
)

plot!(pl_sim[2], t_VB, VB,
    label="Replicate B",
    lw=4,
    color=color_B, markershape=:square, markersize=6,
    alpha=0.8
)

plot!(pl_sim[2], t_VC, VC,
    label="Replicate C",
    lw=3,
    color=color_C, markershape=:square, markersize=5,
    alpha=0.7
)

## ===== VIDEODROP SCATTER =====
scatter!(pl_sim[2], VD_t, VD_A,
    label="Videodrop measures",
    markershape=:+,
    markersize=14,
    markerstrokewidth=3,
    color=vd_color, alpha=0.7
)

scatter!(pl_sim[2], VD_t, VD_B,
    label=nothing,
    markershape=:+,
    markersize=14,
    markerstrokewidth=3,
    color=vd_color, alpha=0.7
)

scatter!(pl_sim[2], VD_t, VD_C,
    label=nothing,
    markershape=:+,
    markersize=14,
    markerstrokewidth=3,
    color=vd_color, alpha=0.7
)

## ===== CYCLE CHANGE VLINES =====
for t_change in cycle_changes
    for i in 1:2
        vline!(pl_sim[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label = t_change == cycle_changes[1] ? "Dilution" : nothing
        )
    end
end

display(pl_sim)