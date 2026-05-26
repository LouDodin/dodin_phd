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


## ===== FUNCTION: LOAD + SHIFT CYCLES =====
function load_rep(replicate, cycles)

    tH_all = Float64[]
    tV_all = Float64[]
    H_all  = Float64[]
    V_all  = Float64[]

    cycle_changes = Float64[]   # <-- store boundaries here
    prev_end = 0.0

    for c in cycles
        df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_all/hostData_hostCondition_Temperature26_Replicate$(replicate)_Cycle$c.csv"), DataFrame)
       
        tH = df_H[:,1] ./ 24
        
        H = df_H[:,2]
        
        # shift time so cycles connect
        #shift = tH[1] - prev_end
        #tH .-= shift

        append!(tH_all, tH)
        append!(H_all, H)

        prev_end = tH[end]

        # store cycle boundary (end of this cycle)
        push!(cycle_changes, prev_end)
    end
    cycle_changes = cycle_changes[1:end-1]

    return tH_all, H_all, cycle_changes
end




## ===== LOAD EACH REPLICATE =====
t_HA, HA, cA = load_rep("A", 1:5)
t_HB, HB, cB = load_rep("B", 1:5)
t_HC, HC, cC = load_rep("C", 1:5)

cycle_changes = unique(sort(vcat(cA, cB, cC)))


## ===== PLOT =====
ytick_vals1   = [10.0^i for i in 0:2:9]
ytick_labels1 = [L"10^{%$i}" for i in 0:2:9]


pl_sim = plot(
    layout=(2,1),
    size=(1500,1000),

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
    legend=:bottomright,
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



## ===== CYCLE CHANGE VLINES =====
for t_change in cycle_changes
    for i in 1:2
        vline!(pl_sim[i], [t_change],
            color=data_color,
            linestyle=:dot,
            lw=2,
            label = t_change == cycle_changes[1] ? "Cycle change" : nothing
        )
    end
end

display(pl_sim)