## Import packages 
using Dates
using CSV
using DataFrames
using Plots
using Measures

## ===== Infos =====
color_A = RGB(0.6, 0.8, 1.0)
color_B = RGB(31/255, 119/255, 180/255)
color_C = RGB(0.0, 0.3, 0.7) 
model_color = RGB(255/255, 127/255, 14/255)

replicate_colors = [color_A, color_B, color_C]
replicates = ["A", "B", "C"]
n_cycles = 5

## ===== Input =====
# Store time and abundance vectors per replicate
t_H_all = Dict{String, Vector{Float64}}()
H_all   = Dict{String, Vector{Float64}}()
t_V_all = Dict{String, Vector{Float64}}()
V_all   = Dict{String, Vector{Float64}}()

for rep in replicates
    t_H_rep = Float64[]
    H_rep   = Float64[]
    t_V_rep = Float64[]
    V_rep   = Float64[]

    t_H_prev_end = nothing
    t_V_prev_end = nothing

    for cycle in 1:n_cycles
        df_H = CSV.read(joinpath(@__DIR__, "input/xp_input_20/hostData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)
        df_V = CSV.read(joinpath(@__DIR__, "input/xp_input_20/virusData_coevoCondition_Temperature20_Replicate$(rep)_Cycle$(cycle).csv"), DataFrame)

        t_H = df_H[:, 1] ./ 24
        H   = df_H[:, 2]
        t_V = df_V[:, 1] ./ 24
        V   = df_V[:, 2]

        # Shift cycles to be contiguous (same logic as before)
        if t_H_prev_end !== nothing
            shift = t_H[1] - t_H_prev_end
            t_H = t_H .- shift
            t_V = t_V .- shift
        end

        t_H_prev_end = t_H[end]

        append!(t_H_rep, t_H)
        append!(H_rep,   H)
        append!(t_V_rep, t_V)
        append!(V_rep,   V)
    end

    t_H_all[rep] = t_H_rep
    H_all[rep]   = H_rep
    t_V_all[rep] = t_V_rep
    V_all[rep]   = V_rep
end

## ===== CHECK THE DATA =====
pl_data = plot(layout=(1,2), size=(900,300), margins=5mm, legend=:topright)

for (i, rep) in enumerate(replicates)
    scatter!(pl_data[1], t_H_all[rep], H_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (cell/ml)",
        yscale=:log10, ylims=(1e2, 1e8))

    scatter!(pl_data[2], t_V_all[rep], V_all[rep],
        color=replicate_colors[i], label="Rep $rep",
        xlabel="time (days)", ylabel="abundances (part/ml)",
        yscale=:log10, ylims=(1e3, 1e10))
end

display(pl_data)