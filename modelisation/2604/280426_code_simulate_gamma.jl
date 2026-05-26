## ===== Packages =====
using CSV
using DataFrames
using DifferentialEquations
using DataInterpolations
using Plots
using Measures
using Sundials

const γ = 1.003676e-05

## ===== Intervalle de simulation =====
const t_start = 0.0    # <-- en jours
const t_end   = 30.0   # <-- en jours

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

t_HA = df_HA[1:end,1]./24; HA = df_HA[1:end,2]
t_VA = df_VA[1:end,1]./24; VA = df_VA[1:end,2]

t_HB = df_HB[1:end,1]./24; HB = df_HB[1:end,2]
ratio = (t_HB[1]-t_HA[end])
t_HB = t_HB .- ratio
t_VB = df_VB[1:end,1]./24; VB = df_VB[1:end,2]
t_VB = t_VB .- ratio

t_HC = df_HC[1:end,1]./24; HC = df_HC[1:end,2]
ratio = (t_HC[1]-t_HB[end])
t_HC = t_HC .- ratio
t_VC = df_VC[1:end,1]./24; VC = df_VC[1:end,2]
t_VC = t_VC .- ratio

t_HD = df_HD[1:end,1]./24; HD = df_HD[1:end,2]
ratio = (t_HD[1]-t_HC[end])
t_HD = t_HD .- ratio
t_VD = df_VD[1:end,1]./24; VD = df_VD[1:end,2]
t_VD = t_VD .- ratio

t_HE = df_HE[1:end,1]./24; HE = df_HE[1:end,2]
ratio = (t_HE[1]-t_HD[end])
t_HE = t_HE .- ratio
t_VE = df_VE[1:end,1]./24; VE = df_VE[1:end,2]
t_VE = t_VE .- ratio


# Vecteurs globaux
t_H = vcat(t_HA, t_HB, t_HC, t_HD, t_HE)
H   = vcat(HA, HB, HC, HD, HE)
t_V = vcat(t_VA, t_VB, t_VC, t_VD, t_VE)
V   = vcat(VA, VB, VC, VD, VE)

# Filtrage sur l'intervalle choisi
mask_H = (t_H .>= t_start) .& (t_H .<= t_end)
mask_V = (t_V .>= t_start) .& (t_V .<= t_end)
t_H_sim = t_H[mask_H];  H_sim = H[mask_H]
t_V_sim = t_V[mask_V];  V_sim = V[mask_V]

## ===== Constants =====
const r = 0.574619342477644
const K = 6.675449070379925e7
const β = 144.0
const δ = 0.02


## ===== Parse & interpolate ϕ(t) =====
function parse_polynome(filepath::String)
    lines = readlines(filepath)

    function extract_float(s::AbstractString)
        m = match(r"([+-])?\s*([0-9]+\.[0-9]+(?:[eE][+-]?[0-9]+)?)", s)
        m === nothing && error("Cannot extract float from: \"$s\"")
        sign = (m.captures[1] == "-") ? "-" : ""
        return parse(Float64, sign * m.captures[2])
    end

    intervals = Vector{NTuple{6,Float64}}()
    i = 1
    while i <= length(lines)
        m_iv = match(r"Interval\s+\[([0-9eE+\-.]+),\s*([0-9eE+\-.]+)\]\s+days:", strip(lines[i]))
        if m_iv !== nothing
            t0 = parse(Float64, m_iv.captures[1])
            t1 = parse(Float64, m_iv.captures[2])
            a  = extract_float(split(strip(lines[i+1]), "=")[2])
            b  = extract_float(split(strip(lines[i+2]), "*")[1])
            c  = extract_float(split(strip(lines[i+3]), "*")[1])
            d  = extract_float(split(strip(lines[i+4]), "*")[1])
            push!(intervals, (t0, t1, a, b, c, d))
            i += 5
            continue
        end
        i += 1
    end

    t_lo = intervals[1][1]
    t_hi = intervals[end][2]

    function phi_raw(t::Real)
        tc  = clamp(Float64(t), t_lo, t_hi)
        idx = length(intervals)
        for k in eachindex(intervals)
            if tc <= intervals[k][2]; idx = k; break; end
        end
        t0, _, a, b, c, d = intervals[idx]
        dt = tc - t0
        return exp(a + b*dt + c*dt^2 + d*dt^3)
    end

    return phi_raw
end

n_intervals  = 4
output_dir   = joinpath(@__DIR__, "240426_output")
poly_file    = joinpath(output_dir, "$(n_intervals)_polynome.txt")

phi_raw = parse_polynome(poly_file)

N_GRID   = 10_000
t_grid   = collect(range(t_start, t_end, length=N_GRID))
phi_grid = phi_raw.(t_grid)
ϕ_interp = LinearInterpolation(phi_grid, t_grid)

## ===== ODE model =====
function model!(dY, Y, p, t)
    γ_local = p[1]
    S, R, V = Y
    ϕt = ϕ_interp(t)
    N  = S + R
    dY[1] = r*S*(1 - N/K) - ϕt*N*V - γ_local*S
    dY[2] = γ_local*S + r*R*(1 - N/K)
    dY[3] = β*ϕt*N*V - δ*V
end

## ===== Solve =====
# Condition initiale = première valeur dans l'intervalle
Y0    = [Float64(H_sim[1]), 0.0, Float64(V_sim[1])]
tspan = (t_start, t_end)

isoutofdomain(u, p, t) = any(x -> x < 0 || !isfinite(x), u)

prob = ODEProblem(model!, Y0, tspan, [γ])
sol  = solve(prob, CVODE_BDF(linear_solver=:Dense);
             reltol        = 1e-8,
             abstol        = 1e-10,
             isoutofdomain = isoutofdomain)

println("Retcode: ", sol.retcode)

## ===== Plot =====
pl = plot(layout=(1,2), size=(1100,380), margins=6mm, dpi=150)

scatter!(pl[1], t_H_sim, H_sim;
    label="data", xlabel="Time (days)", ylabel="Abundance (cell/ml)",
    yscale=:log10, ms=4, title="Host (S+R)   γ=$(γ)")
plot!(pl[1], sol.t, sol[1,:] .+ sol[2,:];
    label="S+R model", lw=2.5, legend=:bottomright)

scatter!(pl[2], t_V_sim, V_sim;
    label="data", xlabel="Time (days)", ylabel="Abundance (virus/ml)",
    yscale=:log10, ms=4, title="Virus")
plot!(pl[2], sol.t, sol[3,:];
    label="model", lw=2.5, legend=:bottomright)

out_png = joinpath(output_dir, "model_SRV_sim.png")
savefig(pl, out_png)
println("Plot saved → $out_png")