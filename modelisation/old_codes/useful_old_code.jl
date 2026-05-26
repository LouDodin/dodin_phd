## ===== Color code =====
col_H = :green
host_palette = cgrad([:darkgreen, :chartreuse])
col_S = host_palette[0.0]
col_I = host_palette[0.5]
col_R = host_palette[1.0]

col_V = :red
virus_palette = cgrad([:darkred, :orangered])
col_Vi   = virus_palette[0.0]
col_Vdp  = virus_palette[0.5]
col_Vdip = virus_palette[1.0]

col_Ev = :blue


## ===== Autorun index =====
output_dir = joinpath(@__DIR__, "output_coculture_cycles")
isdir(output_dir) || mkdir(output_dir)

existing_files = readdir(output_dir)
run_indices = Int[]

for f in existing_files
    m = match(r"run_(\d+)_", f)
    if m !== nothing
        push!(run_indices, parse(Int, m.captures[1]))
    end
end

run_id = isempty(run_indices) ? 1 : maximum(run_indices) + 1
println("Run index = ", run_id)


## ===== Plot Videodrop =====
data = DataFrame(XLSX.readtable(file_path,temp))

for col in names(data)[2:end]
    data[!,col] = parse_number.(data[!,col])
end

t_H_host = data[!,"host t cyto H (day)"]
t_H_coevo = data[!,"coevo t cyto H (day)"]
t_V = data[!,"coevo t cyto V (day)"]
t_partic = data[!,"t videodrop (day)"]

H_host_cols = ["host cyto H $r" for r in replicates]
H_coevo_cols = ["coevo cyto H $r" for r in replicates]
V_cols = ["coevo cyto V $r" for r in replicates]
partic_cols = ["videodrop $r" for r in replicates]

scatter!(plt[i],
    t_partic,
    mean_partic,
    color=color_partic,
    marker=:+,
    markersize=8,
    markerstrokewidth=2,
    label="Virus videodrop"
)


## ===== Plot Ribbon =====

function minmax_ci(vals)
    v = collect(skipmissing(vals))
    if isempty(v)
        return (NaN, NaN, NaN)
    end
    m = mean(v)
    low = m - minimum(v)
    high = maximum(v) - m
    return (m, low, high)
end


mean_H_host = Float64[]
ci_low_H_host = Float64[]
ci_high_H_host = Float64[]

mean_H_coevo = Float64[]
ci_low_H_coevo = Float64[]
ci_high_H_coevo = Float64[]

mean_V = Float64[]
ci_low_V = Float64[]
ci_high_V = Float64[]

mean_partic = Float64[]

for row in eachrow(data)

    m,low,high = minmax_ci(row[H_host_cols])
    push!(mean_H_host,m)
    push!(ci_low_H_host,low)
    push!(ci_high_H_host,high)

    m,low,high = minmax_ci(row[H_coevo_cols])
    push!(mean_H_coevo,m)
    push!(ci_low_H_coevo,low)
    push!(ci_high_H_coevo,high)

    m,low,high = minmax_ci(row[V_cols])
    push!(mean_V,m)
    push!(ci_low_V,low)
    push!(ci_high_V,high)

    partic_vals = collect(skipmissing(row[partic_cols]))
    if isempty(partic_vals)
        push!(mean_partic, NaN)
    else
        push!(mean_partic, mean(partic_vals))
    end
end

plot!(plt[i],
    t_H_host,
    mean_H_host,
    ribbon=(ci_low_H_host,ci_high_H_host),
    yscale=:log10,
    color=color_H,
    linestyle=:dash,
    lw=2,
    fillalpha=0.2,
    grid=true,
    label="Phytoplankton control cytometer",
    ylabel="Concentration (parts/mL)",
    xlabel="Time (day)",
    ylim=ylims_list[i],
    yticks=yticks_list[i],
    legendfontsize=16,
    legend=false
)

plot!(plt[i],
    t_H_coevo,
    mean_H_coevo,
    ribbon=(ci_low_H_coevo,ci_high_H_coevo),
    yscale=:log10,
    color=color_H,
    lw=2,
    fillalpha=0.2,
    grid=true,
    label="Phytoplankton cytometer",
    legendfontsize=16,
    ylabel="Concentration (parts/mL)",
    xlabel="Time (day)",
    legend=false,
    ylim=ylims_list[i],
    yticks=yticks_list[i]
)

plot!(plt[i],
    t_V,
    mean_V,
    ribbon=(ci_low_V,ci_high_V),
    yscale=:log10,
    color=color_V,
    lw=2,
    fillalpha=0.2,
    label="Virus cytometer",
    legendfontsize=16
)


## ===== Unused models =====
S_model = ModelSpec(
    "S",
    [],
    [:μ, :k],
    function (dY, Y, p, t)

        μ, k = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k)
        dI   = zero(I)
        dR   = zero(R)
        dVi  = zero(Vi)
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = dS / max(S, 1e-12)
        dY[2] = zero(Y[2])
        dY[3] = zero(Y[3])
        dY[4] = zero(Y[4])
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

I_model = ModelSpec(
    "I",
    [:η],
    [:μ, :k, :η],
    function (dY, Y, p, t)

        μ, k, η = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = zero(S)
        dI   = - η*I
        dR   = zero(R)
        dVi  = β*η*I
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = zero(Y[1])
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

R_model = ModelSpec(
    "R",
    [:μ_r, :k_r],
    [:μ, :k, :μ_r, :k_r],
    function (dY, Y, p, t)

        μ, k, μ_r, k_r = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = zero(S)
        dI   = zero(I)
        dR   = μ_r*R*(1-H/k_r)
        dVi  = zero(Vi)
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = zero(Ev)

        dY[1] = zero(Y[1])
        dY[2] = zero(Y[2])
        dY[3] = dR / max(R, 1e-12)
        dY[4] = zero(Y[4])
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = zero(Y[7])
    end
)

SIViVdip_model = ModelSpec(
    "SIViVdip",
    [:φi, :β, :δ, :η, :φdip, :εdip, :σdip],
    [:μ, :k, :φi, :β, :δ, :η, :φdip, :εdip, :σdip],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η, φdip, εdip, σdip = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        invHcap = 1.0 - H / k
        V = Vi + Vdp + Vdip

        dS   = μ*S*invHcap - φi*S*Vi - φdip*S*Vdip
        dI   = φi*S*Vi + φdip*S*Vdip - η*I
        dR   = zero(R)
        dVi  = (1-εdip)*β*η*I - φi*H*Vi - φdip*H*Vdip - σdip*Vi
        dVdp = zero(Vdp)
        dVdip   = - φdip*S*Vdip + εdip*β*η*I + σdip*Vi - δ*Vdip

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = dVdp / max(Vdp, 1e-12)
        dY[6] = dVdip / max(Vdip, 1e-12)
    end
)

SIViEv_model = ModelSpec(
    "SIViEv",
    [:φi, :β, :δ, :η, :ν],
    [:μ, :k, :φi, :β, :δ, :η, :ν],
    function (dY, Y, p, t)

        μ, k, φi, β, δ, η, ν = p

        S   = exp(Y[1])
        I   = exp(Y[2])
        R   = exp(Y[3])
        Vi  = exp(Y[4])
        Vdp = exp(Y[5])
        Vdip = exp(Y[6])
        Ev = exp(Y[7])

        H = S + I + R
        V = Vi + Vdp + Vdip

        dS   = μ*S*(1-H/k) - φi*S*Vi - ν*S
        dI   = φi*S*Vi - η*I
        dR   = zero(R)
        dVi  = β*η*I - φi*(H+Ev)*Vi - δ*Vi
        dVdp = zero(Vdp)
        dVdip = zero(Vdip)
        dEv = ν*S

        dY[1] = dS / max(S, 1e-12)
        dY[2] = dI / max(I, 1e-12)
        dY[3] = zero(Y[3])
        dY[4] = dVi / max(Vi, 1e-12)
        dY[5] = zero(Y[5])
        dY[6] = zero(Y[6])
        dY[7] = dEv / max(Ev, 1e-12)
    end
)