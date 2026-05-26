using Plots
using Measures

n_interv_list = [1, 2, 3]

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

    return φ, intervals
end

## ===== PLOT φ(t) =====
p = plot(
    xlabel="t",
    ylabel="φ(t)",
    yscale=:log10,
    title="Comparaison φ(t)",
    legend=:topright,
    size=(900,500),
    margins=10mm
)

colors = [:blue, :red, :green, :purple, :orange]

t_plot = range(0, 67, length=1000)  # adapte si besoin

for (k, n_interv) in enumerate(n_interv_list)

    poly_file = joinpath(@__DIR__, "240426_output/knots_dilutions_$(n_interv)_polynome.txt")
    φ, intervals = parse_polynome(poly_file)

    plot!(
        p,
        t_plot,
        φ.(t_plot),
        lw=3,
        color=colors[k],
        label="n = $n_interv"
    )

    # knots (optionnel mais utile)
    t_knots = [it[1] for it in intervals]
    scatter!(p, t_knots, φ.(t_knots), color=colors[k], markersize=4, label="")
end

savefig(p, joinpath(@__DIR__, "240426_output/phi_comparison.png"))
display(p)