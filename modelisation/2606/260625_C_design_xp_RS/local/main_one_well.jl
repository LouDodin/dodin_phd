using Plots
using Measures
using LaTeXStrings

# Paramètres
const r = 0.574619342477644
const K = 6.675449070379925e7
const α = 2.185e-5
const γ = 9.47815e-3

function simulate_one_well(t_end, prop_R0=1.0-1e-15, H0=1e5; dt=0.1)
    S = (1 - prop_R0) * H0
    R = prop_R0 * H0
    H = S+R

    t = collect(0:dt:t_end)
    S_vec = zeros(length(t))
    R_vec = zeros(length(t))

    S_vec[1] = S
    R_vec[1] = R

    for i in 2:length(t)
        H = S + R

        S_new = max(0.0, S + dt * (r*S*(1 - H/K) - α*S + γ*R))
        R_new = max(0.0, R + dt * (r*R*(1 - H/K) + α*S - γ*R))

        S, R = S_new, R_new

        S_vec[i] = S
        R_vec[i] = R
    end

    return t, S_vec, R_vec
end

t, S, R = simulate_one_well(200)

ytick_vals1   = [10.0^i for i in 0:9]
ytick_labels1 = [L"10^{%$i}" for i in 0:9]

pl = plot(
    layout = (2, 1),
    size = (3000, 2000),
    left_margin = 15mm,
    right_margin = 10mm,
    top_margin = 10mm,
    bottom_margin = 10mm,
    grid = true,
    ytickfontsize = 25,
    legendfontsize = 23,
    guidefontsize = 23,
    xtickfontsize = 23,
    titlefontsize = 23,
    xlabel = "Time [days]",
    legend = :bottomright,
    plot_title = "One well",
    plot_titlefontsize = 28
)

plot!(pl, t, S, subplot=1, label=" S", yscale=:log10, lw=8, ylims=(1e0, 1e9), yticks=(ytick_vals1, ytick_labels1), ylabel="Abundances\n[part/mL]")
plot!(pl, t, R, subplot=1, label=" R", lw=8)

plot!(pl, t, S./(S.+R), subplot=2, label=" prop_S", lw=8, ylabel="Proportions\n[-]")
plot!(pl, t, R./(S.+R), subplot=2, label=" prop_R", lw=8)

display(pl)