using JLD2
using Plots

R = [3]

for r in R
    @load pwd()*"/data/64-$(r).jld" k C lC
    plot(k[2:end], lC[2:end],
            xlabel = "k",
            ylabel = "l Cl",
            yscale = :log10,
            xscale = :log10,
            lw = 2,
            label  = "spherical harmonics spectrum"
            )

    plot!(k[2:end], k[2:end].^(-5/3)./2,
            xlabel = "k",
            ylabel = "l Cl",
            label  = "-5/3 law",
            lw = 2,
            yscale = :log10,
            xscale = :log10, line = :dashdot)

     title!("R = $(r)")

     savefig(pwd()*"/figs/64-harmonics-r-$(r).png")
     savefig(pwd()*"/figs/64-harmonics-r-$(r).pdf")
end
