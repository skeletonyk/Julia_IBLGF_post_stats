using JLD2
using Plots

R = [3]
r = 3
k = []
for r in R
    plot()
    @load pwd()*"/data/64-$(r).jld" k C lC
    for vel_component = 1:3

        plot!(k[2:end], lC[2:end, vel_component],
                xlabel = "k",
                ylabel = "l Cl",
                yscale = :log10,
                xscale = :log10,
                lw = 2,
                label  = "spectrum - u$(vel_component)"
                )


    end
    plot!(k[2:end], k[2:end].^(-5/3),
                xlabel = "k",
                ylabel = "l Cl",
                label  = "-5/3 law",
                lw = 2,
                yscale = :log10,
                xscale = :log10, line = :dashdot)

                title!("R = $(r)")

        #savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel$(vel_component).png")
        #savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel$(vel_component).pdf")
    savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel.pdf")

    tmp = sum(lC,dims=[2])
    plot(k[2:end], tmp[2:end],
                xlabel = "k",
                ylabel = "l Cl",
                yscale = :log10,
                xscale = :log10,
                lw = 2,
                label  = "spherical harmonics spectrum"
                )

            plot!(k[2:end], k[2:end].^(-5/3),
                xlabel = "k",
                ylabel = "l Cl",
                label  = "-5/3 law",
                lw = 2,
                yscale = :log10,
                xscale = :log10, line = :dashdot)

                title!("R = $(r)")

        #savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel$(vel_component).png")
        savefig(pwd()*"/figs/64-harmonics-r-$(r)-sum.pdf")

end
