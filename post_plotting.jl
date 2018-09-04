using JLD2
using Plots

R = [2,3,4,5,6]
r = 4
k = []

ϵν5_1_4 = 2.7888819610767125e-5
η = 0.00990409742169068

const N =64

version = ""
for r in R
    str = "$(N)-$(r)-t0352"
    @load pwd()*"/data/" * str * ".jld" k C lC
    #@load pwd()*"/data/$(N)-$(r).jld" k C lC
    #@load pwd()*"/data/$(N)-$(r)-t350.jld" k C lC
#
    p = plot()
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
                xscale = :log10,
                line = :dashdot)

    title!("R = $(r)")

        #savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel$(vel_component).png")
        #savefig(pwd()*"/figs/64-harmonics-r-$(r)-vel$(vel_component).pdf")
    savefig(pwd()*"/figs/" * str * version * "-vel.png")

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
        savefig(pwd()*"/figs/" * str *version * "-sum.png")

end
