using JLD2
#using GR
using Plots
using LaTeXStrings

R =  [2,3,4,5]
r = 3
k = []

ϵν5_1_4 = 2.7888819610767125e-5
η  = 0.00990409742169068
#η = 10^-1.6

const N = 128

version = "-over-r"
p = plot(
xaxis =  (L" \eta k \sim \eta l/r", :log10),
yaxis =  (L"E(k) \sim l C_l", :log10),
yscale = :log10,
grid = :on,
ylim = (.1, 10^4),
xlim = (0.013, 0.32),
size = (500,360),
dpi = 400)

str_t = "0000"
str0 = "$(N)-t" * str_t
#gr()

str_decl="stampede_"
for r in R
    str = str_decl* "$(N)-$(r)-t" * str_t
    @load pwd()*"/data/" * str * ".jld" k C lC
    #@load pwd()*"/data/$(N)-$(r).jld" k C lC
    #@load pwd()*"/data/$(N)-$(r)-t350.jld" k C lC
#
    k = k * η
    lC = lC /ϵν5_1_4

    tmp = sum(lC, dims=[2]) /4/π / r^3
    plot!(k[2:end-1], tmp[2:end-1],
                lw = 2,
                label  = "r=$(r)"
                )

end

    plot!(k[2:end], k[2:end].^(-5/3) * 3e0,
                label  = "-5/3 law",
                lw = 2,

                line = :dashdot,
                )
    #title!(str0)

    savefig(pwd()*"/scaled_APS/" * str0 *version * "-sum.png")
