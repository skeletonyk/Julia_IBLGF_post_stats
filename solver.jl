using Statistics

include("utils.jl")
include("solver_shell.jl")
include("solver_harmonics_coefficients.jl")
include("solver_spectrum.jl")

function stats_shell(r)
    # method 2 prebuild mesh
    println("working on spherical cap of radius = ", r)

    # find needed block-spherical-shell
    lb=find_left(blocks, r)
    ub=find_right(blocks, r)
    blocks_cap = view(blocks, lb:ub)

    # build τ-ϕ mesh, assuming basicaly equally sized in Ω
    println("building shperical shell field")
    shell = shell_field(r, blocks_cap)

    # calculating flm
    N_l = Int(round(2*pi/block_info.spacing[1] /2)) # number of P_l s
    f_lm = harmonics_coefficients( N_l, shell )

    # C_l spectrum
    k, C, lC = spectrum(N_l, f_lm)

    return shell.f, k, C, lC
end
