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
    refining_parm = 2
    @time shell = shell_field(r, blocks_cap, refining_parm )

    # calculating flm
    N_k = Int(round(2*pi/block_info.spacing[1] /2))  # number of P_l s
    @time f_lm, k = harmonics_coefficients( N_k, shell, r)

    println(k)

    # C_l spectrum
    k, C, lC = spectrum(k, N_k, f_lm, r)

    return shell.u, k, C, lC
end
