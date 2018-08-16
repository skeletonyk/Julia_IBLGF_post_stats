include("solver_shell_intrp.jl")
using Interpolations

struct shell_field_t
    N_τ :: Int64
    N_ϕ :: Int64
    mesh_τ  :: Array{Float64,1}
    mesh_ϕ  :: Array{Float64,1}
    u :: Array{Float64,3}
end

struct shell_field_scalar_t
    N_τ :: Int64
    N_ϕ :: Int64
    mesh_τ  #:: Array{Float64,1}
    mesh_ϕ  #:: Array{Float64,1}
    f
end

function refining(fine :: shell_field_scalar_t, coarse :: shell_field_scalar_t)

    # linear interpolation
    interp_ =LinearInterpolation( (coarse.mesh_τ, coarse.mesh_ϕ), coarse.f)

    for j = 1 : fine.N_ϕ
        for i = 1 : fine.N_τ
            fine.f[i,j] = interp_(fine.mesh_τ[i], fine.mesh_ϕ[j])
        end
    end

    nothing
end

function shell_field(r, blocks_cap, refining_parm)

    println("building shperical shell field")

    # interpolation mesh
    N_planar_2 = 4 * π * r^2 / mean(block_info.spacing)^2
    N_τ = Int(floor(sqrt(N_planar_2/2))+1)
    N_ϕ = N_τ * 2
    mesh_τ = -1: 2/(N_τ-1):1 #collect(LinRange(-1, 1, N_τ))
    mesh_ϕ = -π: 2*π/(N_ϕ-1):π #collect(LinRange(-π, π, N_ϕ))
    f = zeros(N_τ, N_ϕ)

    #spline and refined mesh of u to have better integral
    N_τ_refined = (N_τ-1) * refining_parm * 4 + 1
    N_ϕ_refined = (N_τ-1) * refining_parm  + 1
    mesh_τ_refined = (LinRange(-1, 1, N_τ_refined ))
    mesh_ϕ_refined = (LinRange(-π, π, N_ϕ_refined ))
    u = zeros(N_τ_refined, N_ϕ_refined, 3)

    for vel_component = 1:3
        println("-- interpolating velocity field u$(vel_component) -- ")
        @time let
        Cartesian2sphere_intrp(f, r, N_τ, N_ϕ, mesh_τ, mesh_ϕ, blocks_cap, vel_component)

        refining(
        shell_field_scalar_t(N_τ_refined, N_ϕ_refined, mesh_τ_refined, mesh_ϕ_refined,
                                view(u,:,:, vel_component)),
        shell_field_scalar_t(N_τ, N_ϕ, mesh_τ, mesh_ϕ, f)
        )
        end
    end

    println("building shperical shell field --- end")


    return shell_field_t(N_τ_refined, N_ϕ_refined, mesh_τ_refined, mesh_ϕ_refined, u)
end
