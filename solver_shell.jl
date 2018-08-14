include("solver_shell_intrp.jl")

struct shell_field_t
    N_τ :: Int64
    N_ϕ :: Int64
    mesh_τ  :: Array{Float64,1}
    mesh_ϕ  :: Array{Float64,1}
    f :: Array{Float64,2}
end

function shell_field(r, blocks_cap)
    N_planar_2 = 4 * π * r^2 / mean(block_info.spacing)^2
    N_τ = Int(floor(sqrt(N_planar_2/2))+1)
    N_ϕ = N_τ * 2

    mesh_τ = collect(LinRange(-1, 1, N_τ))
    mesh_ϕ = collect(LinRange(-π, π, N_ϕ))
    f = zeros(N_τ, N_ϕ)
    Cartesian2sphere_intrp(f, r, N_τ, N_ϕ, mesh_τ, mesh_ϕ, blocks_cap)

    return shell_field_t(N_τ, N_ϕ, mesh_τ, mesh_ϕ, f)
end
