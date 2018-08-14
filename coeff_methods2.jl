include("coeff_methods_utils.jl")
using Statistics

function Cartesian2sphere(f, r, N_τ, N_ϕ, mesh_τ, mesh_ϕ, blocks_cap)
    for i = 1 : N_τ
        for j = 1 : N_ϕ
            τ = mesh_τ[i]
            ϕ = mesh_ϕ[j]

            # cal x y z for interpolation
            θ = acos(τ)
            z = r * cos(θ)
            x = r * sin(θ)cos(ϕ)
            y = r * sin(θ)sin(ϕ)

            get = (b, i, j, k)-> get_vel(b, i, j, k, 1)
            f[i,j] = intrpl(blocks_cap, [x,y,z], get)

        end
    end
    f .-= mean(f)
end

function method2(r, Nl)
# method 2 prebuild mesh
        println("working on spherical cap of radius = ", r)
        # find needed block-spherical-shell
        lb=find_left(blocks, r)
        ub=find_right(blocks, r)
        blocks_cap = view(blocks, lb:ub)

        # build τ-ϕ mesh, assuming basicaly equally sized in Ω
        N_planar_2 = 4 * π * r^2 / mean(block_info.spacing)^2
        N_τ = Int(floor(sqrt(N_planar_2/2))+1)
        N_ϕ = N_τ * 2

        # build Ω mesh
        mesh_τ = collect(LinRange(-1, 1, N_τ))
        mesh_ϕ = collect(LinRange(-π, π, N_ϕ))

        println("build Ω mesh")
        f = zeros(N_τ, N_ϕ)
        Cartesian2sphere(f, r, N_τ, N_ϕ, mesh_τ, mesh_ϕ, blocks_cap)

        # calculating flm
        println("f_lm")
        f_lm = zeros(Nl,2*Nl)
        Y_lm = zeros(N_τ, N_ϕ)

        for l = 0:Nl-1
            for m = -l : l
                for i = 1 : N_τ
                    for j = 1 : N_ϕ
                        τ = mesh_τ[i]
                        ϕ = mesh_ϕ[j]
                        # cal x y z for interpolation
                        θ = acos(τ)
                        Y_lm[i, j] = Y(l,m,τ,ϕ)
                    end
                end
                Y_lm .*= f
                f_lm[l+1,m+l+1] = sum(Y_lm ) / N_τ / N_ϕ * 4 * π
            end
        end

    k, C, lC = flm_2_k_C(Nl, f_lm)

    return f, k, C, lC
end
