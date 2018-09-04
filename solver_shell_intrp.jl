using MuladdMacro
include("ddf.jl")
let
xyz = [0.0, 0.0, 0.0]
global function Cartesian2sphere_intrp(f :: Array{Float64,2}, r, N_θ, N_ϕ, mesh_θ, mesh_ϕ, blocks_cap, vel_component)
    for i = 1 : N_θ
        for j = 1 : N_ϕ
            θ = mesh_θ[i]
            ϕ = mesh_ϕ[j]

            # cal x y z for interpolation
            xyz[3] = r * cos(θ)
            xyz[1] = r * sin(θ)cos(ϕ)
            xyz[2] = r * sin(θ)sin(ϕ)

            get = (b, i, j, k)-> get_vel(b, i, j, k, vel_component)
            f[i,j] = intrpl(blocks_cap, xyz, get)

        end
    end
    f .-= mean(f)
end
end

let
    tmp = [0.0, 0.0, 0.0]
    global function block_needed(coord, block, width)
        @. tmp = (block.origin - coord)/block_info.spacing
        #println(tmp)
        if (maximum(tmp) > width) || (minimum(tmp) + block_info.dim[1] < -width)
            return false
        else
            #println(coord, block.origin)
            return true
        end
    end


    vrtx = [0.0, 0.0, 0.0]
    #value ::Float64  =0.0
    #ijk2n = 0
    #n = 0

    global function get_vel(b :: block_t, i :: Int, j :: Int, k :: Int, ind :: Int)
        n = block_info.dim[1]
        @muladd ijk2n = k * n * n + (j * n + i + 1)
        value = b.vel[ijk2n , ind]

        return value
    end
end

let
    δx = [0.0, 0.0, 0.0]
    #sum :: Float64= 0.0
    #weight :: Float64 = 0.0

    global function intrpl(blocks, coord, get_)
        sum :: Float64 = 0.0
        weight :: Float64 = 0.0
        spacing = block_info.spacing[1]
        width = 1 + 1e-5

        for b in blocks
            if block_needed(coord, b, width)
                # loop through elements
                @inbounds for k = 0:block_info.dim[3]-1
                    δx[1] = k + (b.origin[3] - coord[3])./spacing
                    if abs(δx[1]) >= width
                        continue
                    end
                    @inbounds for j = 0:block_info.dim[3]-1
                        δx[2] = j + (b.origin[2] - coord[2])./spacing
                        if abs(δx[2]) >= width
                            continue
                        end
                        @inbounds for i = 0:block_info.dim[3]-1
                            δx[3] = i + (b.origin[3] - coord[3])./spacing
                            if abs(δx[3]) >= width
                                continue
                            end
                            value = get_(b,i,j,k)
                            #@. δx = (vrtx - coord)/block_info.spacing

                            weight = ddf(δx)
                            weight *= value # reuse
                            sum += weight
                        end
                    end
                end
            end
        end

        return sum
    end
end
