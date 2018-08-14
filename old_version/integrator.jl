using MuladdMacro
function Y(l,m, τ, ϕ)
    if m>0
        return GSL.sf_legendre_sphPlm(l, m, τ) * cos(m * ϕ)*sqrt(2)
    elseif m==0
        return GSL.sf_legendre_sphPlm(l, m, τ) * cos(m * ϕ)
    else
        return GSL.sf_legendre_sphPlm(l, -m, τ) * sin(-m * ϕ)*sqrt(2)
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
        #@. vrtx = [i, j, k] * block_info.spacing + b.origin
        @muladd begin
        vrtx[1] = i * block_info.spacing[1] + b.origin[1]
        vrtx[2] = j * block_info.spacing[2] + b.origin[2]
        vrtx[3] = k * block_info.spacing[3] + b.origin[3]
        end

        # using 'F' type of storing order
        n = block_info.dim[1]

        @muladd ijk2n = k * n * n + (j * n + i + 1)
        value = b.vel[ijk2n , ind]

        return vrtx, value
        #nothing
    end
end

let
    δx = [0.0, 0.0, 0.0]
    #sum :: Float64= 0.0
    #weight :: Float64 = 0.0

    global function intrpl(blocks, coord, get_)
        sum :: Float64 = 0.0
        weight :: Float64 = 0.0
        for b in blocks
            if block_needed(coord, b, 3)
                # loop through elements
                @inbounds for k = 0:block_info.dim[3]-1
                    @inbounds for j = 0:block_info.dim[3]-1
                        @inbounds for i = 0:block_info.dim[3]-1
                            vrtx, value = get_(b, i,j,k)
                            #plot!(vrtx[1], vrtx[2])
                            @. δx = (vrtx - coord)/block_info.spacing
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

function integrand(get_,blocks_cap, l,m,r,τ,ϕ)

    θ = acos(τ)
    z = r * cos(θ)
    x = r * sin(θ)cos(ϕ)
    y = r * sin(θ)sin(ϕ)

    f = intrpl(blocks_cap, [x,y,z], get_)
    #println([x,y,z,f])

    return f * Y(l,m,τ,ϕ)
end

function vel_intrp_int_coeff(r, blocks_cap, l, m)
    get_v1 = (b, i, j, k)-> get_vel(b, i, j, k, 1)
    (val,err) =
    @suppress_err hcubature(
    x->integrand(get_v1, blocks_cap, l, m, r, x[1], x[2]),
                [-1, -π], [1, π];
                reltol=1e-4, maxevals = 300)
    return val, err
end
