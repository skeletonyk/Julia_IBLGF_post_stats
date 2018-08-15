@inline function Y_τ(l, m, τ)
    if m>0
        return GSL.sf_legendre_sphPlm(l, m, τ) *sqrt(2)
    elseif m==0
        return GSL.sf_legendre_sphPlm(l, m, τ)
    else
        return GSL.sf_legendre_sphPlm(l, -m, τ)*sqrt(2)
    end
end

@inline function Y_τϕ(m,ϕ)
    if m>=0
        return cos(m * ϕ)
    else
        return sin(-m * ϕ)
    end
end

@inline function Y(l, m, τ, ϕ)
    if m>0
        return GSL.sf_legendre_sphPlm(l, m, τ) * cos(m * ϕ)*sqrt(2)
    elseif m==0
        return GSL.sf_legendre_sphPlm(l, m, τ) * cos(m * ϕ)
    else
        return GSL.sf_legendre_sphPlm(l, -m, τ) * sin(-m * ϕ)*sqrt(2)
    end
end

function harmonics_coefficients(N_l, s)
    println("calculating f_lm")
    f_lm ::Array{Float64,3} = zeros(N_l, 2 * N_l, 3) # for 3 vel vel_component
    Y_lm  ::Array{Float64,2} = zeros(s.N_τ, s.N_ϕ)
    product  ::Array{Float64,2} = zeros(s.N_τ, s.N_ϕ)

    y_τ =0.0
    for l = 0 : N_l-1

        for m = -l : l

            # integrating
            for i = 1 : s.N_τ
                τ = s.mesh_τ[i]
                y_τ =Y_τ(l,m,τ)

                for j = 1 : s.N_ϕ
                    ϕ = s.mesh_ϕ[j]
                    # cal x y z for interpolation
                    #θ = acos(τ)

                    Y_lm[i, j] = y_τ * Y_τϕ(m,ϕ)
                end

            end

            for vel_component = 1 :3

                product .= view(s.u,:,:,vel_component) .* Y_lm
                f_lm[l+1, m+l+1, vel_component] = sum(product) / s.N_τ / s.N_ϕ * 4 * π
            end

        end
    end

    return f_lm
end
