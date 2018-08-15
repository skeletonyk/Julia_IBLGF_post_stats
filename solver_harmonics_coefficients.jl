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

function integrate(f, g, N_τ, N_ϕ)
    sm = 0.0

    # cell #1 = cell #N_ϕ
    for j = 2:N_ϕ-1
        for i = 2:N_τ-1
            sm += f[i,j] * g[i,j] * 0.5
            sm += (f[i-1,j] + f[i-1,j-1] + f[i,j-1] + f[i+1, j] + f[i+1, j+1] +
            f[i, j+1]) * g[i,j] * (1.0/12)
        end
    end

    j = 1
    for i = 2:N_τ-1
        sm += f[i,j] * g[i,j] * 0.5
        sm += (f[i-1,j] + f[i+1, j] + f[i+1, j+1] + f[i, j+1]) * (1.0/12)
    end

    j = N_ϕ
    for i = 2:N_τ-1
        sm += (f[i-1,j-1] + f[i,j-1] ) * (1.0/12)
    end

    i = 1
    for j = 2:N_ϕ-1
        sm += f[i,j] * g[i,j] * 0.25
        sm += (f[i,j-1] + f[i+1, j] + f[i+1, j+1] +
        f[i, j+1]) * g[i,j] * (1.0/12)
    end

    i = N_τ
    for j = 2:N_ϕ-1
        sm += f[i,j] * g[i,j] * 0.25
        sm += (f[i-1,j] + f[i-1,j-1] + f[i,j-1]+
        f[i, j+1]) * g[i,j] * (1.0/12)
    end

    sm += f[1,1] * g[1,1] * (1/12)
    sm += f[N_τ,N_ϕ] * g[N_τ,N_ϕ] * (1/12)
    sm += f[1,N_ϕ] * g[1,N_ϕ] * (1/6)
    sm += f[N_τ,1] * g[N_τ,1] * (1/6)

    return sm
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

                #product .= view(s.u,:,:,vel_component) .* Y_lm
                intgl = integrate(view(s.u,:,:,vel_component), Y_lm, s.N_τ, s.N_ϕ)
                f_lm[l+1, m+l+1, vel_component] = intgl / s.N_τ / s.N_ϕ * 4 * π
            end

        end
    end

    return f_lm
end
