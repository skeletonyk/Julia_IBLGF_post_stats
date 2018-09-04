@inline function Y_θ(l, m, θ)
    if m>0
        return GSL.sf_legendre_sphPlm(l, m, cos(θ)) *sqrt(2) * sin(θ)
    elseif m==0
        return GSL.sf_legendre_sphPlm(l, m, cos(θ))* sin(θ)
    else
        return GSL.sf_legendre_sphPlm(l, -m, cos(θ))*sqrt(2) * sin(θ)
    end
end

@inline function Y_θϕ(m,ϕ)
    if m>=0
        return cos(m * ϕ)
    else
        return sin(-m * ϕ)
    end
end

#@inline function Y(l, m, θ, ϕ)
#    if m>0
#        return GSL.sf_legendre_sphPlm(l, m, θ) * cos(m * ϕ)*sqrt(2) sin
#    elseif m==0
#        return GSL.sf_legendre_sphPlm(l, m, θ) * cos(m * ϕ) sin
#    else
#        return GSL.sf_legendre_sphPlm(l, -m, θ) * sin(-m * ϕ)*sqrt(2) sin
#    end
#end

function integrate(f, g, N_θ, N_ϕ)
    sm = 0.0

    # cell #1 = cell #N_ϕ
    for j = 2:N_ϕ-1
        for i = 2:N_θ-1
            sm += f[i,j] * g[i,j] * 0.5
            sm += (f[i-1,j] + f[i-1,j-1] + f[i,j-1] + f[i+1, j] + f[i+1, j+1] +
            f[i, j+1]) * g[i,j] * (1.0/12)
        end
    end

    j = 1
    for i = 2:N_θ-1
        sm += f[i,j] * g[i,j] * 0.5
        sm += (f[i-1,j] + f[i+1, j] + f[i+1, j+1] + f[i, j+1]) * (1.0/12)
    end

    j = N_ϕ
    for i = 2:N_θ-1
        sm += (f[i-1,j-1] + f[i,j-1] ) * (1.0/12)
    end

    i = 1
    for j = 2:N_ϕ-1
        sm += f[i,j] * g[i,j] * 0.25
        sm += (f[i,j-1] + f[i+1, j] + f[i+1, j+1] +
        f[i, j+1]) * g[i,j] * (1.0/12)
    end

    i = N_θ
    for j = 2:N_ϕ-1
        sm += f[i,j] * g[i,j] * 0.25
        sm += (f[i-1,j] + f[i-1,j-1] + f[i,j-1]+
        f[i, j+1]) * g[i,j] * (1.0/12)
    end

    sm += f[1,1] * g[1,1] * (1/12)
    sm += f[N_θ,N_ϕ] * g[N_θ,N_ϕ] * (1/12)
    sm += f[1,N_ϕ] * g[1,N_ϕ] * (1/6)
    sm += f[N_θ,1] * g[N_θ,1] * (1/6)

    return sm
end

function harmonics_coefficients(N_k, s, r)
    println("calculating f_lm")

    N_l = N_k .* r

    f_lm ::Array{Float64,3} = zeros(N_k, 2 * N_l, 3) # for 3 vel vel_component
    Y_lm  ::Array{Float64,2} = zeros(s.N_θ, s.N_ϕ)
    product  ::Array{Float64,2} = zeros(s.N_θ, s.N_ϕ)

    y_θ =0.0

    for k = 0 : N_k-1
        l = k*r

        println("-- l = $(l) --")
        for m = -l : l
            # integrating
            for i = 1 : s.N_θ
                θ = s.mesh_θ[i]
                y_θ =Y_θ(l,m,θ)

                for j = 1 : s.N_ϕ
                    ϕ = s.mesh_ϕ[j]
                    Y_lm[i, j] = y_θ * Y_θϕ(m,ϕ)
                end

            end

            for vel_component = 1 :3
                #product .= view(s.u,:,:,vel_component) .* Y_lm
                intgl = integrate(view(s.u,:,:,vel_component), Y_lm, s.N_θ, s.N_ϕ)
                f_lm[k+1, m+l+1, vel_component] = intgl / s.N_θ / s.N_ϕ * 4 * π
            end
        end
    end

    println("calculating f_lm --- end")
    k = sqrt.((0:(N_k-1)) .* (1:(N_k)))
    return f_lm, k
end
