function Y(l, m, τ, ϕ)
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
    f_lm = zeros(N_l, 2 * N_l)
    Y_lm = zeros(s.N_τ, s.N_ϕ)

    for l = 0 : N_l-1
        for m = -l : l

            # integrating
            for i = 1 : s.N_τ
                for j = 1 : s.N_ϕ
                    τ = s.mesh_τ[i]
                    ϕ = s.mesh_ϕ[j]
                    # cal x y z for interpolation
                    θ = acos(τ)
                    Y_lm[i, j] = Y(l,m,τ,ϕ)
                end
            end

            Y_lm .*= s.f
            f_lm[l+1,m+l+1] = sum(Y_lm ) / s.N_τ / s.N_ϕ * 4 * π

        end
    end

    return f_lm
end
