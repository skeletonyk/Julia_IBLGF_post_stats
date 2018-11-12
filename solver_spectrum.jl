using Distributions
using MarketTechnicals

function spectrum(K, Nk, f_lm, r)

    Nl = Nk * r

    C = zeros(Nk,3)
    lC = zeros(Nk,3)

    for vel_component = 1 : 3
        for k = 0:Nk-1
            l = k*r
            for m = -l : l
                C[k+1, vel_component] +=( f_lm[k+1, m+l+1, vel_component])^2
            end

            C[k+1, vel_component] = C[k+1, vel_component]/(2*l+1)
            lC[k+1, vel_component] = l*C[k+1, vel_component]
        end
    end

    smooth = 2
    #lC = sma(lC, smooth)
    C = sma(C, smooth)
    #K = sma(collect(K), smooth)

    α = 5/3
    gamma = GSL.sf_gamma
    A = -2^α *( gamma( - α) * gamma( (1+α)/2) ) / (gamma((1 - α)/2) ) * sin(π/2 * α) /r
    C = C/A
    lC = lC/A
    return K, C, lC
end
