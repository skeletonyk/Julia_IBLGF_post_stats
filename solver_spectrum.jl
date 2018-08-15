function spectrum(Nl, f_lm)
    C = zeros(Nl,3)
    lC = zeros(Nl,3)

    for vel_component = 1 : 3
        for l = 0:Nl-1
            for m = -l : l
                C[l+1, vel_component] +=f_lm[l+1, m+l+1, vel_component]^2
            end

            C[l+1, vel_component] = C[l+1, vel_component]/(2*l+1)
            lC[l+1, vel_component] = l*C[l+1, vel_component]
        end
    end

    k = sqrt.((0:(Nl-1)) .* (1:(Nl)))
    return k, C, lC
end
