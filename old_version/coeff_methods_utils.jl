function flm_2_k_C(Nl, f_lm)
    C = zeros(Nl)
    lC = zeros(Nl)
    for l = 0:Nl-1
        for m = -l : l
            C[l+1] +=f_lm[l+1, m+l+1]^2
        end
        C[l+1] = C[l+1]/(2*l+1)
        lC[l+1] = l*C[l+1]
    end
    k = sqrt.((0:(Nl-1)) .* (1:(Nl)))
    return k, C, lC
end
