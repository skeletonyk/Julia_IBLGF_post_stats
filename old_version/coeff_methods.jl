include("coeff_methods_utils.jl")

function adaptive_integration_method(Nl, R)
    f_lm = zeros(Nl,2*Nl)
    for r in R
        println("working on spherical cap of radius = ", r)
        # find needed blocks

        lb=find_left(blocks, r)
        ub=find_right(blocks, r)
        blocks_cap = view(blocks, lb:ub)

        # interpolation sanity check
        #intrp_test(r,blocks_cap)
        # multidimensional integral sanity check
        #@suppress_err integral_test()
        for l = 0:Nl-1
            for m = -l : l
                val, err = vel_intrp_int_coeff(r, blocks_cap, l, m)
                f_lm[l+1,m+l+1] = val
                println([l,m,val,err])
            end
        end
        
    end

    k, C = flm_2_k_C(Nl, f_lm)
    return f_lm, C, k
end
