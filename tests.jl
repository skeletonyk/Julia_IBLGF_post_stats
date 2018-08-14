function intrp_test(r, blocks_cap)
    for i = 1 : 1
        θ = rand()*pi
        ϕ = rand()*2*pi - pi
        z = r * cos(θ)
        x = r * sin(θ)cos(ϕ)
        y = r * sin(θ)sin(ϕ)
        intrpl(blocks_cap, [x,y,z], get_vel)
    end
end

function integral_test()
    (val1,err) =
    @suppress hcubature( x -> Y(1, 0, x[1], x[2]),
                [-1, -π], [1, π];
                reltol=1e-8)

    (val2,err) =
    @suppress hcubature( x -> (Y(2, 1, x[1], x[2]))^2,
                [-1, -π], [1, π];
                reltol=1e-8)

    (val3,err) =
    @suppress hcubature( x -> (Y(2, 1, x[1], x[2]))^2,
                [-1, -π], [1, π];
                reltol=1e-8)

    (val4,err) =
    @suppress hcubature( x -> (Y(2, 0, x[1], x[2]))^2,
                [-1, -π], [1, π];
                reltol=1e-8)

    println([val1, val2, val3, val4])
    nothing

end
