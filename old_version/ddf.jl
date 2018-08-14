using MuladdMacro
function yang3(r)
    rr = abs(r)
    if rr<=1.0
        return @muladd 17/48 + sqrt(3) * pi/108 + rr/4 -
            rr^2/4 + (1-2*rr)/16 * sqrt(-12*rr^2+12*rr+1) -
            sqrt(3)/12 * asin(sqrt(3)/2*(2*rr-1))
    elseif (rr>1.0) && (rr<=2.0)
        return @muladd 55/48-sqrt(3)*pi/108-13*rr/12+
            rr^2/4+(2*rr-3)/48.0*sqrt(-12*rr.^2+36*rr-23) +
            sqrt(3)/36 * asin(sqrt(3)/2*(2*rr-3))
    else return 0.0
    end
end

function ddf(δx)
    return yang3(δx[1])*yang3(δx[2])*yang3(δx[3])
end
