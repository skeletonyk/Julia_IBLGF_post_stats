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

function yang4(r)
# Smoothed 4-point delta function of Yang et al (2009)
    rr = abs(r)
    if rr<=0.5
        return 3/8+pi/32-rr^2/4
    elseif (rr>0.5) && (rr<=1.5)
        return 1/4+(1-rr)/8*sqrt(-2+8*rr-4*rr^2)-1/8*asin(sqrt(2)*(rr-1))
    elseif (rr>1.5) && (rr<=2.5)
        return 17/16-pi/64-3/4*rr+rr^2/8+(rr-2)/16*sqrt(-14+16*rr-4*rr^2)+
             1/16*asin(sqrt(2)*(rr-2))
    else return 0.0
    end

end

function ddf(δx)
    #return yang3(δx[1])*yang3(δx[2])*yang3(δx[3])
    return yang4(δx[1])*yang4(δx[2])*yang4(δx[3])
end
