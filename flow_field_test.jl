field = zeros(400,400)

for bl in blocks
    for k = 0:block_info.dim[3]-1
        for j = 0:block_info.dim[2]-1
            for i = 0:block_info.dim[1]-1
                vrtx, value = get_vel(bl, i, j, k, 1)
                if abs(vrtx[3]) < 1e-3
                    tmp = Int64.(round.((vrtx+center)./block_info.spacing))
                    println(value)
                    field[tmp[1], tmp[2]] = value
                end
            end
        end
    end
end
contour(field,fill=true)
savefig("tau-phi-64.png")
savefig("tau-phi-64.pdf")
