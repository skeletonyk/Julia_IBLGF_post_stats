@inline function appendSpherical(coord)
    @. xy .= coord[:,1].^2 + coord[:,2].^2
    @. coord[:,4] = sqrt.(xy + coord[:,3].^2)
    @. coord[:,5] = atan.(sqrt.(xy), coord[:,3]) # for elevation angle defined from Z-axis down
    @. coord[:,6] = atan.(coord[:,2], coord[:,1])

    nothing
end

function find_left(blocks, r)
    i = 1
    while blocks[i].r < r - block_info.width * sqrt(3)*1.1
        i+=1
    end
    return i
end

function find_right(blocks, r)

    i = size(blocks, 1)
    while blocks[i].r > r + block_info.width *sqrt(3)*1.1
        i =i-1
    end
    return i
end
