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
