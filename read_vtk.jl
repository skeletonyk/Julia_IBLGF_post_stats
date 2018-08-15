#ENV["PYTHON"] = "/Users/yuke1/anaconda2/bin/python"
#ENV["PYTHON"] = "/home/ke/anaconda3/bin/python"
using PyCall

@pyimport read_node

function read_all(directory, center)
    list = filter( x->endswith(x,"vti"), readdir(directory))
    f_count = length(list)
    #coord_all = zeros(size(coord,1)*f_count, 6)

    blocks = Array{block_t, 1}(undef, f_count)

    i = 1
    spacing = []
    dim = []
    dlist :: Array{String,1}= readdir(directory)
    for fn in dlist
        if endswith(fn, "vti")
            #println(directory * fn)
            vel, origin, spacing, dim = read_node.vtiread(directory * fn)
            origin -= center
            r = norm(origin)
            blocks[i]= block_t(vel, origin, r)
            #n = size(coord,1)
            #coord_all[n_start : n+n_start-1, 1:3] .= coord
            i += 1
        end
    end

    blocks_info = blocks_info_t(spacing, dim, (dim[1]+1)*spacing[1])
    sort!(blocks, by=x->x.r)

    return blocks, blocks_info
end
