struct block_t
    vel    :: Array{Float64,2}
    origin :: Array{Float64,1}
    r :: Float64
end

struct blocks_info_t
    spacing :: Array{Float64,1}
    dim :: Array{Int64,1}
    width ::Float64
end
