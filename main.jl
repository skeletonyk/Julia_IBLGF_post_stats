using Pkg
using LinearAlgebra
import GSL
using Cubature
using Suppressor
using JLD2

include("block_type.jl")
include("read_vtk.jl")
include("solver.jl")


const center = [18.8495559215, 18.8495559215, 18.8495559215]
const directory = pwd() * "/256-64-t0044/" #"/flow/"# "/256-64-t0044/"# "/home/kyu/scratch/LES-data-init/comet-256/post-processing/flow/tstep-0000002450/"# pwd() * "/256-64-t0044/"#"/flow/"#"/256-64-t0044/"#

# ---- read all --------------------------------------------------------------
#

println("- reading the file -")
@time b, b_info = read_all(directory,center)
const block_info = b_info
const blocks = b
const N = Int64(round(2*pi/block_info.spacing[1]))
println("- reading the file - end")

# ---- main ------------------------------------------------------------------
#

println()
println("- main - ")
# for all radius choices
R = [3]

for r in R
    @time f, k, C, lC = stats_shell(r)
    @save pwd()*"/data/$(N)-$(r).jld" k C lC
end

println("- main - end")


# ---- end ------------------------------------------------------------------
