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
const directory = pwd() * "/256-64-t0044/"#"/flow/"#"/256-64-t004/"

# ---- read all --------------------------------------------------------------
#

println("- reading the file -")
@time b, b_info = read_all(directory,center)
const block_info = b_info
const blocks = b
println("- reading the file - end")

# ---- main ------------------------------------------------------------------
#

println()
println("- main - ")
# for all radius choices
R = [3]

for r in R
    @time f, k, C, lC = stats_shell(r)
    @save pwd()*"/data/$(64)-$(r).jld" k C lC
end

println("- main - end")


# ---- end ------------------------------------------------------------------
