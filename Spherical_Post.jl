using Pkg
using LinearAlgebra
import GSL
using Cubature
using Suppressor

include("vtk_read.jl")
using JLD2
include("sphere_utils.jl")
include("integrator.jl")
include("tests.jl")
include("ddf.jl")
include("coeff_methods.jl")
include("coeff_methods2.jl")


# Main  ------------------------------------------
const center = [18.8495559215, 18.8495559215, 18.8495559215]
const directory = pwd() * "/256-64-t0044/"#"/flow/"#"/256-64-t004/"

# read all
@time b, b_info = read_all(directory,center)
const block_info = b_info
const blocks = b

# spherical mesh
R = [3]
Nl = Int(round(2*pi/block_info.spacing[1] /2))

for r in R
    @time f, k, C, lC = method2(r, Nl)
    @save "$(Nl*2)-$(r).jld" k C lC
end
