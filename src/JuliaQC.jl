# Description: This file contains the main module for the JuliaQC package.
# The JuliaQC package is a quantum chemistry package written in Julia.
# Author: Chunyu Yang
# Email: chunyu.yang@duke.edu, chyyangustc@outlook.com (permanent)
# Organization: Department of Chemistry, Duke University
# Date: 2024-12-13



module JuliaQC

include("InputIntrepreter.jl")
using .InputIntrepreter

include("HF.jl")
using .HF

include("Integral.jl")
using .Integral

include("CC.jl")
using .CC

include("CI.jl")
using .CI

include("MP.jl")
using .MP

include("DFT.jl")
using .DFT


include("QCMath.jl")
using .QCMath

include("TDDFT.jl")
using .TDDFT

include("GW.jl")
using .GW


end # module JuliaQC
