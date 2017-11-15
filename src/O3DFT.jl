# O3DFT.jl
#  - Orbital 3 Density Functional Theory

# Marder is: "Condensed Mattery Physics", Michael Marder, 2000. 

# Also very useful is Andrei Postnikov's first DFT lecture, on the Thomas-Fermi method 
# http://www.home.uni-osnabrueck.de/apostnik/Lectures/DFT-1.pdf

module O3DFT

println("Orbital 3 Density Functional Theory - Philsophy by Numbers")

include("constants.jl")
include("ThomasFermi.jl")

export ThomasFermi_T, ThomasFermi_T_fnderiv, ThomasFermi_Exc, UAtomic
export ThomasFermi_CoulombPotential
export AtomicTest


end

