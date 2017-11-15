# O3DFT.jl
#  - Orbital 3 Density Functional Theory

push!(LOAD_PATH,"../src/") # load module from local directory
using O3DFT 


function main()
    # Does anyone know what units we should be in?
    for n in 1E25:1E25:1E26 # density /m^3
        @printf("n: %g \t T(n): %g\t T_fnderiv(n): %g\t E_xc(n): %g\n",n,ThomasFermi_T(n),ThomasFermi_T_fnderiv(n),ThomasFermi_Exc(n))
    end

    for n=1:5:20
        AtomicTest(n,verbose=false) #Oh uh huh make it magnificent
    end
end

main() # C-party like the mid-90s
