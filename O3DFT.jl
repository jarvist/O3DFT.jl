# O3DFT.jl
#  - Orbital 3 Density Functional Theory

# Marder is: "Condensed Mattery Physics", Michael Marder, 2000. 

# Constants:
const hbar = 1.054E-34
const h =    6.62606957E-34
const kb =   8.6173324E-5 # in units of eV
const ε_0 =  8.854E-12 #Units: C2N−1m−2, permittivity of free space
const q = eV = 1.60217657E-19
const c = 3E8 # Speed of light in this universe
const me = 9.10938291E-31 # Mass of electron (kg)
const Å=1E-10 # Angstrom

function ThomasFermi_T(n,V)
    # Following eqn. 9.69 Marder p/ 217
    T=V * (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(5/3)
    T
end

function ThomasFermi_T_fnderiv(n,V)
    # Following eqn. 9.76 Marder p/ 217
    T=V * (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(2/3)
    T
end

function ThomasFermi_Exc(n,V)
    # Following eqn. 9.73 Marder p/ 217
    Exc= - 3/4 * (3/pi)^(1/3) * q^2 * n^(5/3)
    Exc
end


function main()
    println("Orbital 3 Density Functional Theory - Philsophy by Numbers")

    const V = 1.0
    # Does anyone know what units we should be in?
    for n in 1E25:1E25:1E26 # density /m^3
        @printf("n: %g \t T(n): %g\t T_fnderiv(n): %g\t E_xc(n): %g\n",n,ThomasFermi_T(n,V),ThomasFermi_T_fnderiv(n,V),ThomasFermi_Exc(n,V))
    end
end

main() # C-party like the mid-90s
