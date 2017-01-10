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
const rbohr = 5.2917721067E-11 # Bohr radius (m)

const Ha=27.2111eV
const Ry=Ha/2

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

"
 AtomicTest(Z,N=100)
 Evaluate Thomas Fermi energy for a spherical atom, atomic charge Z.
 N: number of grid points in 1D density repr.

 c.f. Marder, p.219, Eqn (9.76):
 'The energy of an atomic of nuclear charge Z is approximately -1.5375 Z^(7/3) Ry'
"
function AtomicTest(Z,N=100)
    println("AtomicTest, atomic charge: ",Z)
    density=zeros(N).+Z/N # I know, a bit dirty... Fills density with flat electron density as initial guess.
    gridspacing=N/10Å

    # Do some DFT here...
    V = 0.01 # OK; should be some spherical statement here? In shell? Or volume?

    # Central equation as (9.76) in Marder, dropping Dirac term
    for i in 1:N

        # Coulomb integral
        ThomasFermi_Coulomb=0.0
        for j in 1:N
            if j==i 
                continue 
            end
            ThomasFermi_Coulomb+= q^2*density[j]/(gridspacing*(i-j))
        end

        mu=ThomasFermi_T_fnderiv(density[i],V) + UAtomic(Z,i*gridspacing) + ThomasFermi_Coulomb

        @printf("\t i: %d mu: %g = T_fnderiv %g + UAtomic: %g + Coulomb %g\n",
                i,mu,ThomasFermi_T_fnderiv(density[i],V),UAtomic(Z,i*gridspacing),ThomasFermi_Coulomb)
    end

    # Uhm, calculate total Thomas-Fermi energy here, from density...

    println("Density: ",density)
    println("Marder: E ~= -1.5375.Z^(7/3) = ",-1.5375*Z^(7/3), " Ry") # Nb: do some unit conversions
end

function UAtomic(Z,r)
    U=-Z*q^2/r
end

function main()
    println("Orbital 3 Density Functional Theory - Philsophy by Numbers")

    const V = 1.0
    # Does anyone know what units we should be in?
    for n in 1E25:1E25:1E26 # density /m^3
        @printf("n: %g \t T(n): %g\t T_fnderiv(n): %g\t E_xc(n): %g\n",n,ThomasFermi_T(n,V),ThomasFermi_T_fnderiv(n,V),ThomasFermi_Exc(n,V))
    end

    for n=1:10
        AtomicTest(n) #Oh uh huh make it magnificent
    end
end

main() # C-party like the mid-90s
