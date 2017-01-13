# O3DFT.jl
#  - Orbital 3 Density Functional Theory

# Marder is: "Condensed Mattery Physics", Michael Marder, 2000. 

# Also very useful is Andrei Postnikov's first DFT lecture, on the Thomas-Fermi method 
# http://www.home.uni-osnabrueck.de/apostnik/Lectures/DFT-1.pdf

# Constants:
const hbar = 1.054E-34
const h =    6.62606957E-34
const kb =   8.6173324E-5 # in units of eV
const ε_0 =  8.854E-12 #Units: C2N−1m−2, permittivity of free space
const q = eV = 1.60217657E-19
const c = 3E8 # Speed of light in this universe
const k_e = 8.98755E+9 # Coulomb's constant
const me = 9.10938291E-31 # Mass of electron (kg)
const rbohr = 5.2917721067E-11 # Bohr radius (m)

const Ha=27.2111eV
const Ry=Ha/2

const Å=1E-10 # Angstrom

function ThomasFermi_T(n)
    # Following eqn. 9.69 Marder p/ 217
    T= (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(5/3)
    T
end

function ThomasFermi_T_fnderiv(n)
    # Following eqn. 9.76 Marder p/ 217
#    @printf("ThomasFermi_T_fnderiv(%g,%g)\n",n,V) # to locate NaN error
    T= (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(2/3)
    T
end

function ThomasFermi_Exc(n)
    # Following eqn. 9.73 Marder p/ 217
    Exc= - 3/4 * (3/pi)^(1/3) * q^2 * n^(5/3)
    Exc
end

" Potential due to Atomic charge; simple Coulomb form. "
function UAtomic(Z,r)
    U = -k_e * Z * q^2/r
end


"
 AtomicTest(Z,N=100)
 Evaluate Thomas Fermi energy for a spherical atom, atomic charge Z.
 N: number of grid points in 1D density repr.

 c.f. Marder, p.219, Eqn (9.76):
 'The energy of an atomic of nuclear charge Z is approximately -1.5375 Z^(7/3) Ry'
"
function AtomicTest(Z,N=10)
    println("AtomicTest, atomic charge: ",Z)

    const radius=100Å

    density=zeros(N).+Z/N # I know, a bit dirty... Fills density with flat electron density as initial guess.
    gridspacing=N/radius

    # Do some DFT here...
    V = 4/3*pi*radius^3 # Volume of total sphere of charge density considered. 
    # Nb: Not sure if corect; Kinetic energy T is proport. to V. Defo volume and not potential?

    for n in 1:10
        @printf("SCF Loop: %d\n",n)
        # Central equation as (9.76) in Marder, dropping Dirac term
        for i in 1:N

            # Coulomb integral. 
            # Nb: as spherical atom, probably need some horrible weighting parameters for the area of the spherical shell in the range [i] or [j]
            # i.e. currently it's some 1D pipe full of electrons with a nuclear charge at the end
            ThomasFermi_Coulomb=0.0
            for j in 1:N
                if j==i 
                    continue 
                end
                ThomasFermi_Coulomb+= k_e * q^2*density[j]/(gridspacing*(i-j))
            end

            # mu being the chemical potential; this pulls together (9.76)
            # Mu is the Lagrange multiplier to enforce the constraint that the density is conserved; 
            # helpfully this is also identical to the chemical potential
            mu=ThomasFermi_T_fnderiv(density[i]) + UAtomic(Z,i*gridspacing) + ThomasFermi_Coulomb
            # So I gues we could use the fact that the chem potential is constant, to start moving electron density around?

            # From Postnikov 1.2; mu drops to zero at r=infinity, therefore mu is zero everywhere

            # TODO: Insert self-consistency here...

            @printf("\t i: %d density: %g T: %g \n\t\tmu %g = T_fnderiv %g + UAtomic: %g + Coulomb %g\n",
            i,density[i],
            ThomasFermi_T(density[i]),
            mu,ThomasFermi_T_fnderiv(density[i]),UAtomic(Z,i*gridspacing),ThomasFermi_Coulomb)

            # Nb: horrid hack :^)
            density[i]-=mu*10E35 # vary density based on how much chemical potential mu exceeds 0
            if density[i]<0.0; density[i]=0.0; end
        end
        # Impose constraint sum. density = Z
        @printf("Sum of density pre normalisatio: %f\n",sum(density))
        density=density.*Z/sum(density)
 
    end

    # TODO: calculate total Thomas-Fermi energy here, from density...

    println("Density: ",density)
    println("Marder: E ~= -1.5375.Z^(7/3) = ",-1.5375*Z^(7/3), " Ry") # Nb: do some unit conversions
end

function main()
    println("Orbital 3 Density Functional Theory - Philsophy by Numbers")

    # Does anyone know what units we should be in?
    for n in 1E25:1E25:1E26 # density /m^3
        @printf("n: %g \t T(n): %g\t T_fnderiv(n): %g\t E_xc(n): %g\n",n,ThomasFermi_T(n),ThomasFermi_T_fnderiv(n),ThomasFermi_Exc(n))
    end

    for n=9:10
        AtomicTest(n) #Oh uh huh make it magnificent
    end
end

main() # C-party like the mid-90s
