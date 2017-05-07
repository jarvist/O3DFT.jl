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

"""
ThomasFermi_T(n)

Thomas Fermi kinetic energy (T).
 Following eqn. 9.69 Marder p/ 217
"""
function ThomasFermi_T(n)
    T= (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(5/3)
    T
end

"
ThomasFermi_T_fnderiv(n)

Thomas Fermi kinetic energy (T) as a functional derivative.
 Following eqn. 9.76 Marder p/ 217
"
function ThomasFermi_T_fnderiv(n)
#    @printf("ThomasFermi_T_fnderiv(%g,%g)\n",n,V) # to locate NaN error
    T= (hbar^2)/(2*me) * 3/5 * (3*pi^2)^(2/3) * n^(2/3)
    T
end

" 
 ThomasFermi_Exc(n)

 Thomas Fermi exchange and correlation energy; takes electron density, returns energy.
 Following eqn. 9.73 Marder p/ 217
"
function ThomasFermi_Exc(n)
    Exc= - 3/4 * (3/pi)^(1/3) * q^2 * n^(5/3)
    Exc
end

" 
UAtomic(Z,r)

Potential due to Atomic charge; simple bare Coulomb form. 
"
function UAtomic(Z,r)
    U = -k_e * Z * q^2/r
end

function ThomasFermi_CoulombPotential(density,i,gridspacing,N)
    # Coulomb integral. 
    # Nb: as spherical atom, probably need some horrible weighting parameters for the area of the spherical shell in the range [i] or [j]
    # i.e. currently it's some 1D pipe full of electrons with a nuclear charge at the end
    Potential=0.0
    for j in 1:N
        if j==i 
            continue 
        end
        Potential+= k_e * q^2*density[j]/(gridspacing*(i-j))
    end
    Potential
end


"
 AtomicTest(Z,N=100)
 Evaluate Thomas Fermi energy for a spherical atom, atomic charge Z.
 N: number of grid points in 1D density repr.

 c.f. Marder, p.219, Eqn (9.76):
 'The energy of an atomic of nuclear charge Z is approximately -1.5375 Z^(7/3) Ry'
"
function AtomicTest(Z; N=10,verbose::Bool=false)
    println("AtomicTest, atomic charge: ",Z)

    const radius=25Å

    density=zeros(N).+Z/N # I know, a bit dirty... Fills density with flat electron density as initial guess.
    gridspacing=N/radius

    # Do some DFT here...
    V = 4/3*pi*radius^3 # Volume of total sphere of charge density considered. 
    # Nb: Not sure if corect; Kinetic energy T is proport. to V. Defo volume and not potential?

    sumE=0.0
    for n in 1:10
        if verbose; @printf("SCF Loop: %d\n",n); end;
        sumE=0.0
        # Central equation as (9.76) in Marder, dropping Dirac term
        for i in 1:N

            # mu being the chemical potential; this pulls together (9.76)
            # Mu is the Lagrange multiplier to enforce the constraint that the density is conserved; 
            # helpfully this is also identical to the chemical potential
            mu=ThomasFermi_T_fnderiv(density[i]) + UAtomic(Z,i*gridspacing) + ThomasFermi_CoulombPotential(density,i,gridspacing,N)
            # So I gues we could use the fact that the chem potential is constant, to start moving electron density around?

            # From Postnikov 1.2; mu drops to zero at r=infinity, therefore mu is zero everywhere

            # TODO: Insert self-consistency here...

            if verbose
                @printf("\t i: %d density: %g T: %g \n\t\tmu %g = T_fnderiv %g + UAtomic: %g + Coulomb %g\n",
                i,density[i],
                ThomasFermi_T(density[i]),
                mu,ThomasFermi_T_fnderiv(density[i]),UAtomic(Z,i*gridspacing),ThomasFermi_CoulombPotential(density,i,gridspacing,N))
            end

            # OK; calculate total energy
            E=ThomasFermi_T(density[i]) + density[i]*UAtomic(Z,i*gridspacing) + density[i]*ThomasFermi_CoulombPotential(density,i,gridspacing,N)
            if verbose
                @printf("\t\tE %g = T %g + U %g + Coulomb %g\nTotal E: %g J = %g eV\n",
                E,ThomasFermi_T(density[i]), density[i]*UAtomic(Z,i*gridspacing), density[i]*ThomasFermi_CoulombPotential(density,i,gridspacing,N),
                sumE,sumE/q)
            end
            sumE+=E

            # Nb: horrid hack :^)
            density[i]-=mu*10E35 # vary density based on how much chemical potential mu exceeds 0
            if density[i]<0.0; density[i]=0.0; end
        end
        # Impose constraint sum. density = Z
        if verbose
            @printf("Sum of density pre normalisation: %f\n",sum(density))
        end
        density=density.*Z/sum(density)
 
    end

    # TODO: calculate total Thomas-Fermi energy here, from density...

    println("Density: ",density)
    @printf("Total E: %g J = %g Ry = %g eV\n",sumE,sumE/Ry,sumE/q)
    println("Marder analytic reference: E ~= -1.5375.Z^(7/3) = ",-1.5375*Z^(7/3), " Ry") # Nb: do some unit conversions
end

function main()
    println("Orbital 3 Density Functional Theory - Philsophy by Numbers")

    # Does anyone know what units we should be in?
    for n in 1E25:1E25:1E26 # density /m^3
        @printf("n: %g \t T(n): %g\t T_fnderiv(n): %g\t E_xc(n): %g\n",n,ThomasFermi_T(n),ThomasFermi_T_fnderiv(n),ThomasFermi_Exc(n))
    end

    for n=1:5:20
        AtomicTest(n,verbose=false) #Oh uh huh make it magnificent
    end
end

main() # C-party like the mid-90s
