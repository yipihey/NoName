using Cosmology
using Roots

global const gyr_in_s = 1.e9*31_556_952

# recover time from redshift by rootfinding on the age function
function z_from_t(cosmo, t)
    # this function will not work if we simulate the Universe for ages older than 500 billion years. So 35 times into the future... (not likely)
    res = fzero(x -> age_gyr(cosmo, x) - t, [-1.+1e-16, 1e90]) # between 0 and 1e3 Gyr ...
    return res
end

function a_from_t(cosmo, t; zwherea1=0.0)
    z = z_from_t(cosmo, t)
    a = (1.0 + zwherea1)/(1.0 + z)
    return a
end

function a_from_t_internal(cosmo, t, zinit; zwherea1=0.0)
    u = get_units(cosmo, 0.0, zinit=zinit)
    time = t*(u.T/gyr_in_s)
    a_from_t(cosmo, time, zwherea1=zwherea1)
end

function z_from_t_internal(cosmo, t, zinit; zwherea1=0.0)
    a = a_from_t_internal(cosmo, t, zinit)
    z_from_a(cosmo, a)
end

function dadt(cosmo, a; zwherea1=0.0)
    tv = a/(1.0 + zwherea1) 
    O_k = 1.-cosmo.Ω_m-cosmo.Ω_Λ-cosmo.Ω_r # curvature today
    dadt = sqrt(2./(3.*cosmo.Ω_m * a) *    # (Peebles93, eq. 13.3)
                (cosmo.Ω_m + O_k * tv +
                 cosmo.Ω_Λ * tv^3))
    return dadt
end

function dadt_from_t_internal(cosmo, t, zinit; zhwerea1=0.0)
    a = a_from_t_internal(cosmo, t, zinit, zwherea1=zinit)
    dadt(cosmo, a; zwherea1=zinit)
end

function z_from_a(cosmo, a; zwherea1=0.0) 
    (1.0 + zwherea1)/a - 1.0
end

type CodeUnits
    ρ
    L
    T
    Temp
    V
end

doc"""
Follow the unit system used in Enzo: http://enzo-project.org
"""
function get_units(cosmo, z; zinit=0., boxsize=1.0)
#    a = a_from_t(cosmo, t, zwherea1=zinit)
    CurrentRedshift = z # z_from_a(cosmo, a, zwherea1=zinit)
    u      = CodeUnits(1.0, 1.0, 1.0, 1.0, 1.0)
    u.ρ    = 1.8788e-29*cosmo.Ω_m*cosmo.h^2*(1.0 + CurrentRedshift)^3
    u.L    = 3.085678e24*boxsize/cosmo.h/(1.0 + CurrentRedshift)
    # temperature unit assumes mean molecular weight μ = 1 
    u.Temp = 1.81723e6*boxsize^2 * cosmo.Ω_m * (1 + zinit) 
    u.T    = 2.519445e17/sqrt(cosmo.Ω_m)/cosmo.h/(1 + zinit)^1.5
    u.V    = 1.22475e7*boxsize*sqrt(cosmo.Ω_m)*sqrt(1 + zinit)
    return u
end
