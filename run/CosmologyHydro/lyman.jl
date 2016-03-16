using HydroEOS, GridTypes, PyPlot
include("cosmology.jl")

λ_from_z(z) = 1216.*(z+1)

function LineOfSight(gp, gd, p, z, cosmo) # only works for 2D case now
    g = gp[1]
    d = gd[1]
    ρ  = d.d["ρ"] 
    v1 = d.d["V₁"]
    v2 = d.d["V₂"]
    u = get_units(cosmo, z)

    size = 500
    db = zeros(2,size)
    τ = zeros(size)
    velocity = zeros(size)
    loc = rand(5:g.dim[2]-5)
    mean_rho = mean(ρ)
    for i=1:500
        x = 2*i+loc
        y = pi*i
        db[1,i] = a_from_t(cosmo, age_gyr(cosmo,0)+sqrt(4+pi^2)*i*u.L/3e10/gyr_in_s)
        db[2,i] = ρ[x%g.dim[1]+1,floor(Int64,y%g.dim[1])+1,1]/mean_rho - 1
        vx = v1[x%g.dim[1]+1,floor(Int64,y%g.dim[1])+1,1]
        vy = v2[x%g.dim[1]+1,floor(Int64,y%g.dim[1])+1,1]
        velocity[i] = (dadt(cosmo,db[1,i])*sqrt(4+pi^2)*i+vx*2/sqrt(4+pi^2)+vy*pi/sqrt(4+pi^2))*u.V/1e5/g.dim[1] # note peculiar velocity is added
        τ[i] = τ_from_delta_b(cosmo,db[2,i], db[1,i])
    end
    return velocity, τ, db
end

function τ_from_delta_b(cosmo, db, a)
	Ω_b = 0.2 * cosmo.Ω_m*cosmo.h^2
	T0 = 8300. #8300, 11200
	Γ_HI = 1.5e-12
	z = 1./a -1
	H = dadt(cosmo,a)^2/a^2
	coeff = 0.433*((1)/3.5)^6*(Ω_b/0.02)^2*
			(T0/6000)^(-0.7)/(cosmo.h/0.65)/(H/cosmo.h/100/3.68)/(Γ_HI/1.5e-12)
	return coeff*(1+db)^1.7
end

v, t, db = LineOfSight(gp, gd, p, 3, NoName.cosmo)
writedlm("velocity100GeV.dat",v)
writedlm("tau100GeV.dat",t)
plot(v, e.^(-t))
#=
fig, ax = subplots()
ax[:plot](l,e.^(-t))
title("Lyman-α forest")
xlabel("λ")
ylabel("Flux = exp(-τ)")
ax[:legend]()=#