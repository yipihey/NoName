## Support for cosmological calculations
3/2016-

Integration in Comoving Coordinates is supported to carry out Computational Cosmology simulations

We follow the unit system of [enzo](http://enzo-project.org). 
The implementation is found in the `src/cosmology.jl` file. 
To specify the cosmology and calculate redshifts and times we use the [Cosmology.jl package](https://github.com/JuliaAstro/Cosmology.jl). In our `.conf` files we can specify e.g. a flat ΛCDM model with `Cosmology = cosmology(OmegaM=0.27)` or also specify the the hubble constant in `100km/s/Mpc` with `Cosmology = cosmology(OmegaM=0.27,h=0.69)`. 

If we specified a Cosmology the code will use `InitialRedshift` and `FinalRedshift` as specified in the `.conf` file to set the start and end time overriding whatever we may have put as `StartTime` and `StopTime`. 


