module HydroEOS

export EOS, EOSIsothermal

export ε_from_P, P_from_ε, cs_from_ε, cs_from_P

type EOS
    γ::Real
end

ε_from_P(P, ρ, eos::EOS) = P/(ρ * (eos.γ-1))
P_from_ε(ε, ρ, eos::EOS) = (ρ * (eos.γ-1)) * ε
cs_from_ε(ε, eos::EOS) = sqrt(eos.γ * (eos.γ-1) * ε)
cs_from_P(P, ρ, eos::EOS) = sqrt(eos.γ * P/ρ)

type EOSIsothermal
    cs::Real
end

ε_from_P(P, ρ, eos::EOSIsothermal) = eos.cs^2
P_from_ε(ε, ρ, eos::EOSIsothermal) = ρ * eos.cs^2
cs_from_ε(ε, eos::EOSIsothermal) = eos.cs
cs_from_P(P, ρ, eos::EOSIsothermal) = eos.cs

end
