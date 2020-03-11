
abstract type AbstractIntegratorSPARK{dType, tType} <: IntegratorPRK{dType, tType} end
abstract type AbstractIntegratorHSPARK{dType, tType} <: AbstractIntegratorSPARK{dType, tType} end
abstract type AbstractIntegratorVSPARK{dType, tType} <: AbstractIntegratorSPARK{dType, tType} end

abstract type AbstractParametersSPARK{DT,TT} <: Parameters{DT,TT} end

@define ParametersSPARK begin
    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}
end

function update_params!(params::AbstractParametersSPARK, sol::AtomicSolutionPDAE)
    # set time for nonlinear solver and copy previous solution
    params.t  = sol.t
    params.q .= sol.q
    params.p .= sol.p
    params.λ .= sol.λ
end
