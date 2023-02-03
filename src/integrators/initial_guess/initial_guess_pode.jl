"""
`InitialGuessPODE`: Initial guess for partitioned ordinary differential equations

### Fields

* `int`: interpolation structure
* `v`:   vector field for q
* `f`:   vector field for p
* `Δt`:  time step
* `s`:   number of extrapolation stages (for initialisation)
"""
struct InitialGuessPODE{TT, VT, FT, IT <: Extrapolation}
    int::IT
    v::VT
    f::FT
    Δt::TT
    s::Int

    function InitialGuessPODE(int::IT, v::VT, f::FT, Δt::TT) where {TT <: Real, VT <: Function, FT <: Function, IT <: Extrapolation}
        new{TT,VT,FT,IT}(int, v, f, Δt, get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessPODE(interp::Type{<:Extrapolation}, v::Function, f::Function, Δt)
    int = interp()
    InitialGuessPODE(int, v, f, Δt)
end

function InitialGuessPODE(interp::Type{<:Extrapolation}, problem::Union{PODEProblem,PDAEProblem,HODEProblem,HDAEProblem}, Δt = tstep(problem))
    InitialGuessPODE(interp, _get_v̄(equation(problem), parameters(problem)), _get_f̄(equation(problem), parameters(problem)), Δt)
end

InitialGuess(interp, problem::Union{PODEProblem,PDAEProblem,HODEProblem,HDAEProblem}, Δt = tstep(problem)) = InitialGuessPODE(interp, problem, Δt)
InitialGuess(problem::Union{PODEProblem,PDAEProblem,HODEProblem,HDAEProblem}, Δt = tstep(problem)) = InitialGuessPODE(get_config(:ig_extrapolation), problem, Δt)


Base.:(==)(ig1::InitialGuessPODE{TT1}, ig2::InitialGuessPODE{TT2}) where {TT1,TT2}= (
                                TT1 == TT2
                             && ig1.int == ig2.int
                             && ig1.v   == ig2.v
                             && ig1.f   == ig2.f
                             && ig1.Δt  == ig2.Δt
                             && ig1.s   == ig2.s)


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessPODE{TT},
                t₁::TT,
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                t₀::TT,
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT}) where {DT,TT}

    extrapolate_pode!(t₀, q₀, p₀, t₁, q₁, p₁, ig.v, ig.f, MidpointExtrapolation(ig.s))

    ig.v(v₀, t₀, q₀, p₀)
    ig.v(v₁, t₁, q₁, p₁)
    ig.f(f₀, t₀, q₀, v₀)
    ig.f(f₁, t₁, q₁, v₁)
end


function Solutions.update_vector_fields!(ig::InitialGuessPODE{TT}, t₁::TT,
                               q₁::SolutionVector{DT},
                               p₁::SolutionVector{DT},
                               v₁::Vector{DT},
                               f₁::Vector{DT}) where {DT,TT}
    ig.v(v₁, t₁, q₁, p₁)
    ig.f(f₁, t₁, q₁, v₁)
end


function GeometricBase.evaluate!(ig::InitialGuessPODE{TT},
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                guess_q::SolutionVector{DT},
                c_q::TT) where {DT,TT}

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁."
        guess_q .= q₁
    else
        extrapolate!((one(TT)+c_q)*ig.Δt, guess_q, zero(TT), q₀, v₀, ig.Δt, q₁, v₁, ig.int)
    end
end


function GeometricBase.evaluate!(ig::InitialGuessPODE{TT},
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                guess_q::SolutionVector{DT},
                guess_v::SolutionVector{DT},
                c_q::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_v)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0."
        guess_q .= q₁
        guess_v .= 0
    else
        extrapolate!((one(TT)+c_q)*ig.Δt, guess_q, guess_v, zero(TT), q₀, v₀, ig.Δt, q₁, v₁, ig.int)
    end
end


function GeometricBase.evaluate!(ig::InitialGuessPODE{TT},
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                guess_q::SolutionVector{DT},
                guess_p::SolutionVector{DT},
                c_q::TT, c_p::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁."
        guess_q .= q₁
    else
        extrapolate!((one(TT)+c_q)*ig.Δt, guess_q, zero(TT), q₀, v₀, ig.Δt, q₁, v₁, ig.int)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁."
        guess_p .= p₁
    else
        extrapolate!((one(TT)+c_p)*ig.Δt, guess_p, zero(TT), p₀, f₀, ig.Δt, p₁, f₁, ig.int)
    end
end


function GeometricBase.evaluate!(ig::InitialGuessPODE{TT},
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT},
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                guess_q::SolutionVector{DT},
                guess_p::SolutionVector{DT},
                guess_v::SolutionVector{DT},
                guess_f::SolutionVector{DT},
                c_q::TT, c_p::TT) where {DT,TT}
    @assert length(guess_q) == length(guess_v)
    @assert length(guess_p) == length(guess_f)
    @assert length(guess_q) == length(guess_p)

    if q₀ == q₁
        @warn "q₀ and q₁ in initial guess are identical! Setting q=q₁ and v=0."
        guess_q .= q₁
        guess_v .= 0
    else
        extrapolate!((one(TT)+c_q)*ig.Δt, guess_q, guess_v, zero(TT), q₀, v₀, ig.Δt, q₁, v₁, ig.int)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁ and f=0."
        guess_p .= p₁
        guess_f .= 0
    else
        extrapolate!((one(TT)+c_p)*ig.Δt, guess_p, guess_f, zero(TT), p₀, f₀, ig.Δt, p₁, f₁, ig.int)
    end
end