"""
`InitialGuessIODE`: Initial guess for partitioned ordinary differential equations

### Fields

* `int`: interpolation structure
* `v`:   vector field for q
* `f`:   vector field for p
* `Δt`:  time step
* `s`:   number of extrapolation stages (for initialisation)
"""
mutable struct InitialGuessIODE{DT, TT, VT, FT, IT <: Interpolator}
    int::IT
    v::VT
    f::FT
    Δt::TT
    s::Int

    function InitialGuessIODE{DT,TT,VT,FT,IT}(interp, v, f, Δt) where {DT,TT,VT,FT,IT}
        new(interp, v, f, Δt, get_config(:ig_extrapolation_stages))
    end
end

function InitialGuessIODE{DT,D}(interp, v::VT, f::FT, Δt::TT) where {D, DT, TT, VT <: Function, FT <: Function}
    InitialGuessIODE{DT,TT,VT,FT,interp}(interp(zero(DT), one(DT), Δt, D), v, f, Δt)
end

function InitialGuess(interp, equ::IODE{DT,TT,ϑT,FT,GT,V̄T,F̄T}, Δt::TT) where {DT,TT,ϑT,FT,GT,V̄T,F̄T}
    InitialGuessIODE{DT,TT,V̄T,F̄T,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v̄, equ.f̄, Δt)
end

function InitialGuess(interp, equ::VODE{DT,TT,ϑT,FT,GT,V̄T,F̄T}, Δt::TT) where {DT,TT,ϑT,FT,GT,V̄T,F̄T}
    InitialGuessIODE{DT,TT,V̄T,F̄T,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v̄, equ.f̄, Δt)
end

function InitialGuess(interp, equ::IDAE{DT,TT,ϑT,FT,UT,GT,ϕT,V̄T,F̄T}, Δt::TT) where {DT,TT,ϑT,FT,UT,GT,ϕT,V̄T,F̄T}
    InitialGuessIODE{DT,TT,V̄T,F̄T,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v̄, equ.f̄, Δt)
end

function InitialGuess(interp, equ::VDAE{DT,TT,θT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T}, Δt::TT) where {DT,TT,θT,FT,GT,G̅T,ϕT,ψT,V̄T,F̄T}
    InitialGuessIODE{DT,TT,V̄T,F̄T,interp}(interp(zero(DT), one(DT), Δt, equ.d), equ.v̄, equ.f̄, Δt)
end


Base.:(==)(ig1::InitialGuessIODE{DT1,TT1}, ig2::InitialGuessIODE{DT2,TT2}) where {DT1,DT2,TT1,TT2}= (
                                DT1 == DT2
                             && TT1 == TT2
                             && ig1.int == ig2.int
                             && ig1.v   == ig2.v
                             && ig1.f   == ig2.f
                             && ig1.Δt  == ig2.Δt
                             && ig1.s   == ig2.s)


"Initialise initial guess, i.e., given t₀, t₁, q₁ compute q₀, v₀, v₁."
function initialize!(ig::InitialGuessIODE{DT,TT,VT,IT},
                t₁::TT,
                q₁::SolutionVector{DT},
                p₁::SolutionVector{DT},
                v₁::SolutionVector{DT},
                f₁::SolutionVector{DT},
                t₀::TT,
                q₀::SolutionVector{DT},
                p₀::SolutionVector{DT},
                v₀::SolutionVector{DT},
                f₀::SolutionVector{DT}) where {DT,TT,VT,IT}

    midpoint_extrapolation((t,q,p,v)->ig.v(t,q,v), ig.f, t₁, t₀, q₁, q₀, p₁, p₀, ig.s)

    ig.v(t₀, q₀, v₀)
    ig.v(t₁, q₁, v₁)
    ig.f(t₀, q₀, v₀, f₀)
    ig.f(t₁, q₁, v₁, f₁)
end


function update_vector_fields!(ig::InitialGuessIODE{DT,TT}, t₁::TT,
                               q₁::SolutionVector{DT},
                               p₁::SolutionVector{DT},
                               v₁::Vector{DT},
                               f₁::Vector{DT}) where {DT,TT}
    ig.v(t₁, q₁, v₁)
    ig.f(t₁, q₁, v₁, f₁)
end


function CommonFunctions.evaluate!(ig::InitialGuessIODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q)
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessIODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q, guess_v)
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessIODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁."
        guess_p .= p₁
    else
        evaluate!(ig.int, p₀, p₁, f₀, f₁, one(TT)+c_p, guess_p)
    end
end


function CommonFunctions.evaluate!(ig::InitialGuessIODE{DT,TT},
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
        evaluate!(ig.int, q₀, q₁, v₀, v₁, one(TT)+c_q, guess_q, guess_v)
    end

    if p₀ == p₁
        @warn "p₀ and p₁ in initial guess are identical! Setting p=p₁ and f=0."
        guess_p .= p₁
        guess_f .= 0
    else
        evaluate!(ig.int, p₀, p₁, f₀, f₁, one(TT)+c_p, guess_p, guess_f)
    end
end
