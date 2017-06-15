
using DecFP: Dec128


macro define(name, definition)
    quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


macro reexport(ex)
    isa(ex, Expr) && (ex.head == :module ||
                      ex.head == :using ||
                      ex.head == :importall ||
                      (ex.head == :toplevel &&
                       all(e->isa(e, Expr) && (e.head == :using || e.head == :importall), ex.args))) ||
        error("@reexport: syntax error")

    if ex.head == :module
        modules = Any[ex.args[2]]
        ex = Expr(:toplevel, ex, Expr(:using, :., ex.args[2]))
    elseif ex.head == :using || ex.head == :importall
        modules = Any[ex.args[end]]
    else
        modules = Any[e.args[end] for e in ex.args]
    end

    esc(Expr(:toplevel, ex,
             [:(eval(Expr(:export, names($mod)...))) for mod in modules]...))
end


function _dec128(x::Expr)
    y = x
    y.args .= _dec128.(y.args)
    return  y
end

_dec128(x::Int) = Dec128(x)

_dec128(x) = x

macro dec128(x)
    return esc(_dec128(x))
end
