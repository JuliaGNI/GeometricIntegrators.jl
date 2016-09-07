
macro define(name, definition)
    return quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end
