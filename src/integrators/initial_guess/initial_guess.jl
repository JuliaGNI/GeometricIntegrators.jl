
abstract type InitialGuess end

struct NoInitialGuess <: InitialGuess end

const OptionalInitialGuess = Union{NoInitialGuess, InitialGuess}
