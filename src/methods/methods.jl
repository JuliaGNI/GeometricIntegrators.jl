
abstract type Method end

abstract type ODEMethod <: Method end
abstract type PODEMethod <: Method end
abstract type HODEMethod <: PODEMethod end
abstract type IODEMethod <: Method end
abstract type LODEMethod <: IODEMethod end
abstract type SODEMethod <: Method end

abstract type DAEMethod <: Method end
abstract type PDAEMethod <: Method end
abstract type HDAEMethod <: PDAEMethod end
abstract type IDAEMethod <: Method end
abstract type LDAEMethod <: IDAEMethod end
