module TimestepperTheory

export vertical_modes, averaging_filter, barotropic_response, stage_matrices, evolution_matrix,
       exact_evolution, uncoupled_matrix, coupling_corrections, frozen_predictions, progressive_predictions, global_error,
       scheme_configuration, stability_polynomial, amplification_factor

export NoFilter, F1Filter, F2Filter, F3Filter,
       RungeKutta2, WickerSkamarock3, ShuOsher3, ModifiedRungeKutta4, AdamsBashforth2,
       ForwardBackward, RungeKutta3, ExactSolver, ImplicitSolver,
       FrozenForcing, FinalStageQuadraticForcing, ProgressiveForcing, SSPForcing,
       MidpointSampling, ContinuousSampling

using DocStringExtensions: TYPEDSIGNATURES
using GenericSchur: GenericSchur
using LinearAlgebra: I, Diagonal, eigen, eigvals, norm, opnorm, dot

include("algorithms.jl")
include("vertical_modes.jl")
include("averaging_filters.jl")
include("barotropic_response.jl")
include("stage_matrices.jl")
include("amplification_factors.jl")
include("diagnostics.jl")
include("progressive_predictions.jl")
include("rotation.jl")

export rotating_evolution_matrix, FrozenRotation, ProgressivePressure, ProgressiveRotation, SSPRotation, SSPSubstepRotation

end
