####
#### Outer time discretizations
####

abstract type AbstractRungeKuttaScheme end

struct RungeKutta2            <: AbstractRungeKuttaScheme end
struct WickerSkamarock3       <: AbstractRungeKuttaScheme end
struct ShuOsher3              <: AbstractRungeKuttaScheme end
struct ModifiedRungeKutta4{T} <: AbstractRungeKuttaScheme
    a :: T
end

ModifiedRungeKutta4(; a=23 // 100) = ModifiedRungeKutta4(a)

struct AdamsBashforth2 end

####
#### Barotropic solvers
####

abstract type AbstractBarotropicSolver end

struct ForwardBackward <: AbstractBarotropicSolver end
struct RungeKutta3     <: AbstractBarotropicSolver end
struct ExactSolver     <: AbstractBarotropicSolver end
struct ImplicitSolver  <: AbstractBarotropicSolver end

####
#### Averaging filters for the Split--Explicit Barotropic solver
####

abstract type AbstractAveragingFilter end

struct NoFilter <: AbstractAveragingFilter end
struct F1Filter <: AbstractAveragingFilter end
struct F2Filter <: AbstractAveragingFilter end
struct F3Filter <: AbstractAveragingFilter end

####
#### Sampling of the slow forcing reconstructions
####

abstract type AbstractForcingSampling end

struct MidpointSampling   <: AbstractForcingSampling end
struct ContinuousSampling <: AbstractForcingSampling end

####
#### Slow forcing reconstructions
####

abstract type AbstractForcingReconstruction end

struct FrozenForcing <: AbstractForcingReconstruction end
struct SSPForcing    <: AbstractForcingReconstruction end

struct FinalStageQuadraticForcing{S <: AbstractForcingSampling} <: AbstractForcingReconstruction
    sampling :: S
end

struct ProgressiveForcing{S <: AbstractForcingSampling} <: AbstractForcingReconstruction
    sampling :: S
end

FinalStageQuadraticForcing(; sampling=MidpointSampling()) = FinalStageQuadraticForcing(sampling)
ProgressiveForcing(; sampling=MidpointSampling()) = ProgressiveForcing(sampling)

forcing_sampling(::Union{FrozenForcing, SSPForcing}) = MidpointSampling()
forcing_sampling(forcing::Union{FinalStageQuadraticForcing, ProgressiveForcing}) = forcing.sampling
