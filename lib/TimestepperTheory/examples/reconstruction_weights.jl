# # Progressive forcing weights for MRK4
#
# For $a=0.23$, the stage durations are $\gamma_m\Delta t$, with
# ```math
# (\gamma_1,\gamma_2,\gamma_3,\gamma_4)=(0.23,1/3,1/2,1).
# ```
# Stage $m$ has forcing samples $\overline{\mathbf G}^{\,j}$ at $c_j\Delta t$, where
# $(c_0,c_1,c_2,c_3)=(0,0.23,1/3,1/2)$ and $j=0,\ldots,m-1$.
# The implementation uses $s=t'/(\gamma_m\Delta t)$. Thus the sample locations are $s_j=c_j/\gamma_m$.
# We write
# ```math
# \widehat{\mathbf G}^{\,m}(s)=\sum_{p=0}^{\min(m-1,2)}\mathbf C_p^m s^p,
# \qquad V_{j,p}=s_j^p,\qquad
# \mathbf C_p^m=\sum_{j=0}^{m-1}C^m_{p,j}\overline{\mathbf G}^{\,j}.
# ```
# The matrix $C^m$ solves $VC^m=I$ for stages 1–3. At stage 4, the same linear solve gives
# the least-squares quadratic through the four samples. Its rows are the constant, linear, and quadratic weights.

using LinearAlgebra
using TimestepperTheory: ModifiedRungeKutta4, stage_fractions, reconstruction_weights

a = 0.23
γ = stage_fractions(ModifiedRungeKutta4(; a), Float64)

for m in eachindex(γ)
    C = reconstruction_weights(m, γ)
    println("Stage $m: rows p = 0, 1, 2; columns j = 0, …, $(m - 1)")
    display(round.(C; digits=8))
end

# At the final stage $s=t'/\Delta t$, so the matrix being fitted is explicitly
# ```math
# V=\begin{bmatrix}
# 1&0&0\\1&0.23&0.23^2\\1&1/3&1/9\\1&1/2&1/4
# \end{bmatrix},\qquad C^4=\underset{C}{\operatorname{argmin}}\,\|VC-I\|_F^2.
# ```
# Compute these weights directly and compare with the routine used by the evolution matrices:

c = [0.0, a, 1 / 3, 1 / 2]
V = [cⱼ^p for cⱼ in c, p in 0:2]
C = V \ Matrix{Float64}(I, 4, 4)
@assert C ≈ reconstruction_weights(4, γ)

# A quadratic must be recovered exactly. For example, samples of $2-3s+4s^2$ give coefficients $(2,-3,4)$:

G = @. 2 - 3c + 4c^2
@assert C * G ≈ [2, -3, 4]
display(C * G)

# To obtain the manuscript coefficients multiplying physical time $t'^p$, divide row $p$ by the stage duration
# raised to power $p$:
# ```math
# \mathbf G_p^m=\frac{\sum_j C^m_{p,j}\overline{\mathbf G}^{\,j}}{(\gamma_m\Delta t)^p}.
# ```
# The weights depend only on the stage fractions, so they are computed once and applied to each set of stage values.
