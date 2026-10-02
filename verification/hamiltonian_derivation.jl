using Symbolics, Plots, LaTeXStrings, Latexify, SymbolicNumericIntegration

# ==============================================================================
# 1. Variable and Parameter Definitions
# ==============================================================================
N_nodes = 2
@variables t

# Define states, costates, and control
@variables s(t) b(t) u(t)
@variables λ_s(t) λ_b(t)

# Define array states (q) and costates (lambda_q)
@variables q[1:N_nodes] λ_q[1:N_nodes]

# Define parameters
@variables P_in_avg α t_f
@variables γ[1:N_nodes]

# Define functions for P_out and S
@variables P_out(..) S(..)[1:N_nodes]

# ==============================================================================
# 2. Objective Function & System Dynamics
# ==============================================================================
println("Objective function:")
g_x = (1 / t_f) * sum(γ[i] * q[i] for i in 1:N_nodes)
println(latexify(g_x))

ṡ = u
ḃ = P_in_avg - P_out(u)
q̇ = [S(s)[i] * (1 - q[i])^2 - α * q[i]^2 for i in 1:N_nodes]

println("q̇ = ")
println(latexify(q̇))
display(latexify(q̇))

# ==============================================================================
# 3. Hamiltonian Construction
# ==============================================================================
inner_product = λ_s * ṡ + λ_b * ḃ + sum(λ_q[i] * q̇[i] for i in 1:N_nodes)
H = g_x + inner_product

println("=== Hamiltonian ===")
println("H = ", H)
println()
println("H = ", latexify(H))
display(latexify(H))

# ==============================================================================
# 4. First-Order Necessary Conditions (Costates & Control Optimality)
# ==============================================================================
λ_ṡ = -Symbolics.derivative(H, s)
λ_ḃ = -Symbolics.derivative(H, b)
λ_q̇ = [-Symbolics.derivative(H, q[i]) for i in 1:N_nodes]

# Control Optimality: dH/du = 0
dH_du = Symbolics.derivative(H, u)

println("=== Costate Dynamics ===")
println("d(λ_s)/dt = ", λ_ṡ)
display(latexify(λ_ṡ))
println("d(λ_s)/dt = ", latexify(λ_ṡ))

println("d(λ_b)/dt = ", λ_ḃ)   # Outputs 0 -> lambda_b is constant
display(latexify(λ_ḃ))
println("d(λ_b)/dt = ", latexify(λ_ḃ))

for i in 1:N_nodes
    println("d(λ_q_$i)/dt = ", λ_q̇[i])
    display(latexify(λ_q̇[i]))
    println("d(λ_q_$i)/dt = ", latexify(λ_q̇[i]))
end
println()

println("=== Control Optimality Condition (dH/du = 0) ===")
println("dH/du = ", dH_du)
display(latexify(dH_du))
println("dH/du = ", latexify(dH_du))
println()


# ==============================================================================
# 5. Symbolic Solution of λ_q Costate ODEs via Integrating Factors
# ==============================================================================
println("\n=== 5. Symbolic Solution of λ_q Costate ODEs ===")

@variables τ ξ
@variables λ_q_0[1:N_nodes]

# Declare uninterpreted symbolic integral operator (version-agnostic)
@variables Integral(..)
symbolic_int(expr, var, lower, upper) = Integral(expr, var, lower, upper)

# Extract linear ODE parameters: dλ_q_i/dt = A_i(t)*λ_q_i(t) + B_i
# A_i(t) = 2 * (S_i(s)*(1 - q_i) + α*q_i)
A = [2 * (S(s)[i] * (1 - q[i]) + α * q[i]) for i in 1:N_nodes]
B = [-γ[i] / t_f for i in 1:N_nodes]

λ_q_sol = Vector{Any}(undef, N_nodes)
λ_q_0_periodic = Vector{Any}(undef, N_nodes)
Φ_q = Vector{Any}(undef, N_nodes)

for i in 1:N_nodes
    # 1. Integrating factor: μ_i(t) = exp(- ∫_0^t A_i(τ) dτ)
    int_A_t = symbolic_int(A[i], t, 0, t)
    μ_t = exp(-int_A_t)
    
    # 2. General Initial-Value Solution: λ_q_i(t) = (λ_q_i(0) + ∫_0^t exp(-∫_0^τ A_i(ξ)dξ)*B_i dτ) / μ_i(t)
    int_A_τ = symbolic_int(A[i], t, 0, τ)
    forcing_integrand = exp(-int_A_τ) * B[i]
    int_forcing_t = symbolic_int(forcing_integrand, τ, 0, t)
    
    λ_q_sol[i] = (λ_q_0[i] + int_forcing_t) / μ_t
    
    # 3. Periodic Boundary Condition: λ_q_i(0) = λ_q_i(t_f)
    int_A_tf = symbolic_int(A[i], t, 0, t_f)
    μ_tf = exp(-int_A_tf)
    int_forcing_tf_norm = symbolic_int(exp(-int_A_τ) * (-1 / t_f), τ, 0, t_f)
    
    # Initial state λ_q_i(0) under periodicity
    λ_q_0_periodic[i] = γ[i] * (int_forcing_tf_norm / (μ_tf - 1))
    
    # 4. Factorized Trajectory: λ_q_i(t) = γ_i * Φ_i(t)
    int_forcing_t_norm = symbolic_int(exp(-int_A_τ) * (-1 / t_f), τ, 0, t)
    Φ_q[i] = ((int_forcing_tf_norm / (μ_tf - 1)) + int_forcing_t_norm) / μ_t
    
    # Output LaTeX expressions
    println("\n--- Node $i Costate Solution ---")
    println("General IVP Solution λ_q_$i(t):")
    display(latexify(λ_q_sol[i]))
    
    println("Periodic Initial Costate λ_q_$i(0):")
    display(latexify(λ_q_0_periodic[i]))
    
    println("Proportional Factorization λ_q_$i(t) = γ_$i * Φ_$i(t):")
    display(latexify(γ[i] * Φ_q[i]))
end