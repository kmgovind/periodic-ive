import Mathlib

open Real Filter Topology

noncomputable section

theorem battery_trajectory_preservation
    (b₁ b₂ f : ℝ → ℝ)
    (hb₁ : ∀ t, HasDerivAt b₁ (f t) t)
    (hb₂ : ∀ t, HasDerivAt b₂ (f t) t) :
    ∀ t, b₁ t - b₂ t = b₁ 0 - b₂ 0 := by
  have hdiff : ∀ t, HasDerivAt (fun t => b₁ t - b₂ t) 0 t := by
    intro t
    have h1 := hb₁ t
    have h2 := hb₂ t
    have h_sub := h1.sub h2
    rwa [Pi.sub_def, sub_self] at h_sub
  have hconst : ∀ t, (fun t => b₁ t - b₂ t) t = (fun t => b₁ t - b₂ t) 0 :=
    fun t => is_const_of_deriv_eq_zero
      (fun x => (hdiff x).differentiableAt)
      (fun x => (hdiff x).deriv) t 0
  intro t
  exact hconst t

theorem battery_trajectory_preservation_strict
    (b₁ b₂ f : ℝ → ℝ)
    (hb₁ : ∀ t, HasDerivAt b₁ (f t) t)
    (hb₂ : ∀ t, HasDerivAt b₂ (f t) t)
    (h0 : b₂ 0 < b₁ 0) :
    ∀ t, b₂ t < b₁ t := by
  intro t
  have h := battery_trajectory_preservation b₁ b₂ f hb₁ hb₂ t
  linarith [h, h0]
