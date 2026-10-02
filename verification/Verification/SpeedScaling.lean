import Mathlib

/-!
Speed-scaling certificate. No axioms or sorry.
The analytic interface is explicit: positive periodic trajectories, their linear
sensitivity ODE, and differentiation of the reward integral. Existence and smooth
dependence of the nonlinear periodic ODE are proved in the companion note, not here.
-/
open Set MeasureTheory intervalIntegral
noncomputable section
namespace SpeedScaling

def gap (z y : ℝ) : ℝ :=
  (y-z)^2 * (1 + 2*(y+z) + y^2 + y*z + z^2) /
    (y*(1+z)^2*(1+y)^2)

def primitive (y : ℝ) : ℝ := Real.log y - Real.log (1+y) + 1/(1+y)

theorem gap_identity {z y : ℝ} (hz : 0 < z) (hy : 0 < y) :
    (y-z)/(1+z)^2 - z*(y-z)/(y*(1+y)^2) = gap z y := by
  unfold gap
  field_simp
  ring

theorem gap_nonneg {z y : ℝ} (hz : 0 < z) (hy : 0 < y) : 0 ≤ gap z y := by
  unfold gap
  positivity

theorem gap_pos {z y : ℝ} (hz : 0 < z) (hy : 0 < y) (hne : y ≠ z) :
    0 < gap z y := by
  unfold gap
  have hs : 0 < (y-z)^2 := sq_pos_of_ne_zero (sub_ne_zero.mpr hne)
  positivity

theorem primitive_derivative {y : ℝ} (hy : 0 < y) :
    HasDerivAt primitive (1/(y*(1+y)^2)) y := by
  have h1 : (1:ℝ)+y ≠ 0 := ne_of_gt (by linarith)
  have ha := (hasDerivAt_id y).const_add 1
  have h := ((Real.hasDerivAt_log (ne_of_gt hy)).sub (ha.log h1)).add
    ((hasDerivAt_const y (1:ℝ)).div ha h1)
  convert! h using 1
  simp only [id_eq]
  field_simp
  ring

/-- Algebraic consequence of differentiating k z' = S - α z². -/
theorem sensitivity_filter {k α z v zp vp : ℝ} (hk : k ≠ 0)
    (h : k*vp + 2*α*z*v = -zp) :
    zp + k*vp = (2*α/k)*z*(z-(z+k*v)) := by
  have he : zp + k*vp = -2*α*z*v := by linarith
  rw [he]
  field_simp
  ring

/-- Periodic total derivatives integrate to zero (FTC, not an assumed integral identity). -/
theorem periodic_primitive_integral
    (T k α : ℝ) (_hT : 0 < T) (hk : 0 < k) (hα : 0 < α)
    (z y : ℝ → ℝ) (hz : Continuous z) (hy : Continuous y)
    (hyp : ∀ x, 0 < y x) (hper : y T = y 0)
    (hode : ∀ x, HasDerivAt y ((2*α/k)*z x*(z x-y x)) x) :
    (∫ x in 0..T, z x*(y x-z x)/(y x*(1+y x)^2)) = 0 := by
  let f := fun x => z x*(y x-z x)/(y x*(1+y x)^2)
  have hf : Continuous f := by
    unfold f
    exact (hz.mul (hy.sub hz)).div (hy.mul ((continuous_const.add hy).pow 2))
      (fun x => ne_of_gt (by have hp := hyp x; positivity))
  have hd : ∀ x, HasDerivAt (fun t => primitive (y t)) ((-(2*α/k))*f x) x := by
    intro x
    convert! (primitive_derivative (hyp x)).comp x (hode x) using 1
    dsimp [f]
    ring
  have hi := integral_eq_sub_of_hasDerivAt (fun x _ => hd x)
    ((continuous_const.mul hf).intervalIntegrable 0 T)
  rw [hper, sub_self, intervalIntegral.integral_const_mul] at hi
  have hn : -(2*α/k) ≠ 0 := neg_ne_zero.mpr (ne_of_gt (by positivity))
  exact (mul_eq_zero.mp hi).resolve_left hn

/-- Exact nonnegative-integrand identity for the sensitivity of clarity. -/
theorem reward_identity
    (T k α : ℝ) (hT : 0 < T) (hk : 0 < k) (hα : 0 < α)
    (z y : ℝ → ℝ) (hz : Continuous z) (hy : Continuous y)
    (hzp : ∀ x, 0 < z x) (hyp : ∀ x, 0 < y x) (hper : y T = y 0)
    (hode : ∀ x, HasDerivAt y ((2*α/k)*z x*(z x-y x)) x) :
    (∫ x in 0..T, (y x-z x)/(1+z x)^2) = ∫ x in 0..T, gap (z x) (y x) := by
  have hf : Continuous (fun x => (y x-z x)/(1+z x)^2) :=
    (hy.sub hz).div ((continuous_const.add hz).pow 2)
      (fun x => ne_of_gt (by have hp := hzp x; positivity))
  have hg : Continuous (fun x => z x*(y x-z x)/(y x*(1+y x)^2)) :=
    (hz.mul (hy.sub hz)).div (hy.mul ((continuous_const.add hy).pow 2))
      (fun x => ne_of_gt (by have hp := hyp x; positivity))
  have he : (fun x => gap (z x) (y x)) =
      (fun x => (y x-z x)/(1+z x)^2 - z x*(y x-z x)/(y x*(1+y x)^2)) := by
    funext x
    exact (gap_identity (hzp x) (hyp x)).symm
  rw [he, integral_sub (hf.intervalIntegrable 0 T) (hg.intervalIntegrable 0 T),
    periodic_primitive_integral T k α hT hk hα z y hz hy hyp hper hode, sub_zero]

/-- Strictness requires a genuine mismatch at some phase, not merely S ≠ 0. -/
theorem reward_sensitivity_positive
    (T k α dR : ℝ) (hT : 0 < T) (hk : 0 < k) (hα : 0 < α)
    (z y : ℝ → ℝ) (hz : Continuous z) (hy : Continuous y)
    (hzp : ∀ x, 0 < z x) (hyp : ∀ x, 0 < y x) (hper : y T = y 0)
    (hode : ∀ x, HasDerivAt y ((2*α/k)*z x*(z x-y x)) x)
    (hne : ∃ x ∈ Icc 0 T, y x ≠ z x)
    (hdR : k*T*dR = ∫ x in 0..T, (y x-z x)/(1+z x)^2) : 0 < dR := by
  have hg : Continuous (fun x => gap (z x) (y x)) := by
    unfold gap
    fun_prop (disch := intro x; have hp := hyp x; have hq := hzp x; exact ne_of_gt (by positivity))
  have hp := integral_lt_integral_of_continuousOn_of_le_of_exists_lt hT
    continuous_const.continuousOn hg.continuousOn
    (fun x _ => gap_nonneg (hzp x) (hyp x))
    (by obtain ⟨x,hx,hn⟩ := hne; exact ⟨x,hx,gap_pos (hzp x) (hyp x) hn⟩)
  simp only [intervalIntegral.integral_zero] at hp
  rw [← reward_identity T k α hT hk hα z y hz hy hzp hyp hper hode, ← hdR] at hp
  exact pos_of_mul_pos_right hp (by positivity : 0 ≤ k*T)

/-- A derivative theorem yields finite, not just infinitesimal, comparisons. -/
theorem finite_scaling {R : ℝ → ℝ}
    (hc : ContinuousOn R (Ioi 0)) (hd : ∀ k ∈ Ioi (0:ℝ), 0 < deriv R k) :
    StrictMonoOn R (Ioi 0) := by
  apply strictMonoOn_of_deriv_pos (convex_Ioi 0) hc
  simpa only [interior_Ioi] using hd

/-- IVT plus strict monotonicity gives the unique energy scaling.
The companion note proves these hypotheses for the stated propulsion model. -/
theorem energy_scale_exists_unique (D : ℝ → ℝ)
    (hc : ContinuousOn D (Ici 1)) (hm : StrictMonoOn D (Ici 1))
    (hbase : D 1 = 0) (hunbounded : ∀ E : ℝ, ∃ K ≥ 1, E ≤ D K)
    (Δ : ℝ) (hΔ : 0 < Δ) : ∃! k : ℝ, 1 < k ∧ D k = Δ := by
  obtain ⟨K, hK, hDK⟩ := hunbounded Δ
  have hcK : ContinuousOn D (Icc 1 K) := hc.mono (fun _ hx => hx.1)
  have hmem : Δ ∈ Icc (D 1) (D K) := ⟨by rw [hbase]; exact hΔ.le, hDK⟩
  obtain ⟨k, hk, he⟩ := intermediate_value_Icc hK hcK hmem
  have hkgt : 1 < k := by
    have hn : k ≠ 1 := by intro heq; rw [heq, hbase] at he; linarith
    exact lt_of_le_of_ne hk.1 (Ne.symm hn)
  refine ⟨k, ⟨hkgt, he⟩, ?_⟩
  intro l hl
  exact hm.injOn hl.1.le hk.1 (hl.2.trans he.symm)

/-- Endpoint accounting: extra starting energy equals extra net lap expenditure. -/
theorem battery_endpoint {b₀ bf D₁ D₂ Δ : ℝ}
    (hbase : b₀-D₁=bf) (hscale : D₂=D₁+Δ) : b₀+Δ-D₂=bf := by linarith

/-- Optimizer attainment is essential in the strict value comparison used here. -/
theorem strict_value_transfer {V₁ V₂ R₁ R₂ : ℝ}
    (hopt : R₁=V₁) (himprove : R₁<R₂) (hfeasible : R₂≤V₂) : V₁<V₂ := by linarith

/-- Constant nonzero sensing is an explicit obstruction to unconditional strictness. -/
theorem constant_sensing_counterexample (k : ℝ) :
    (0:ℝ) = (1/k)*((1-(1/2:ℝ))^2-(1/2:ℝ)^2) := by ring

#print axioms reward_sensitivity_positive
#print axioms finite_scaling
#print axioms energy_scale_exists_unique
end SpeedScaling
