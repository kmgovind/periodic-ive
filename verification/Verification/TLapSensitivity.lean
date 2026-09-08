import Mathlib

open MeasureTheory intervalIntegral Filter Topology

noncomputable section

variable (u δu : ℝ → ℝ) (L T₀ : ℝ)

def A (T : ℝ) : ℝ := ∫ t in (0 : ℝ)..T, u t
def B (T : ℝ) : ℝ := ∫ t in (0 : ℝ)..T, δu t

theorem lap_time_sensitivity
    (T : ℝ → ℝ)
    (hT₀ : T 0 = T₀)
    (hconstraint : ∀ᶠ ε in 𝓝 (0 : ℝ), A u (T ε) + ε * B δu (T ε) = L)
    (hu_cont : Continuous u) (hδu_cont : Continuous δu)
    (T' : ℝ) (htDeriv : HasDerivAt T T' 0)
    (hu_ne : u T₀ ≠ 0) :
    T' = -(B δu T₀) / (u T₀) := by
  have huInt : IntervalIntegrable u volume 0 T₀ := hu_cont.intervalIntegrable 0 T₀
  have hδuInt : IntervalIntegrable δu volume 0 T₀ := hδu_cont.intervalIntegrable 0 T₀
  -- FTC facts, stated plainly at the constant T₀ -- no rewriting needed here at all
  have hA' : HasDerivAt (A u) (u T₀) T₀ :=
    intervalIntegral.integral_hasDerivAt_right huInt
      (hu_cont.stronglyMeasurableAtFilter volume (𝓝 T₀))
      hu_cont.continuousAt
  have hB' : HasDerivAt (B δu) (δu T₀) T₀ :=
    intervalIntegral.integral_hasDerivAt_right hδuInt
      (hδu_cont.stronglyMeasurableAtFilter volume (𝓝 T₀))
      hδu_cont.continuousAt
  -- Compose via the `_of_eq` chain-rule variant, which takes the point
  -- mismatch (T 0 vs T₀) as an explicit *equality hypothesis* rather than
  -- requiring syntactic/defeq agreement -- this sidesteps both the `▸`
  -- rewrite-direction issue and the instance-path diamond noise from before.
  -- Signature pattern (from HasStrictDerivAt.comp_of_eq / .scomp_of_eq):
  --   HasDerivAt.comp_of_eq (x) (hg : HasDerivAt g g' y) (hf : HasDerivAt f f' x)
  --                          (hy : y = f x) : HasDerivAt (g ∘ f) (g' * f') x
  have hAT : HasDerivAt (fun ε => A u (T ε)) (u T₀ * T') 0 := by
    have h := HasDerivAt.comp_of_eq 0 hA' htDeriv hT₀.symm
    -- h : HasDerivAt (A u ∘ T) (u T₀ * T') 0 -- unfold `∘` explicitly via `show`
    show HasDerivAt (fun ε => (A u ∘ T) ε) (u T₀ * T') 0
    exact h
  have hBT : HasDerivAt (fun ε => B δu (T ε)) (δu T₀ * T') 0 := by
    have h := HasDerivAt.comp_of_eq 0 hB' htDeriv hT₀.symm
    show HasDerivAt (fun ε => (B δu ∘ T) ε) (δu T₀ * T') 0
    exact h
  -- product rule for ε ↦ ε * B δu (T ε): build it as `(fun ε => ε) * (fun ε => B δu (T ε))`
  -- explicitly rather than relying on `id * f` normalizing to the fun-binder shape.
  have hεBT : HasDerivAt (fun ε => ε * B δu (T ε)) (B δu T₀) 0 := by
    have hid : HasDerivAt (fun ε : ℝ => ε) 1 0 := hasDerivAt_id 0
    have h := hid.mul hBT
    -- h : HasDerivAt ((fun ε => ε) * fun ε => B δu (T ε))
    --       (1 * B δu (T 0) + 0 * (δu T₀ * T')) 0
    -- two gaps vs. the goal: (a) T 0 vs T₀ in the value, (b) Pi.mul-form vs
    -- fun-binder-form in the domain -- close (a) with hT₀, bridge (b) with `show`
    rw [hT₀] at h
    have hval : (1:ℝ) * B δu T₀ + 0 * (δu T₀ * T') = B δu T₀ := by ring
    rw [hval] at h
    show HasDerivAt (fun ε => ε * B δu (T ε)) (B δu T₀) 0
    exact h
  have hSum : HasDerivAt (fun ε => A u (T ε) + ε * B δu (T ε))
      (u T₀ * T' + B δu T₀) 0 := hAT.add hεBT
  have hEq : (fun ε => A u (T ε) + ε * B δu (T ε)) =ᶠ[𝓝 (0 : ℝ)] (fun _ => L) :=
    hconstraint
  have hZero : HasDerivAt (fun ε => A u (T ε) + ε * B δu (T ε)) 0 0 :=
    (hEq.hasDerivAt_iff).mpr (hasDerivAt_const 0 L)
  have hEqDeriv : u T₀ * T' + B δu T₀ = 0 := hSum.unique hZero
  field_simp
  linarith [hEqDeriv]

end
