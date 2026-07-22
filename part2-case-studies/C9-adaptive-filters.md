---
title: "C9 · Adaptive filters"
parent: Case studies
nav_order: 9
---

# C9 — Cerebellum-like circuits → adaptive filters and negative images
{: .no_toc }

> **Circuit.** The mormyrid electrosensory lateral line lobe (ELL) — a cerebellum-like structure in which molecular-layer parallel fibres, carrying corollary discharge and proprioceptive signals, synapse on the apical dendrites of principal (medium ganglion) cells that also receive raw electroreceptor input.
> **Primitive extracted.** Least-mean-squares adaptive filtering: the circuit learns a *negative image* of its own predictable input and subtracts it, leaving an innovation. The anti-Hebbian plasticity rule **is** stochastic gradient descent on the squared postsynaptic response.
> **Status.** As settled as anything in systems neuroscience. The learning rule was measured directly, the algorithm predicts the rule's *sign* and *timing asymmetry*, and the predictions held.
> **Structures thread.** Projection and least squares; orthogonality principles and innovations; stochastic approximation. The basis-span argument of §2.2 is the same question asked in `../structures/S-11-expanders-and-optimal-degree.md` (what a random expansion can represent). Same threads as `./C7-dopamine-and-td-learning.md`, arrived at from the opposite direction; index in `../structures/README.md`.

---

## Contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## 1. The phenomenon

A mormyrid electric fish emits a brief electric organ discharge (EOD) several times per second and senses the resulting field with skin electroreceptors. This is a signal-detection nightmare of its own making: the reafferent field produced by its own discharge is enormous compared with the perturbation caused by a nearby water flea, and it recurs on every command. The animal is shouting into a microphone and trying to hear a whisper.

Curtis Bell showed in 1981 that ELL principal cells solve this by *learning* the reafference. Pair the electric organ corollary discharge (EOCD) with an arbitrary, experimenter-imposed sensory pattern for a few minutes and the cell's response to that pattern shrinks to nothing. Then remove the imposed pattern and deliver the command alone: the cell produces a clean **negative image** — a response with the same time course as the cancelled input and the opposite sign, decaying over minutes. Nothing about the imposed pattern was natural. The circuit had built a bespoke predictor.

Bell, Han, Sugawara & Grant (1997) then went after the synapse. Pairing a parallel-fibre EPSP with a postsynaptic dendritic broad spike produces **depression** if the EPSP precedes the spike by tens of milliseconds, and no depression (or relative potentiation) if the order is reversed or the inputs are unpaired. A non-associative potentiation slowly restores baseline. This is spike-timing-dependent plasticity with the sign flipped relative to hippocampal STDP — hence *anti-Hebbian*.

The temptation is to file this as "an interesting inverted plasticity rule". Resist it. The sign is not a curiosity. It is a theorem.

## 2. The conversion

### 2.1 Set-up

Let the command cycle define a time variable $t$ (time since the last EOCD). Parallel fibres provide a $K$-dimensional temporal basis $p(t)=(p_1(t),\dots,p_K(t))^\top$, locked to the command — granule cells in the eminentia granularis, which Kennedy et al. (2014) showed tile the post-command interval with delayed, diverse responses. The principal cell receives raw electrosensory drive $R(t)$ and the weighted parallel-fibre input, so its output (relative to baseline) is

$$E(t) = R(t) + w^\top p(t).$$

Two conventions matter. First, $w$ is *signed*: parallel fibres are excitatory, but they also drive feed-forward inhibitory interneurons, so the net effect of a fibre on the principal cell can be either sign, and depression of the excitatory limb below its inhibitory counterpart yields a net negative weight. Second — and this is the whole trick — **the cell's own output is the error signal.** There is no teacher. The target is zero.

Decompose the input into a command-locked (reafferent) part and everything else:

$$R(t) = R_{\text{re}}(t) + R_{\text{ex}}(t), \qquad \mathbb{E}[R_{\text{ex}}(t)\,p_j(t)] = 0 \;\;\forall j,$$

the orthogonality holding because prey are not correlated with the fish's motor commands. This one statistical assumption is doing all the work, and it is worth naming: *self-generated input is exactly the input predictable from a copy of the command; everything else is news.*

### 2.2 The objective and the Wiener solution

Take the objective to be the mean square output,

$$J(w) = \tfrac12\big\langle E(t)^2\big\rangle,$$

where $\langle\cdot\rangle$ averages over $t$ within the cycle and over trials. Then

$$\frac{\partial J}{\partial w_j} = \big\langle E(t)\,p_j(t)\big\rangle,\qquad \nabla_w J = C w + q,\quad C = \langle p p^\top\rangle,\; q = \langle p R\rangle .$$

Setting $\nabla_w J = 0$ gives the Wiener–Hopf solution

$$\boxed{\,w^* = -C^{-1}q\,}$$

and the **orthogonality principle** $\langle E^* p_j\rangle = 0$ for all $j$: at the optimum the residual is uncorrelated with every parallel fibre. Writing $\Pi$ for orthogonal projection onto $\operatorname{span}\{p_j\}$ under this inner product,

$$w^{*\top}p = -\Pi R, \qquad E^* = (I-\Pi)R = R_{\text{ex}} + (I-\Pi)R_{\text{re}} .$$

If the granule-cell basis spans the reafference, $(I-\Pi)R_{\text{re}}=0$ and the cell's output *is* the exafferent signal — the whisper, with the shout removed. The negative image $-\Pi R$ is not a stored copy of anything: it is a projection, the best approximation to the reafference available in the span of the temporal basis.

### 2.3 The learning rule is gradient descent, and it must be anti-Hebbian

Replace the ensemble average by its instantaneous sample (Robbins–Monro / Widrow–Hoff, 1960):

$$\Delta w_j = -\eta\,\frac{\partial J}{\partial w_j} \;\longrightarrow\; \boxed{\,\Delta w_j = -\eta\,E(t)\,p_j(t)\,}$$

Read that off the page slowly. Presynaptic activity $p_j$ times postsynaptic activity $E$, with a **minus** sign. Coincidence of pre and post *weakens* the synapse. The anti-Hebbian rule is not an exotic alternative to Hebbian plasticity; it is what LMS looks like when the postsynaptic cell's own output plays the role of the error. Any circuit whose job is to null a predictable input, and whose error is its own output, is forced into an anti-Hebbian rule. Bell's 1997 measurement is the *sign of a gradient*, measured with a patch pipette.

Two refinements bring the mathematics into contact with the biology.

**Timing.** Plasticity depends on the pre–post interval through a kernel, so the realistic rule is

$$\dot w_j(t) = -\eta \int_0^\infty\! K(\tau)\,p_j(t-\tau)\,\Phi[E(t)]\,d\tau \;+\; \beta\,(w_0 - w_j),$$

with $\Phi$ a (rectifying, broad-spike-mediated) readout of the postsynaptic response and $\beta$ the non-associative recovery. Convolving with $K$ replaces the basis $p$ by a filtered basis $\tilde p = K*p$; the fixed point becomes $\langle E\,\tilde p_j\rangle = 0$, which coincides with the true optimum iff $\tilde p$ spans the same subspace as $p$ (Exercise 5). The finite width of the STDP window is therefore not a biological imperfection — it sets *which* temporal features the circuit is able to cancel. Roberts & Bell (2000) worked out these consequences in detail and showed that the measured window shape yields stable, correctly-signed image formation.

**Decay.** The non-associative term makes the fixed point a ridge solution,

$$w^\dagger = \Big(C + \tfrac{\beta}{\eta}I\Big)^{-1}\Big(\tfrac{\beta}{\eta}w_0 - q\Big),$$

which regularises the inverse when $C$ is ill-conditioned (many granule cells, correlated basis functions) and lets the circuit forget stale images. Continual learning, implemented as weight decay, discovered by evolution before it was rediscovered by machine learning.

### 2.4 Dynamics: what LMS predicts that "it learns to cancel" does not

Averaging the update, with $v = w - w^*$:

$$\langle v_{n+1}\rangle = (I-\eta C)\,\langle v_n\rangle .$$

Three consequences, all measurable:

1. **Stability requires $0<\eta<2/\lambda_{\max}(C)$.** Too-fast plasticity oscillates and diverges. Predicts a ceiling on learning rate set by the *variance* of parallel-fibre activity.
2. **Learning is exponential with mode-specific rates.** In the eigenbasis of $C$, error components decay as $(1-\eta\lambda_i)^n$, i.e. with time constants $\tau_i = 1/(\eta\lambda_i)$. Components of the reafference aligned with high-variance directions of the granule-cell code are cancelled *fast*; components aligned with low-variance directions crawl. This is a strong, non-obvious prediction: cancellation speed for an imposed pattern should be predictable from the pattern's overlap with the principal components of parallel-fibre activity. Bell's original observation that some imposed waveforms are cancelled far faster than others is exactly this eigenvalue spread.
3. **Steady state is not zero error.** Stochastic gradient descent with fixed step size rattles around the optimum. The standard LMS misadjustment result gives excess mean-square error $J_{\text{exc}} \approx \tfrac{\eta}{2}\operatorname{tr}(C)\,J_{\min}$, i.e. a fractional overshoot $M=\tfrac{\eta}{2}\operatorname{tr}(C)$. **The circuit cannot be both fast and precise.** Predicts a quantitative trade-off between how quickly a fish adapts to a new load on its electric organ and the residual reafferent noise floor it tolerates.

### 2.5 The cerebellar reading, and its complication

Fujita (1982) proposed exactly this architecture for the cerebellum: mossy fibre → granule cell expansion as a basis, Purkinje cell as the summing junction, climbing fibre as the error, parallel-fibre LTD as the descent step. Dean, Porrill, Ekerot & Jörntell (2010) assembled the modern evidence. The formal skeleton is identical to ELL. The difference is where the error comes from: in ELL the principal cell *is* the error unit (target = 0, self-supervised); in cerebellum the error arrives on a separate wire (the climbing fibre), which raises a problem the ELL never faces. The climbing fibre reports a *sensory* error, but gradient descent on motor performance requires the error in *motor* coordinates — you need to back-propagate through the plant. Porrill, Dean & Stone (2004) showed that a recurrent architecture, in which the adaptive filter sits in a feedback loop around the plant, converts the available sensory error into the right thing: the correct rule becomes **decorrelation** of the filter's output from the sensory error, which is again $\langle E p\rangle \to 0$, no Jacobian required. Worth internalising: whenever you cannot get the true gradient, look for an architecture that makes decorrelation sufficient.

### 2.6 Why this is one of the most reusable readings in neuroscience

"Circuit subtracts a learned prediction of its own input" is a template, and once you have it you see it everywhere: corollary-discharge inhibition of cricket auditory neurons during singing (Poulet & Hedwig, 2006); motor-cortex-driven suppression of self-generated sounds in mouse auditory cortex, which is learned and sound-specific (Schneider, Nelson & Mooney, 2014; Schneider et al., 2018); the mormyrid ELL's own multi-layer, continually-learning version (Muller, Zadina, Abbott & Sawtell, 2019); and, at the level of pure architecture, predictive coding. What makes the ELL the *canonical* instance is that here the learning rule was independently measured, so the conversion is not an interpretation laid over a phenomenon — it is a match between an equation and a pipette.

The link to `../part1-foundations/06-control-and-filtering.md` is direct. There, the Kalman filter's innovation $\nu_t = y_t - H\hat x_{t|t-1}$ is white and orthogonal to the past by construction, and the gain is computed from a known model. Here the same orthogonality is *achieved by learning* rather than derived from a model, and the price is misadjustment (§2.4) and a basis-limited hypothesis class (§2.2). LMS is the model-free, biologically implementable cousin of the Kalman gain; Exercise 6 makes the correspondence quantitative.

---

## 3. Model to build

Forty lines that reproduce Bell's central experiment, including the afterimage.

```python
import numpy as np
rng = np.random.default_rng(1)

T, K = 200, 30                                  # 200 1-ms bins per command cycle; 30 PFs
t   = np.arange(T)
ctr = np.linspace(0, T, K)
wid = 2.0 * T / K
P   = np.exp(-0.5 * ((t[None, :] - ctr[:, None]) / wid) ** 2)   # K x T temporal basis

reaff = 1.2*np.exp(-0.5*((t-40)/8.0)**2) - 0.6*np.exp(-0.5*((t-90)/25.0)**2)

w, eta, beta = np.zeros(K), 2e-2, 1e-5
resid, images = [], []

for trial in range(4000):
    probe = trial >= 3000                       # after learning: silence the reafference
    prey = np.zeros(T)
    if rng.random() < 0.3:                      # command-independent "exafferent" input
        t0 = rng.integers(20, T - 20)
        prey += 0.8 * np.exp(-0.5 * ((t - t0) / 6.0) ** 2)
    R = (0.0 if probe else 1.0) * reaff + prey + 0.05 * rng.standard_normal(T)

    E = R + w @ P                               # postsynaptic response == error signal
    if not probe:
        w += (-eta * (P @ E) + beta * (0.0 - w)) / T   # anti-Hebbian + decay to baseline
        resid.append(np.mean((E - prey) ** 2))
    else:
        images.append(E)

negative_image = np.mean(images, axis=0)
print("residual power, first/last 100 learning trials:",
      np.mean(resid[:100]).round(4), np.mean(resid[-100:]).round(4))
print("corr(negative image, -reafference):",
      np.corrcoef(negative_image, -reaff)[0, 1].round(3))
```

**What to look for.**
- Residual power falls roughly an order of magnitude (≈0.085 → ≈0.009 with these settings), on a visibly multi-exponential trajectory — plot $\log$ residual and look for the *distinct* time constants of §2.4.
- The probe response is an inverted copy of the reafference: correlation with $-\texttt{reaff}$ of about 0.98. It is not 1.000, and the shortfall is not a bug — it is $(I-\Pi)R_{\text{re}}$, the part of the reafference outside the span of 30 Gaussian bumps of width $2T/K$. Narrow the bumps and the correlation rises; that is Exercise 4 in miniature. You never told the circuit what the reafference was.
- The prey pulses survive intact. Compute $\langle E\,p_j\rangle$ after learning and confirm it is ≈ 0 for every $j$ while $\langle E \cdot \text{prey}\rangle$ is large: orthogonality principle, measured.
- Now break it. Set `wid = 0.4*T/K` (a basis of narrow, non-overlapping bumps that undersamples the slow component) and re-run: the fast Gaussian is cancelled, the slow one is not. The circuit can only cancel what its basis spans.
- Then flip the sign of `eta` to make the rule Hebbian and watch the weights diverge in the direction of the largest eigenvalue of $C$ — gradient *ascent* on output power. There is no version of this circuit that works with Hebbian plasticity.

---

## 4. Exercises

**★ 1. The sign is forced.**
Derive $w^*=-C^{-1}q$ and the orthogonality principle, and state precisely the assumption under which the residual equals the exafferent signal. Then show that a Hebbian rule $\Delta w = +\eta E p$ has no stable fixed point when $C\succ0$.

<details markdown="1"><summary>Solution</summary>

$J(w)=\tfrac12\langle (R+w^\top p)^2\rangle = \tfrac12\langle R^2\rangle + w^\top q + \tfrac12 w^\top C w$ with $q=\langle pR\rangle$, $C=\langle pp^\top\rangle$. Then $\nabla J = q + Cw = 0 \Rightarrow w^*=-C^{-1}q$. $J$ is convex ($C\succeq0$) so this is the global minimum, unique if $C\succ0$. Orthogonality: $\langle E^*p\rangle = \langle (R + w^{*\top}p)p\rangle = q + Cw^* = 0$.

Since $w^{*\top}p$ is the element of $\operatorname{span}\{p_j\}$ minimising $\langle (R+\cdot)^2\rangle$, it equals $-\Pi R$, so $E^*=(I-\Pi)R$. If (i) $\langle R_{\text{ex}}p_j\rangle=0$ for all $j$ (exafference uncorrelated with the command) and (ii) $R_{\text{re}}\in\operatorname{span}\{p_j\}$ (the basis is rich enough), then $\Pi R = R_{\text{re}}$ and $E^*=R_{\text{ex}}$. Both assumptions are empirical claims about the fish's world and its granule cells respectively — and (ii) is exactly what Kennedy et al. (2014) went looking for.

Hebbian: $\langle\Delta w\rangle = +\eta(q+Cw)$, so $v=w-w^*$ obeys $\langle v_{n+1}\rangle = (I+\eta C)\langle v_n\rangle$, whose eigenvalues all exceed 1 for $C\succ0$. Divergence along the top eigenvector, at rate $(1+\eta\lambda_{\max})^n$. The fixed point exists but is maximally unstable — the rule maximises output power.
</details>

**★★ 2. Learning rates and the eigenvalue spread.**
Give the exact stability condition on $\eta$ and the per-mode time constants. Suppose parallel-fibre activity has eigenvalues spanning two decades. What is the ratio of slowest to fastest cancellation, and what experiment tests it?

<details markdown="1"><summary>Solution</summary>

From $\langle v_{n+1}\rangle=(I-\eta C)\langle v_n\rangle$, convergence requires $|1-\eta\lambda_i|<1$ for every eigenvalue, i.e. $0<\eta<2/\lambda_{\max}$. Mode $i$ decays as $(1-\eta\lambda_i)^n \approx e^{-\eta\lambda_i n}$, so $\tau_i = 1/(\eta\lambda_i)$ update steps (= command cycles).

With $\lambda_{\max}/\lambda_{\min}=100$, the slowest mode takes 100× as many trials as the fastest, and $\eta$ is capped by $\lambda_{\max}$, so you cannot fix it by turning up the learning rate — the condition number of $C$ is a hard limit. (This is why practical LMS uses normalised or whitened variants; a plausible role for the ELL's feed-forward inhibition and for granule-cell decorrelation.)

Test: measure parallel-fibre population covariance $C$ (or estimate its principal directions from granule-cell recordings), then impose reafferent waveforms deliberately aligned with the $i$-th eigenvector and measure cancellation time constants. Prediction: $\tau_i \propto 1/\lambda_i$, with a *single* fitted $\eta$ across all waveforms. That is a one-parameter fit to a whole family of curves — the shape of prediction that can actually be wrong.
</details>

**★★ 3. Speed–precision trade-off.**
Derive the LMS misadjustment $M=\tfrac{\eta}{2}\operatorname{tr}(C)$ (small-$\eta$, independence-assumption version) and interpret it as a design constraint on the fish.

<details markdown="1"><summary>Solution</summary>

Write $v_n=w_n-w^*$ and $E = E^* + v_n^\top p$ with $E^*$ the optimal residual, which is orthogonal to $p$. Then
$$v_{n+1} = v_n - \eta\,p\,(E^* + p^\top v_n) = (I-\eta pp^\top)v_n - \eta p E^*.$$
Taking the covariance $V_n=\langle v_nv_n^\top\rangle$, assuming $p$ independent of $v_n$ (the standard independence assumption) and $\langle pE^*\rangle=0$:
$$V_{n+1} = V_n - \eta(CV_n+V_nC) + \eta^2\langle p p^\top v v^\top pp^\top\rangle + \eta^2 J_{\min}' C,$$
where $J_{\min}'=\langle E^{*2}\rangle$. To first order in $\eta$, drop the third term; the steady state satisfies $CV+VC = \eta J_{\min}' C$, giving $V = \tfrac{\eta}{2}J_{\min}' I$. The excess mean-square error is $\langle (v^\top p)^2\rangle = \operatorname{tr}(CV) = \tfrac{\eta}{2}J_{\min}'\operatorname{tr}(C)$, so the fractional excess (misadjustment) is $M=\tfrac{\eta}{2}\operatorname{tr}(C)$.

Interpretation: $\operatorname{tr}(C)$ is total parallel-fibre power. Adaptation speed $\propto\eta$; residual noise floor $\propto\eta$. A fish that must track a rapidly changing body geometry pays with a permanently noisier percept. Two testable corollaries: (i) manipulations that increase granule-cell drive should *worsen* steady-state cancellation at fixed $\eta$; (ii) species or life stages with more stable electrogenesis should show slower, more precise adaptation.
</details>

**★★ 4. What the circuit cannot cancel.**
Let the basis be the first $m$ Fourier components on the command cycle and let the reafference contain a component at a higher frequency. Compute the residual exactly, and state the general principle.

<details markdown="1"><summary>Solution</summary>

With $p_j$ orthonormal Fourier modes ($C=I$), $\Pi$ is truncation to the first $m$ modes. If $R_{\text{re}}(t) = \sum_k a_k\phi_k(t)$, then $w^{*\top}p = -\sum_{k\le m}a_k\phi_k$ and

$$E^* = \sum_{k>m}a_k\phi_k, \qquad \langle E^{*2}\rangle = \sum_{k>m}a_k^2 .$$

The unrepresentable power passes through untouched — it is *neither* cancelled *nor* amplified, and it is indistinguishable, to the circuit, from prey.

General principle: an adaptive filter's residual is the projection of its input onto the orthogonal complement of the span of its basis. This turns a structural fact about granule cells into a functional prediction: characterise the granule-cell temporal basis, compute its span, and you can predict *which* artificial reafferent waveforms the ELL will fail to cancel — and those failures should look, behaviourally, like phantom prey. That is a rare thing: a prediction of a specific, engineered illusion.
</details>

**★★★ 5. A non-delta plasticity kernel.**
The measured STDP window has width $\sim50$ ms. Model the rule as $\dot w_j = -\eta\,(K * p_j)(t)\,E(t)$ with kernel $K$. When is this still exact gradient descent on $J$? When is it not, and what goes wrong?

<details markdown="1"><summary>Solution</summary>

Averaging, $\langle\dot w_j\rangle = -\eta\langle \tilde p_j E\rangle$ with $\tilde p = K*p$. The fixed point is $\langle \tilde p_j E\rangle = 0$ for all $j$, i.e. $\tilde{q} + \tilde{C}_{\times}w = 0$ where $\tilde q = \langle\tilde p R\rangle$ and $(\tilde C_\times)_{jk}=\langle \tilde p_j p_k\rangle$. This is a *cross*-correlation normal equation, not the symmetric one.

Two clean cases. (a) If $K$ acts as a delay-and-scale on each basis function, $\tilde p_j = c\,p_{\sigma(j)}$ for a permutation $\sigma$ (true if the basis is a delay line and $K$ is narrow), then $\operatorname{span}\{\tilde p\}=\operatorname{span}\{p\}$ and the fixed point is the true optimum: the residual is still orthogonal to the span, so $E^*=(I-\Pi)R$ regardless. Timing offsets change the *route* to the solution, not the solution. (b) If $K$ smears across basis functions so that $\operatorname{span}\{\tilde p\}\subsetneq\operatorname{span}\{p\}$, the fixed point solves a projected problem in the smaller span, and residual power increases; and if $\tilde C_\times$ is non-normal or has eigenvalues with negative real part, the dynamics can be unstable or oscillatory even for small $\eta$ — the "wrong-sign" region of the STDP window matters, which is precisely why Roberts & Bell had to check stability numerically rather than assert it.

Practical upshot for the reader: measure the STDP kernel and the granule basis, compute $\tilde C_\times$, and check its spectrum. If it has eigenvalues in the left half plane the circuit as described cannot work, and something else (e.g. the non-associative term, or interneuron dynamics) is doing stabilising work you have not modelled.
</details>

**★★★ 6. LMS as a degraded Kalman filter.**
Suppose the true reafference drifts: $w^{\text{true}}_{n+1}=w^{\text{true}}_n+\xi_n$, $\xi\sim\mathcal N(0,\sigma_w^2 I)$, with observation noise variance $\sigma_\varepsilon^2$. In the scalar case ($K=1$, $p\equiv1$), derive the steady-state Kalman gain and hence the optimal constant learning rate. Interpret.

<details markdown="1"><summary>Solution</summary>

Kalman recursion for the scalar random walk: prior variance $S_n = \Sigma_{n-1}+\sigma_w^2$; gain $k_n = S_n/(S_n+\sigma_\varepsilon^2)$; posterior $\Sigma_n = (1-k_n)S_n = S_n\sigma_\varepsilon^2/(S_n+\sigma_\varepsilon^2)$. At steady state $S=\Sigma+\sigma_w^2$:

$$S-\sigma_w^2 = \frac{S\sigma_\varepsilon^2}{S+\sigma_\varepsilon^2} \;\Longrightarrow\; S^2-\sigma_w^2S-\sigma_w^2\sigma_\varepsilon^2=0 \;\Longrightarrow\; S=\frac{\sigma_w^2+\sqrt{\sigma_w^4+4\sigma_w^2\sigma_\varepsilon^2}}{2}.$$

Let $\rho=\sigma_w^2/\sigma_\varepsilon^2$ (drift-to-noise ratio). For $\rho\ll1$, $S\approx\sigma_w\sigma_\varepsilon$ and

$$\eta^* = k = \frac{S}{S+\sigma_\varepsilon^2} \approx \frac{\sigma_w}{\sigma_\varepsilon} = \sqrt{\rho}.$$

So the optimal LMS rate is the **square root** of the drift-to-noise ratio — sub-linear, which is why well-tuned adaptive filters look sluggish relative to naive intuition. LMS with fixed $\eta$ is exactly the constant-gain approximation to this filter; it is optimal when the environment's non-stationarity is itself stationary, and it is the price of not maintaining a covariance matrix. (The multivariate version replaces the scalar by mode-wise gains $\propto\sqrt{\sigma_{w,i}^2/\sigma_\varepsilon^2}$, which is a normalised-LMS rule — another thing the circuit could plausibly implement via granule-cell gain control.)

Prediction: perturb the fish's electrogenesis with a drift of controlled magnitude and measure the effective learning rate; it should scale as the square root of the drift variance, not linearly. Nobody has done this. It is a good experiment.
</details>

---

## 5. Reading path

- **Bell (1981)**, *An efference copy which is modified by reafferent input* — read it for: the original negative-image demonstration, and the fact that the cancelled pattern can be arbitrary.
- **Widrow & Hoff (1960)**, *Adaptive switching circuits* (and Widrow & Stearns, *Adaptive Signal Processing*, 1985) — read it for: LMS, convergence in the mean, and misadjustment, derived by engineers with no interest in fish.
- **Fujita (1982)**, *Adaptive filter model of the cerebellum* — read it for: the first clean statement that a cerebellum-like circuit is a filter with a learned weight vector.
- **Bell, Han, Sugawara & Grant (1997)**, *Synaptic plasticity in a cerebellum-like structure depends on temporal order* — read it for: the anti-Hebbian STDP window measured directly; the sign of the gradient.
- **Roberts & Bell (2000)**, *Computational consequences of temporally asymmetric learning rules: II. Sensory image cancellation* — read it for: what the finite STDP kernel does to stability and to the fixed point.
- **Porrill, Dean & Stone (2004)**, *Recurrent cerebellar architecture solves the motor-error problem* — read it for: how decorrelation substitutes for a gradient you cannot compute.
- **Bell, Han & Sawtell (2008)**, *Cerebellum-like structures and their implications for cerebellar function* — read it for: the comparative anatomy that licenses generalising from mormyrids to cerebellum.
- **Dean, Porrill, Ekerot & Jörntell (2010)**, *The cerebellar microcircuit as an adaptive filter: experimental and computational evidence* — read it for: the best single audit of the adaptive-filter reading of cerebellum.
- **Kennedy, Wayne, Kaifosh, Alviña, Abbott & Sawtell (2014)**, *A temporal basis for predicting the sensory consequences of motor commands in an electric fish* — read it for: the granule-cell basis, i.e. the $p(t)$ in every equation above, measured.
- **Muller, Zadina, Abbott & Sawtell (2019)** — read it for: a multi-layer, continually-learning extension of the ELL model and what changes when the filter is deep.
- **Poulet & Hedwig (2006)**, *The cellular basis of a corollary discharge* — read it for: the same template in an insect, with a single identified interneuron.
- **Schneider, Nelson & Mooney (2014)**, *A synaptic and circuit basis for corollary discharge in the auditory cortex* — read it for: the mammalian cortical instance, learned and sound-specific.

---

## 6. Open problems and what would settle them

**Is the granule-cell basis adapted to the reafference statistics, or generic?** §2.4 says learning speed is governed by the spectrum of $C$; an evolved system should shape $C$ to match the reafference it must cancel. Settling move: compare the measured principal directions of granule-cell activity with the principal components of the natural reafferent waveform ensemble. Strong alignment is evidence for an optimised basis; a generic (e.g. approximately isotropic, random-expansion) basis is evidence that the ELL follows the same design as the mushroom body — which would be a beautiful, and unifying, negative result. Note the connection to `./C2-transient-synchrony.md`: random expansion versus tuned basis is the same argument in a different animal.

**How is the effective learning rate set, and is it optimal?** Exercise 6 gives a $\sqrt{\rho}$ prediction that has never been tested. Settling move: controlled perturbation of reafference variability plus measurement of adaptation time constants across conditions.

**Does the deep (multi-layer) version still compute a projection?** Once the filter has hidden layers, the residual is no longer an orthogonal projection and the "innovation" interpretation weakens. Settling move: in the Muller et al. architecture, ask whether the measured residual is still orthogonal to the parallel-fibre basis; if it is not, the algorithm has changed and needs a new level-2 statement.

**What reads the innovation?** The negative-image story explains what the ELL removes but not what the next stage does with what remains. Until someone shows a downstream population whose behaviour depends on the residual in the way the theory requires, the conversion is a complete account of an encoder and no account at all of a computation. This is the same gap as in `./C7-dopamine-and-td-learning.md` §6 ("is the distributional code read out?"), and it is the most common way a level-2 claim in neuroscience is quietly incomplete.

*Next:* `./C10-evidence-accumulation.md`, where the algorithm is fixed by an optimality theorem rather than by a plasticity rule — and where the causal evidence for the circuit implementing it turns out to be much weaker than the textbooks say.
