---
title: "C3 · Coincidence detection"
parent: Case studies
nav_order: 3
---

# C3 — Coincidence detection + delay lines → cross-correlation

> **Circuit.** The barn owl's nucleus laminaris: bilateral axonal delay lines from nucleus magnocellularis interdigitate along an array of bipolar coincidence-detecting neurons, so that each NL neuron fires maximally for one interaural time difference, and ITD is read out as a place code.
> **Primitive extracted.** A running cross-correlation of the two ears' spike trains, sampled in parallel at a bank of internal delays — the Jeffress model, and one of the very few cases in neuroscience where an algorithm was proposed from first principles and the anatomy later found to match.
> **Status.** Well established in birds; **partly wrong as a general theory**, because the mammalian medial superior olive appears not to use a place code at all. This unit is as much about the epistemology of that failure as about the algorithm.
> **Structures thread.** `structures/S-04-manifold-capacity.md` — the Fisher-information argument in §2.5 is the low-dimensional ancestor of the population-geometry tools developed there; the Harper–McAlpine result is a worked example of a coding scheme being optimal only on part of a parameter manifold. See `structures/README.md`.

## 1. The phenomenon

A barn owl striking a mouse in total darkness localises it to about 1–2 degrees in azimuth. Azimuth is carried almost entirely by the interaural time difference: at the owl's head width the full physiological ITD range is roughly $\pm 250\ \mu$s, and behavioural and neural thresholds are on the order of $10\ \mu$s or less. Neurons have millisecond membrane time constants and spike-timing jitter in the hundreds of microseconds. A system built from sloppy components resolves a fraction of the sloppiness — *hyperacuity*, and the first thing any account has to explain.

Jeffress (1948) proposed the mechanism before any of the relevant physiology existed: an array of coincidence detectors, each receiving input from the two ears through axons of *different lengths*, so that each detector is "tuned" to the ITD that exactly compensates its internal delay. ITD becomes a position along the array — a place code — and the operation is a cross-correlation computed in parallel.

Carr & Konishi (1988, 1990) then found it. In nucleus laminaris, magnocellularis axons from the contralateral side enter from below and run dorsally while ipsilateral axons enter from above and run ventrally; conduction delay increases monotonically along each, and the two gradients oppose. NL neurons along the dorsoventral axis are tuned to systematically varying ITDs. Knudsen & Konishi (1978) had already shown a map of auditory space in the inferior colliculus, and Moiseff & Konishi (1981) that ITD is the cue. It is the cleanest circuit-to-algorithm correspondence in the field.

Then the trouble. In small mammals — gerbil, guinea pig — McAlpine, Jiang & Palmer (2001) found that the ITDs to which MSO/IC neurons are best tuned mostly lie *outside* the animal's physiological ITD range, clustering near a quarter of the characteristic-frequency period. There is no place code; ITD is read from the *slope* of broadly tuned responses, and plausibly from the difference between two large hemispheric populations. Brand et al. (2002) showed that precisely timed glycinergic inhibition is necessary for microsecond ITD tuning in the gerbil MSO, a mechanism with no counterpart in the Jeffress picture.

## 2. The conversion

### 2.1 From coincidence detectors to a correlator: do it properly

Model each ear's input to a given frequency channel as an inhomogeneous Poisson process. Let $\lambda_L(t)$, $\lambda_R(t)$ be the instantaneous rates of the left and right afferent populations after cochlear filtering at characteristic frequency $f$, and let $dN_L, dN_R$ be the corresponding counting measures.

A coincidence detector with internal delay $\delta$ and a symmetric window $w(\cdot)$ (normalised, width $\sim$ its membrane integration time) computes

$$
C_\delta \;=\; \iint w\big((u-\delta) - v\big)\, dN_L(u)\, dN_R(v),
$$

i.e. it counts pairs of spikes — one from each side — that arrive within a window of each other *after* the ipsilateral line has been delayed by $\delta$. Take expectations. Because the two ears' spike trains are conditionally independent given the rates,

$$
\mathbb{E}[C_\delta] \;=\; \iint w(u-\delta-v)\,\lambda_L(u)\lambda_R(v)\,du\,dv
\;=\; \int w(\tau - \delta)\, R_{LR}(\tau)\, d\tau,
$$

where

$$
R_{LR}(\tau) \;=\; \int \lambda_L(u)\,\lambda_R(u-\tau)\, du
$$

is the **cross-correlation function of the two ears' drive**. So:

$$
\boxed{\ \mathbb{E}[C_\delta] = \big(w * R_{LR}\big)(\delta)\ }
$$

The array of NL neurons evaluates the cross-correlation function, smoothed by the coincidence window, at the discrete set of internal delays $\{\delta(x)\}$ realised by the delay lines. That is the entire conversion, and it is worth noticing how little it assumes: only conditional independence and a multiplicative-then-integrate detector. It does not assume linearity, tuning curves, or Gaussianity.

Now let the sound source impose an ITD $\Delta$, so $\lambda_R(t) = \lambda_L(t-\Delta)$. Then $R_{LR}(\tau) = R_{LL}(\tau - \Delta)$, the *autocorrelation* of the monaural drive shifted to $\Delta$, and

$$
\mathbb{E}[C_\delta] = (w * R_{LL})(\delta - \Delta),
$$

peaked at $\delta = \Delta$. The position of the maximally active neuron reports the ITD: place code. The width of the peak is the width of $R_{LL}$ convolved with $w$ — the stimulus autocorrelation width, not the neuron's time constant.

### 2.2 Why the peak is periodic, and what the owl does about it

Cochlear filtering makes each channel narrowband. For a channel with centre frequency $f$ and bandwidth $B$, the monaural autocorrelation is

$$
R_{LL}(\tau) \;\propto\; e^{-\pi^2 B^2\tau^2}\cos(2\pi f\tau)
$$

for a Gaussian-shaped filter — an oscillation at $f$ under an envelope of width $\sim 1/B$. Therefore $\mathbb{E}[C_\delta]$ has side peaks at $\delta = \Delta \pm k/f$, all of comparable height whenever $B \ll f$ (which is exactly the case: auditory filters have roughly constant $Q$, so $1/B \approx Q/f$ with $Q \sim 5$–$10$, giving 5–10 visible side peaks).

The ambiguity is real whenever the physiological ITD range exceeds one period:

$$
2\,\Delta_{\max} > 1/f \quad\Longleftrightarrow\quad f > \frac{1}{2\Delta_{\max}} .
$$

For the owl, $\Delta_{\max} \approx 250\ \mu$s, so ambiguity begins at $2$ kHz — and the owl uses ITD up to $8$ kHz. Every high-frequency channel is therefore *individually* ambiguous. The resolution is cross-frequency convergence: the side peaks sit at $\Delta \pm k/f$, whose positions depend on $f$, while the true peak sits at $\Delta$ regardless. Summing (or multiplying) across frequency channels leaves only the true peak standing. This is the ICc→ICx transformation in the owl (Konishi, 2003), and Peña & Konishi (2001) showed that ICx space-specific neurons' receptive fields are well described as a *multiplication* of ITD and interaural-level tuning — multiplicative combination, which suppresses side peaks far more effectively than addition. §3 demonstrates the effect in ten lines.

### 2.3 Sampling: what the delay-line spacing has to satisfy

The array samples $(w*R_{LL})(\delta-\Delta)$ at discrete $\delta$. Treated as a signal in $\delta$, its bandwidth is set by the cochlear channel: the correlation function contains frequencies up to about $f + B/2$. Nyquist requires

$$
\Delta\delta \;<\; \frac{1}{2(f + B/2)} \;\approx\; \frac{1}{2f},
$$

i.e. **the delay-line spacing must be finer than half the period of the channel's characteristic frequency.** At $f = 6$ kHz that is $83\ \mu$s — an easy anatomical requirement, corresponding to a fraction of a millimetre of unmyelinated axon. Undersampling would not merely blur the estimate, it would *alias*: a true ITD would masquerade as a different one, producing systematic localisation errors rather than noise. That the owl makes no such errors is a (weak) confirmation of adequate sampling.

Note carefully what Nyquist does *not* say: it says nothing about the precision with which $\Delta$ can be estimated. That is the next section, and confusing the two is the standard error in discussions of hyperacuity.

### 2.4 Hyperacuity: the calculation

The claim to be defeated is "neurons with $100\ \mu$s jitter cannot resolve $10\ \mu$s." The answer is that estimation precision is not sampling precision: with a smooth tuning curve and many neurons and spikes, the population can localise the *peak* far better than the peak is *wide*.

Take a single NL/IC neuron in a channel at frequency $f$ with cyclic ITD tuning

$$
\bar r(\Delta) \;=\; r_0\big[1 + \cos\big(2\pi f(\Delta - \delta)\big)\big],
$$

and Poisson spike counts over a window $T$, so the variance equals the mean, $\sigma^2 = \bar r T$. The Fisher information about $\Delta$ carried by this neuron is

$$
I_F(\Delta) \;=\; \frac{\big(\partial_\Delta \bar r\, T\big)^2}{\bar r T}
\;=\; \frac{T\big(2\pi f\, r_0 \sin\phi\big)^2}{r_0(1+\cos\phi)}, \qquad \phi = 2\pi f(\Delta-\delta).
$$

Average over an array whose best delays $\delta$ tile the cycle uniformly:

$$
\big\langle I_F\big\rangle_\phi = (2\pi f)^2 r_0 T \cdot \frac{1}{2\pi}\int_0^{2\pi}\frac{\sin^2\phi}{1+\cos\phi}d\phi
= (2\pi f)^2 r_0 T \cdot \frac{1}{2\pi}\int_0^{2\pi}(1-\cos\phi)\,d\phi
= (2\pi f)^2\, r_0 T ,
$$

using $\sin^2\phi/(1+\cos\phi) = 1-\cos\phi$. So each neuron contributes

$$
\boxed{\ I_F = (2\pi f)^2 \times (\text{spikes counted})\ }
$$

and $M$ independent neurons give $\sigma_\Delta \ge 1/\sqrt{M I_F}$ by Cramér–Rao. Put numbers in: $f = 5$ kHz, $r_0 T = 20$ spikes per neuron.

$$
I_F = (2\pi \cdot 5000)^2 \times 20 = 1.97\times10^{10}\ \mathrm{s^{-2}},\qquad
\sigma_\Delta = \frac{1}{\sqrt{M\,I_F}} =
\begin{cases}
7.1\ \mu\mathrm{s}, & M=1\\
2.3\ \mu\mathrm{s}, & M=10\\
0.71\ \mu\mathrm{s}, & M=100 .
\end{cases}
$$

Hyperacuity is not mysterious; it falls out of $(2\pi f)^2$. The single most important structural fact is that **Fisher information scales with $f^2$**, so precision improves quadratically with characteristic frequency. This is why the owl pushed ITD processing to 8 kHz — and, as we are about to see, why mammals that cannot phase-lock above ~1.5 kHz had to solve the problem differently.

### 2.5 Where the beautiful story is species-specific: Grothe, McAlpine, and the slope code

The Jeffress model makes three linked predictions: (i) internal delays tile the physiological range; (ii) best ITDs form an orderly map; (iii) the ITD estimate is the *location* of the peak of population activity. In small mammals, all three appear to fail.

McAlpine, Jiang & Palmer (2001) found that best delays in guinea pig IC cluster around $\delta_{\text{best}} \approx 1/(8f)$ — a quarter cycle — which for low CFs lies *beyond* the physiological range. The physiological ITDs therefore fall on the steep flank of the tuning curve, not near its peak. Harper & McAlpine (2004) then asked the normative question properly: for a population of $M$ neurons with cyclic tuning, what arrangement of best delays maximises Fisher information about ITD, given the head size (which sets $\Delta_{\max}$) and the CF?

The answer is a phase transition in parameter space. Redo §2.4's calculation but now let $\bar r(\Delta) = r_0[1+\cos(2\pi f(\Delta-\delta))]$ over a *restricted* range $|\Delta| \le \Delta_{\max}$. Fisher information at $\Delta$ from a neuron with best delay $\delta$ is $\propto \sin^2\phi/(1+\cos\phi) = 1-\cos\phi$, maximised at $\phi = \pi$ — i.e. the *most informative* neuron for a given ITD is the one whose best delay is half a cycle away. Two regimes follow:

- **$2\Delta_{\max} f \gg 1$** (large head relative to the period — the owl at high CF). Half-cycle offsets are available *within* the physiological range for every $\Delta$, and the informative population is different for each $\Delta$. Tiling the range with best delays is optimal: **place code**.
- **$2\Delta_{\max} f \lesssim 1$** (small head, low CF — gerbil, human at low frequency). A half-cycle offset lies outside the physiological range. The best you can do is push all best delays to the edge, $\delta \approx \pm 1/(8f)$ to $\pm 1/(4f)$, and read the ITD from the *difference* between the two hemispheric populations: **two-channel slope code**.

This is the lesson of the unit, and it is more important than the algorithm. Both codes solve the *same* problem — maximise Fisher information about ITD — and the optimal implementation flips as a function of $2\Delta_{\max}f$, a parameter that has nothing to do with the computation and everything to do with the animal's head and its cochlea. The Jeffress model is not "wrong"; it is the large-$2\Delta_{\max}f$ limit of a more general normative theory, and it was mistaken for a universal only because the first system examined happened to be at that limit.

The mechanistic differences follow. In mammals, internal delays need not come from axonal length at all: candidates include cochlear-place disparities between the two ears, asymmetric synaptic dynamics (Jercog et al., 2010), and — the best-supported — precisely timed contralateral glycinergic inhibition that shifts the effective coincidence time (Brand et al., 2002; Pecka et al., 2008). Franken and colleagues' in-vivo whole-cell work in the gerbil MSO further showed that the internal delays behave more like *phase* delays than pure time delays, which is a qualitatively different object. Joris & Yin (2007) is the honest survey of what "internal delay" can even mean.

### 2.6 Connection to C2: the same primitive, two ways of using it

Strip the delay lines away and what remains is: *multiply two spike trains and integrate over a short window*. That is precisely the Kenyon cell's operation in C2's hypothesis (d), and the mathematics is shared. Compare the two SNR expressions:

- Owl NL, §2.1: $\mathbb{E}[C_\delta] = (w * R_{LR})(\delta)$ — the detector reports the correlation *at a specific lag*, and the lag is the variable of interest.
- Locust KC, C2 §2.4: $d' \approx K/\sqrt{B\nu\sigma}$ — the detector reports *whether* a synchronous packet exists at lag zero, and the packet's membership is the variable of interest.

Both exploit the same fact: coincidence within a window $w$ suppresses uncorrelated background as $\sqrt{w}$ while preserving correlated input, so the achievable SNR improves as $w$ shrinks toward the jitter. The owl externalises the lag as anatomy; the locust externalises the *identity* of the coincident set as which KC fires. One circuit motif, two computational uses, distinguished entirely by what the array of detectors is indexed by. If you take one structural idea from Part 2, take this one: **the primitive is coincidence; the algorithm is determined by what you vary across the detector array.**

## 3. Worked example / model to build

Phase-locked binaural spike trains, an array of internal delays, and the frequency-ambiguity problem with its resolution.

```python
import numpy as np
dt, T = 5e-6, 0.20                       # 5 us bins, 200 ms
n  = int(T/dt); tt = np.arange(n)*dt
ITD = 80e-6                              # true interaural time difference (s)
CFs = [2000., 3000., 4000., 5000., 6000.]
NF, RATE, VS, W = 30, 250., 0.95, 30e-6  # fibres/side/CF, rate, vector strength, window

def fibres(cf, shift, seed):
    r = np.random.default_rng(seed)
    lam = RATE*(1 + VS*np.cos(2*np.pi*cf*(tt - shift)))
    return (r.random((NF, n)) < lam*dt).sum(0).astype(float)

def xcorr(a, b, maxlag):
    L = int(maxlag/dt); m = 1 << int(np.ceil(np.log2(n + L)))+1
    A, B = np.fft.rfft(a, m), np.fft.rfft(b, m)
    c = np.fft.irfft(A*np.conj(B), m)
    return np.concatenate([c[-L:], c[:L+1]])          # lags -L..L

lags = np.arange(-int(300e-6/dt), int(300e-6/dt)+1)*dt
win  = np.ones(int(W/dt)); win /= win.sum()
tot  = np.zeros_like(lags)
for i, cf in enumerate(CFs):
    l = fibres(cf, 0.0, 10+i); r = fibres(cf, ITD, 100+i)
    c = np.convolve(xcorr(r, l, 300e-6), win, mode='same')
    tot += c
    print("CF %4.0f Hz: single-channel peak at %+6.1f us  (period %5.1f us)"
          % (cf, lags[np.argmax(c)]*1e6, 1e6/cf))
print("\nsummed across CFs: peak at %+6.1f us   (true ITD %+6.1f us)"
      % (lags[np.argmax(tot)]*1e6, ITD*1e6))

# ---- hyperacuity: Fisher information of a single 5 kHz channel ---------------
f, r0Tspk = 5000., 20.                   # 20 spikes in the analysis window
IF = (2*np.pi*f)**2 * r0Tspk             # per neuron, s^-2  (see Section 2.4)
for Mn in [1, 10, 100]:
    print("M=%3d neurons -> ITD discrimination limit %.2f us" % (Mn, 1e6/np.sqrt(Mn*IF)))
```

**What to look for.** Three of the five single-frequency channels put their peak in the *wrong* place — at $2$ kHz you should get $+75\ \mu$s (correct), at $3$ kHz $-255\ \mu$s, at $4$ kHz $+80$ (correct), at $5$ kHz $+285$, at $6$ kHz $-80$. Check each wrong answer against the period printed alongside: every error is the true ITD plus or minus an integer number of periods. That is §2.2's ambiguity, seen directly. The summed correlogram lands on $+80\ \mu$s exactly. Then the Fisher calculation reproduces $7.1 / 2.3 / 0.71\ \mu$s for $M = 1/10/100$.

**Things to do to it.** (i) Replace the sum across CFs with a product and compare side-lobe suppression — this is the Peña & Konishi multiplication. (ii) Shrink the number of CFs to one and try to recover ITD anyway; you will find you can, provided $|{\rm ITD}| < 1/(2f)$, which is exactly the $\pi$-limit. (iii) Restrict the delay array to $|\delta| \le 130\ \mu$s (a gerbil-sized head) at CF $= 500$ Hz and see that the peak never falls inside the array — then implement the two-channel slope readout and check it still works. That is §2.5 in code, and it is the most valuable of the three.

## 4. Exercises

**E1 (★).** Show that a coincidence detector receiving *independent* Poisson inputs at rates $\lambda_L, \lambda_R$ with window $w$ produces coincidences at rate $2\lambda_L\lambda_R w$, and use this to compute the "chance" floor of the correlogram in §3.

<details markdown="1"><summary>Solution</summary>

Fix a left spike at time $u$. The expected number of right spikes within $\pm w/2$ of it is $\lambda_R w$ (window of total width $w$; if instead you define $|u-v|<w$ the width is $2w$, which is the source of the factor 2). Over a duration $T$ there are $\lambda_L T$ left spikes, so the expected coincidence count is $\lambda_L T \cdot \lambda_R \cdot 2w$, i.e. a rate

$$\mathcal{R}_{\text{chance}} = 2\lambda_L\lambda_R w .$$

In §3, per channel: $\lambda = 30\ \text{fibres}\times 250\ \mathrm{s^{-1}} = 7500\ \mathrm{s^{-1}}$ per side, $w = 30\ \mu$s, $T = 0.2$ s. Chance coincidences $= 2\times 7500^2 \times 30\times10^{-6}\times 0.2 = 675$ per lag bin's worth of window. The *modulated* part is what carries the signal: with vector strength $\mathrm{VS}=0.95$ the rate is $\lambda(1+0.95\cos)$, so the correlogram's oscillatory component has amplitude $\tfrac12 (0.95)^2 \approx 0.45$ of the flat term. Signal-to-chance is therefore $\approx 0.45$, and the fluctuation on the flat term is $\sqrt{675}\approx 26$, i.e. 4% — comfortably detectable. This is why the demo works with only 30 fibres.

The general lesson: the correlogram sits on a large pedestal set by the product of the rates, and the ITD information is in the modulation on top of it. Any mechanism that reduces the pedestal without reducing the modulation (e.g. a threshold, or inhibition) improves the code — which is one reading of the glycinergic inhibition in MSO.
</details>

**E2 (★★).** Derive the "$\pi$-limit." Show that for a neuron with cyclic ITD tuning at CF $f$, the Fisher information per neuron, averaged over a *restricted* ITD range $|\Delta|\le\Delta_{\max}$, is maximised by a best delay near a quarter cycle when $2\Delta_{\max}f \lesssim 1$.

<details markdown="1"><summary>Solution</summary>

From §2.4, the Fisher information at ITD $\Delta$ from a neuron with best delay $\delta$ is
$$I_F(\Delta;\delta) \propto 1 - \cos\big(2\pi f(\Delta-\delta)\big),$$
increasing monotonically in $|\Delta - \delta|$ up to half a cycle.

Suppose ITDs are uniform on $[-\Delta_{\max},\Delta_{\max}]$. The neuron's average information is
$$\bar I(\delta) \propto \frac{1}{2\Delta_{\max}}\int_{-\Delta_{\max}}^{\Delta_{\max}}\!\!\big[1-\cos(2\pi f(\Delta-\delta))\big]d\Delta
= 1 - \frac{\sin(2\pi f\Delta_{\max})}{2\pi f\Delta_{\max}}\cos(2\pi f\delta).$$

Write $a = 2\pi f\Delta_{\max}$. Then $\bar I(\delta) = 1 - \mathrm{sinc}(a)\cos(2\pi f\delta)$, where $\mathrm{sinc}(a)=\sin a/a$.

- If $a < \pi$ (i.e. $2\Delta_{\max}f < 1$): $\mathrm{sinc}(a) > 0$, so $\bar I$ is maximised by making $\cos(2\pi f\delta)$ as negative as possible, i.e. $\delta = \pm 1/(2f)$ — half a cycle, *outside* the range. Since best delays are physically constrained (and since one also needs a *monotone*, unambiguous readout over the range, which fails once $|\delta|$ exceeds a quarter cycle in the wrong direction), the constrained optimum sits at the edge of the usable region, near $|\delta| \approx 1/(4f)$ — and empirically McAlpine et al. report a clustering nearer $1/(8f)$, which is what you get when you additionally require the *whole* physiological range to lie on a monotone flank rather than optimising pointwise.
- If $a$ is large ($2\Delta_{\max}f \gg 1$): $\mathrm{sinc}(a)\to 0$, so $\bar I(\delta)$ is nearly independent of $\delta$ — every best delay is equally good *on average*, and the population should instead be arranged to maximise information *at each particular $\Delta$*, which requires tiling. Place code.

So the transition is governed by $\mathrm{sinc}(2\pi f\Delta_{\max})$ becoming negligible, i.e. by $2\Delta_{\max}f$ crossing order unity. Owl at 5 kHz: $2\times250\times10^{-6}\times5000 = 2.5$. Gerbil at 500 Hz with $\Delta_{\max}=130\ \mu$s: $0.13$. Human at 500 Hz with $\Delta_{\max}=700\ \mu$s: $0.7$ — intermediate, which is why the human literature is contested.
</details>

**E3 (★★).** The owl's NL neurons have membrane time constants of order a millisecond, yet the coincidence window must be tens of microseconds. Resolve the apparent contradiction, and state the two mechanisms that make it work.

<details markdown="1"><summary>Solution</summary>

The coincidence window is not the membrane time constant; it is the time over which the *summed* EPSP exceeds threshold, and that can be much shorter than $\tau_m$ if (i) the EPSPs are fast and (ii) the neuron is operating far up a steep threshold nonlinearity.

Quantitatively: let each input produce an EPSP $g(t)$ of half-width $\theta$, and let $K$ synchronous inputs be required to reach threshold out of $N$ available. The probability of reaching threshold as a function of relative arrival jitter falls off over a timescale $\sim\theta/\sqrt{K}$ — *narrower* than the EPSP itself, because $K$ near-simultaneous events summing to a sharp peak lose height rapidly as they spread. The effective window therefore shrinks as the required coincidence order grows.

Two mechanisms in the real system:
1. **Very fast synaptic and membrane kinetics.** NL/MSO neurons express low-voltage-activated potassium currents (Kv1) that produce a strongly sublinear, differentiating response: the cell effectively computes $d/dt$ of its input, so slow summation is actively cancelled and only fast coincident input survives. The effective $\tau$ is an order of magnitude below the passive $\tau_m$.
2. **Bipolar dendrites.** Agmon-Snir, Carr & Rinzel (1998) showed that segregating the two ears' inputs onto separate dendrites *improves* coincidence detection: inputs from one ear alone saturate locally (driving the local membrane toward the synaptic reversal potential and reducing driving force), so monaural input is penalised relative to binaural input. The dendritic architecture implements a soft multiplication rather than a sum. This is a beautiful case of morphology being the algorithm.

Note the shape of the answer: the "contradiction" dissolves once you stop identifying the computational window with a passive time constant. The same move is needed in C2, where KC integration windows are far shorter than KC membrane time constants for exactly the same reasons.
</details>

**E4 (★★, computational).** Using the §3 code, implement the gerbil case: CF $=500$ Hz, $\Delta_{\max}=130\ \mu$s, and a two-channel readout. Compare its Fisher information to a Jeffress array with the same number of neurons.

<details markdown="1"><summary>Solution</summary>

Setup. Two populations, "left-preferring" with best delay $\delta_+ = +1/(8f) = +250\ \mu$s and "right-preferring" with $\delta_- = -250\ \mu$s, $M/2$ neurons each. Readout $\hat\Delta \propto (r_+ - r_-)$.

Fisher information of the two-channel scheme. Using $I_F \propto 1-\cos(2\pi f(\Delta-\delta))$ per neuron with $f=500$ Hz, $2\pi f = 3142\ \mathrm{s^{-1}}$: over $|\Delta|\le130\ \mu$s, $2\pi f(\Delta - \delta_+)$ ranges over $[-1.19, -0.38]$ rad, so $1-\cos \in [0.07, 0.63]$, mean $\approx 0.33$. The two populations contribute symmetrically, so
$$I_F^{\text{2ch}} \approx M\times 0.33 \times (2\pi f)^2 r_0T .$$

Jeffress array with best delays uniform on $[-130, +130]\ \mu$s. Then $2\pi f(\Delta-\delta)$ ranges over $[-0.82, 0.82]$ rad and $\langle 1-\cos\rangle$ over the joint uniform distribution is $1 - \mathrm{sinc}^2(0.82)/1 \approx 1 - (0.731)^2$… computing directly: $\langle\cos(2\pi f(\Delta-\delta))\rangle = \mathrm{sinc}(a)^2$ with $a = 2\pi f\Delta_{\max}=0.408$ rad for each of the two independent uniforms — $\mathrm{sinc}(0.408)=0.9725$, squared $=0.9458$. So $\langle 1-\cos\rangle = 0.054$ and
$$I_F^{\text{Jeffress}} \approx M\times 0.054 \times (2\pi f)^2 r_0T .$$

Ratio: $0.33/0.054 \approx 6$. **The two-channel arrangement carries about six times more Fisher information per neuron than a Jeffress array at gerbil parameters.** Numerically, with $r_0T=20$ and $M=100$: $\sigma_\Delta^{\text{2ch}} = 1/\sqrt{100\times0.33\times(3142)^2\times20} = 22\ \mu$s versus $54\ \mu$s for Jeffress — and gerbil behavioural thresholds are in the tens of microseconds, consistent with the former.

The code check: simulate both, decode by maximum likelihood, and plot RMS error against $\Delta$. You should find the Jeffress array is not merely worse on average but *catastrophically* worse near $\Delta = 0$, where every neuron sits near its tuning peak and the derivative vanishes. That is the intuitive statement of the result: a place code wastes its most numerous neurons at exactly the ITD the animal encounters most often.
</details>

**E5 (★★).** Show that the "$f^2$ scaling" of Fisher information implies a specific relation between an animal's head width and the frequency range it should use for ITD, and check it against owl and human.

<details markdown="1"><summary>Solution</summary>

Two competing pressures.

*Upward:* $I_F \propto f^2$ per spike, so higher CF is better — precision improves linearly in $f$.

*Downward:* phase-locking degrades above a species-specific cutoff $f_\ell$ (about $9$ kHz in the owl, ~$1.4$ kHz in humans and most mammals), because synaptic and hair-cell jitter $\sigma_j$ smears the phase. Vector strength falls roughly as $\mathrm{VS}(f) = \exp(-2\pi^2f^2\sigma_j^2)$ for Gaussian jitter, and the usable information scales as $I_F \propto f^2\,\mathrm{VS}(f)^2 = f^2 e^{-4\pi^2f^2\sigma_j^2}$.

Maximise: $d/df[2\ln f - 4\pi^2f^2\sigma_j^2] = 2/f - 8\pi^2f\sigma_j^2 = 0 \Rightarrow f^\star = 1/(2\pi\sigma_j)$. So the optimal ITD frequency is set entirely by temporal jitter, *not* by head size. Owl: $\sigma_j \approx 25\ \mu$s gives $f^\star \approx 6.4$ kHz — right in the owl's ITD band. Human: $\sigma_j \approx 100\ \mu$s gives $f^\star \approx 1.6$ kHz — right at the human phase-locking limit.

Head size enters separately, through $\Delta_{\max}$, and it decides only the *code* (§2.5), via $2\Delta_{\max}f^\star$. Owl: $2\times250\ \mu\text{s}\times6400 = 3.2 \gg 1$ → place code. Human: $2\times700\ \mu\text{s}\times1600 = 2.2$ → marginal, which is exactly why human ITD coding is the most contested case in the literature. Gerbil: $2\times130\ \mu\text{s}\times1200 = 0.31 \ll 1$ → two-channel.

The clean conclusion, worth stating as a principle: **jitter sets the frequency; head size sets the code.** Two independent parameters, two independent design decisions, and the Jeffress-versus-slope debate concerned only the second.
</details>

**E6 (★★★).** Jeffress-style cross-correlation is a *linear* operation on the two ears' rate functions. Show that an optimal (Bayesian) ITD estimator is generally *not* the correlation peak, and identify the conditions under which they coincide.

<details markdown="1"><summary>Solution</summary>

Set up the estimation problem. Observations are the two spike trains $N_L, N_R$; the parameter is $\Delta$. Under conditionally-Poisson afferents with drive $\lambda(t)$ and $\lambda(t-\Delta)$, the log-likelihood is

$$\ell(\Delta) = \int \log\lambda(u)\,dN_L(u) + \int\log\lambda(v-\Delta)\,dN_R(v) - \int\lambda - \int\lambda(\cdot-\Delta).$$

The last two terms are $\Delta$-independent for a long stationary stimulus, so the ML estimate maximises $\int\log\lambda(v-\Delta)dN_R(v)$ — a correlation of the right ear's spikes with the *log* of the drive, not with the left ear's spikes.

When do these coincide? If $\lambda$ is unknown, the natural estimator substitutes the left-ear spike train (or its smoothed version) for $\lambda$. Then ML becomes $\max_\Delta \int \log\hat\lambda_L(v-\Delta)\,dN_R(v)$, and this reduces to plain cross-correlation $\max_\Delta \int \hat\lambda_L(v-\Delta)dN_R(v)$ only when $\log \hat\lambda_L \approx \text{const} + c\,\hat\lambda_L$, i.e. **when the modulation depth of the drive is small**. For strongly modulated (high vector strength) inputs the log compresses peaks, and ML weights the *troughs* of the drive more heavily than correlation does — because a right-ear spike arriving where the left drive is near zero is very strong evidence against that $\Delta$.

Three practical consequences:
1. Cross-correlation is the small-modulation limit of ML. At $\mathrm{VS}=0.95$ (§3) we are far from that limit, and a log-compressive nonlinearity before coincidence detection should improve performance measurably. Try it in the §3 code: replace the fibre counts by $\log(1+c\cdot\text{count})$ before correlating and compare peak sharpness.
2. Real coincidence detectors are *not* linear multipliers; the dendritic saturation of E3 is compressive, and the Kv1 differentiation is expansive on onsets. Whether the net operation is closer to correlation or to ML is an empirically open question and, to my knowledge, has not been asked in these terms.
3. Under prior information — most sounds are in front, head-shadow gives a correlated ILD cue — the posterior mean differs from the likelihood peak by a prior-dependent shift toward the midline. The Peña–Konishi multiplication of ITD and ILD tuning is precisely a factorised likelihood product, so the owl's ICx is doing Bayesian cue combination, not correlation-peak-picking. That reframing makes the space map a *posterior* map, which is a strictly stronger and more interesting level-2 claim than the Jeffress model, and it is the direction I would take this system if I worked on it.
</details>

## 5. Reading path

- **Jeffress (1948)**, *A place theory of sound localization* (J. Comp. Physiol. Psychol.) — read it for: three pages of pure algorithmic reasoning with no data, and a reminder of what a level-2 proposal looks like when made honestly.
- **Knudsen & Konishi (1978)**, *A neural map of auditory space in the owl* (Science) — read it for: the map, and for the standard of evidence needed to claim one.
- **Moiseff & Konishi (1981)**, *Neuronal and behavioral sensitivity to binaural time differences in the owl* (J. Neurosci.) — read it for: the establishment of ITD as the cue and the measurement of behavioural acuity that §2.4 has to explain.
- **Carr & Konishi (1988)**, *Axonal delay lines for time measurement in the owl's brainstem* (PNAS), and **Carr & Konishi (1990)**, *A circuit for detection of interaural time differences in the brain stem of the barn owl* (J. Neurosci.) — read them for: the anatomy that matches Jeffress, and for the recordings from the delay-line axons themselves.
- **Agmon-Snir, Carr & Rinzel (1998)**, *The role of dendrites in auditory coincidence detection* (Nature) — read it for: morphology as algorithm, and for E3.
- **Konishi (2003)**, *Coding of auditory space* (Annu. Rev. Neurosci.) — read it for: the mature synthesis of the owl system, including cross-frequency convergence and the resolution of phase ambiguity.
- **Peña & Konishi (2001)**, *Auditory spatial receptive fields created by multiplication* (Science) — read it for: multiplicative cue combination, and for E6's Bayesian reframing.
- **McAlpine, Jiang & Palmer (2001)**, *A neural code for low-frequency sound localization in mammals* (Nat. Neurosci.) — read it for: the observation that broke the universality of the Jeffress model. Read it slowly; the argument is subtler than its slogan.
- **Brand, Behrend, Marquardt, McAlpine & Grothe (2002)**, *Precise inhibition is essential for microsecond interaural time difference coding* (Nature) — read it for: a mechanism for internal delay with no counterpart in the Jeffress circuit.
- **Harper & McAlpine (2004)**, *Optimal neural population coding of an auditory spatial cue* (Nature) — read it for: the normative unification. This is the single most important paper in this unit, because it converts an apparent species disagreement into a parameter regime.
- **Grothe (2003)**, *New roles for synaptic inhibition in sound localization* (Nat. Rev. Neurosci.) — read it for: the mammalian counter-story argued by its principal advocate.
- **Joris & Yin (2007)**, *A matter of time: internal delays in binaural processing* (Trends Neurosci.) — read it for: a careful taxonomy of what an "internal delay" can physically be, which is where most confusion in this literature lives.
- **Grothe, Pecka & McAlpine (2010)**, *Mechanisms of sound localization in mammals* (Physiol. Rev.) — read it for: the comprehensive mammalian account, and for the evolutionary argument about why the two lineages diverged.
- **Ashida & Carr (2011)**, *Sound localization: Jeffress and beyond* (Curr. Opin. Neurobiol.) — read it for: the even-handed summary; read it last, after you have formed your own view.

## 6. Open problems and what would settle them

**Are the mammalian and avian systems different algorithms or different implementations of one?** Harper & McAlpine's answer — same objective, different optimum — is elegant, but it has been tested mainly by comparing across species with many confounded differences. *Settling it:* the decisive test is *within* an animal, across CF. In a species whose ITD range spans the transition (a large mammal, or human psychophysics), the theory predicts a *crossover*: low-CF channels should show slope coding and high-CF channels (if phase-locked at all) should show something more Jeffress-like, in the same brain, with the transition at $2\Delta_{\max}f \sim 1$. Nobody has looked for the crossover; a species where it should occur mid-range is the right preparation.

**What actually generates internal delay in mammals?** Axonal length, cochlear disparity, synaptic dynamics, and inhibitory timing are all supported by some evidence and none is established as dominant. *Settling it:* the mechanisms make different predictions under manipulation — cochlear disparity is frequency-locked and immutable, inhibitory timing is pharmacologically removable (Brand et al. did this), and axonal delay is temperature-dependent with a specific $Q_{10}$. A combined pharmacology-plus-temperature experiment on the same MSO neurons would separate them.

**Is the readout correlation or likelihood?** §E6 shows that cross-correlation is the small-modulation approximation to optimal estimation, and real inputs are strongly modulated. *Settling it:* measure NL/MSO input–output functions with binaural stimuli of controlled modulation depth and ask whether the effective operation is a product (correlation) or a product-of-logs (likelihood). The prediction is a systematic, modulation-dependent shift in the neuron's best delay that pure correlation cannot produce.

**The methodological question, which is the real content of this unit.** Jeffress's model was proposed from first principles, confirmed in exquisite anatomical detail, and taught as a universal for fifty years — and it is a special case. What generalises is not the delay line but the objective (maximise Fisher information about the cue) and the primitive (coincidence detection). The transferable discipline: when you extract an algorithm from a circuit, immediately ask *which parameters of the animal the extraction depended on*, and then predict what the same objective would build in an animal with different parameters. If you cannot answer that, you have described one circuit, not extracted an algorithm. Apply this test to C1's LSH (what if $N$ or the noise level were different?) and to C2's four readings (what if the oscillation frequency or the plume statistics were different?) — in both cases it generates experiments.
