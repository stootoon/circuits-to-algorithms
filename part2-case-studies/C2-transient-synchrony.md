---
title: "C2 · Transient synchrony"
parent: Case studies
nav_order: 2
---

# C2 — Transient synchrony and trajectories → reservoir computation and dynamic decorrelation

> **Circuit.** The locust antennal lobe: ~830 cholinergic projection neurons and ~300 non-spiking GABAergic local neurons in a densely reciprocal loop, driven by receptor neurons, producing a 20–30 Hz field oscillation, cycle-by-cycle transient synchrony of shifting PN subsets, and odour-specific spatiotemporal trajectories lasting a second or more.
> **Primitive extracted.** Contested. Four candidates: a fixed nonlinear dynamical expansion read out linearly (reservoir); a progressive whitening transform implemented by the transient itself; a stable heteroclinic sequence in which the trajectory *is* the code; and a synchrony-based gating signal that makes downstream coincidence detection selective. They are not the same claim and they are not equally supported.
> **Status.** Competing readings, unresolved — and this is the most instructive unresolved case in sensory neuroscience, because all four readings are consistent with the standard measurements and are discriminated only by analyses nobody has yet run on the right data.
> **Structures thread.** `structures/S-03-combinatorial-threshold-linear-networks.md` is the direct formal partner: CTLNs are the cleanest existing framework in which specific inhibitory connectivity graphs *provably* generate sequential attractors, which is hypothesis (c) with the mathematics done properly. `structures/S-04-manifold-capacity.md` gives the tools to ask whether the transient increases separability in the sense that matters. `structures/S-06-hyperbolic-odor-space.md` is what the decorrelation hypothesis is implicitly assuming about the target geometry. See `structures/README.md`.

## 1. The phenomenon

Six measurements. Learn them as measurements, separately from any story about them.

**1. Odour-evoked oscillation.** Odour presentation induces a 20–30 Hz oscillation in the locust mushroom body local field potential (Laurent & Davidowitz, 1994). It is not a property of individual neurons; it is a network state.

**2. Transient synchrony with shifting membership.** Individual PNs phase-lock to the oscillation during some cycles and not others. The *identity* of the synchronised assembly changes from cycle to cycle in an odour-specific, reproducible way (Laurent, Wehr & Davidowitz, 1996; Wehr & Laurent, 1996). A PN's spike is therefore labelled twice: by which neuron fired, and by which cycle it fired in.

**3. Two dissociable codes.** Picrotoxin blocks fast GABA_A-mediated inhibition in the AL. MacLeod & Laurent (1996) showed this abolishes the fast synchronisation while leaving the *slow* odour-specific temporal patterning of PN firing rates essentially intact. This dissociation is the single most useful fact in this unit, because it means synchrony and slow patterning can be manipulated separately, and therefore tested separately.

**4. Synchrony matters behaviourally, and selectively.** Stopfer, Bhagavan, Smith & Laurent (1997) desynchronised the AL with picrotoxin and tested honeybees behaviourally: discrimination of *molecularly similar* odorants was impaired; discrimination of dissimilar odorants was not. Rates were preserved; the fine-discrimination deficit came from the loss of temporal structure.

**5. Trajectories, not fixed points.** Mazor & Laurent (2005) recorded PN populations during 1-s odour pulses and analysed the population vector as a trajectory. The response has three phases: a fast transient after onset, an approach to a quasi-fixed point during sustained stimulation, and a distinct off-transient. A classifier trained on the population state performs *best during the transients* and worse at the fixed point — odour identity is carried by the trajectory, not by the steady state. This is the observation that anchors every reading below.

**6. Progressive decorrelation.** Representations of similar odours start close and diverge over the transient. In locust this appears as concentration and identity being progressively separated over the response (Stopfer, Jayaraman & Laurent, 2003) and as sequence-dependent trajectory separation (Broome, Jayaraman & Laurent, 2006). The cleanest version of the phenomenon is in zebrafish olfactory bulb: Friedrich & Laurent (2001) showed mitral cell ensemble patterns for chemically similar amino acids are highly correlated at odour onset and decorrelate over ~800 ms, while total activity is roughly conserved.

And downstream: Kenyon cells are sparse, oscillation-locked coincidence detectors (Perez-Orive et al., 2002), with intrinsic and circuit properties that favour coincidence detection over integration (Perez-Orive, Bazhenov & Laurent, 2004), and KC→β-lobe synapses show spike-timing-dependent plasticity referenced to the oscillation cycle (Cassenaer & Laurent, 2007).

## 2. The conversion

### 2.0 Stating the level-2 question properly

Everyone agrees on the level-1 description: a driven, recurrent, inhibition-dominated network with a fast oscillatory instability and slow adaptation, producing a stimulus-specific trajectory. The disagreement is about which *feature* of that trajectory is the computation, and hence what would count as damage.

Set up notation once. Let $s \in \mathbb{R}^{N_R}$ be the receptor drive, $x(t) \in \mathbb{R}^{N_P}$ the PN state. The network is a filter $x = \mathcal{F}[s]$ and the downstream reader applies some functional $g$. The four hypotheses are four claims about the joint design of $(\mathcal{F}, g)$:

| | what is load-bearing | what $\mathcal{F}$ is optimised for | falsified by |
|---|---|---|---|
| (a) reservoir | high dimension + fading memory | nothing task-specific | AL outperforming a matched random network on odour tasks specifically |
| (b) whitening | the *change* in the covariance spectrum over time | input statistics | spectrum not flattening; decorrelation not tracking natural odour statistics |
| (c) heteroclinic | the *order and timing* of metastable states | reproducible sequence generation | timing jitter that does not accumulate; no metastability |
| (d) gating | *which* PNs are coincident within a cycle | making a coincidence detector selective | desynchronisation without rate change leaving downstream selectivity intact |

### 2.1 Reading (a): reservoir / liquid-state computing

**The claim.** The AL is a fixed, non-adapted nonlinear dynamical system whose only job is to lift a low-dimensional, slowly-varying input into a high-dimensional trajectory with fading memory; all task specificity lives in a memoryless readout, which for the mushroom body is the KC layer plus KC→MBON weights.

**The mathematics.** Maass, Natschläger & Markram (2002) prove a universality theorem. Let $\mathbb{B}$ be a class of filters (the "liquid") and $\mathbb{F}$ a class of memoryless readouts. Say $\mathbb{B}$ has the *pointwise separation property* if for any two input histories $u \ne v$ there is some $B \in \mathbb{B}$ with $(Bu)(0) \ne (Bv)(0)$, and $\mathbb{F}$ has the *approximation property* if it is universal on compacts. Then the class of filters realisable as $\{g(B_1u, \dots, B_ku)\}$ is dense in the class of time-invariant filters with **fading memory**:

$$
\forall \varepsilon>0\ \exists \delta, T:\quad \sup_{t\in[-T,0]}\|u(t)-v(t)\| < \delta \;\Longrightarrow\; |F[u](0) - F[v](0)| < \varepsilon .
$$

By Boyd & Chua, fading-memory time-invariant filters are exactly those approximable by finite Volterra series. So: a generic recurrent network plus a trained linear readout is a universal approximator for the class of computations a nervous system could plausibly want.

The relevant quantitative object is the **separation ratio**. For two odours,

$$
\mathrm{Sep}(t) = \frac{\|x_{u}(t) - x_{v}(t)\|}{\|x_u(t)\| + \|x_v(t)\|},
$$

and the design tension is that a network which separates inputs strongly also amplifies noise, so a network at the edge of chaos (largest Lyapunov exponent $\lambda_{\max} \approx 0$) trades these off optimally (Legenstein & Maass, 2007).

**Why I am unimpressed as stated.** The theorem is an existence result about a *class* of networks; it says almost nothing about any *particular* network. "The AL is a reservoir" predicts: high dimension (true of essentially any recurrent circuit), fading memory (true unless there are attractors), linear decodability (true of most neural data). It is compatible with every measurement in §1 and forbids almost nothing. As usually deployed in the olfaction literature it is not a hypothesis, it is a reassurance.

It becomes a real hypothesis with one added word: *generic*. The sharp version is **"the AL is no better than a random network with matched size, timescales, and connection statistics, at any olfactory task."** That is falsifiable, and §2.5 says how.

### 2.2 Reading (b): progressive decorrelation as a normative objective

**The claim.** The transient is not incidental to representation-building; it *is* the algorithm. The AL performs a whitening transform, and it performs it *progressively*, so the amount of decorrelation applied is a function of elapsed time.

**The derivation, which is the best thing in this unit.** Take the linearised recurrent circuit with lateral inhibition $W \succeq 0$:

$$
\tau \dot r = -r - W r + s \quad\Longrightarrow\quad r(t) = \underbrace{(I+W)^{-1}\Big(I - e^{-(I+W)t/\tau}\Big)}_{\textstyle M(t)}\, s .
$$

Suppose $W$ shares eigenvectors with the input covariance $C_s$ — which is exactly what anti-Hebbian plasticity on the inhibitory weights, $\Delta W_{ij} \propto r_i r_j$ (Földiák, 1990; Atick & Redlich, 1990), delivers at convergence. In the shared eigenbasis, with $C_s$ eigenvalue $\lambda_k$ and $W$ eigenvalue $w_k$:

$$
m_k(t) = \frac{1 - e^{-(1+w_k)t/\tau}}{1+w_k}, \qquad \big[C_r(t)\big]_k = \lambda_k\, m_k(t)^2 .
$$

Read off the two limits.

*Early*, $t \ll \tau/(1+w_{\max})$: $m_k(t) \approx t/\tau$ for **every** $k$, so $C_r(t) \approx (t/\tau)^2 C_s$. The output covariance is a scaled copy of the input covariance — **no decorrelation at all at odour onset**.

*Late*, $t \gg \tau/(1+w_{\min})$: $m_k \to 1/(1+w_k)$, so $[C_r]_k = \lambda_k/(1+w_k)^2$. Choosing $1+w_k = \sqrt{\lambda_k}/c$ — i.e. $W = c^{-1}C_s^{1/2} - I$ — gives $[C_r]_k = c^2$ for all $k$: **exact whitening**, $M(\infty) = c\,C_s^{-1/2}$.

Now the punchline. The convergence rate of mode $k$ is $(1+w_k)/\tau = \sqrt{\lambda_k}/(c\tau)$. **High-variance modes converge fastest.** So as $t$ increases, the modes are suppressed in decreasing order of their input variance, and the output spectrum flattens *monotonically*:

$$
\boxed{\ \frac{d}{dt}\left[\frac{\lambda_{\max}m_{\max}^2}{\lambda_{\min}m_{\min}^2}\right] < 0 \quad\text{— the covariance spectrum flattens over the transient.}\ }
$$

Time is the whitening knob. The system is not "settling into" its representation; the representation *is* the settling, and where you sample the trajectory determines how much redundancy has been removed. This immediately explains Friedrich & Laurent's zebrafish result — patterns start correlated and decorrelate over hundreds of milliseconds — without any extra assumptions, and it explains why the *transient* is more informative than the fixed point (Mazor & Laurent) whenever the fixed point has been driven into adaptation-limited compression.

It also makes a prediction that is *directly measurable on existing datasets and has not been reported*: compute the eigenvalue spectrum of the PN population covariance across odours as a function of time from odour onset. Hypothesis (b) says the ratio $\lambda_1/\lambda_D$ must decrease monotonically, with the leading eigenvalue decaying first and fastest, and it says the transient duration should scale as $\tau c/\sqrt{\lambda_{\min}}$, i.e. odour panels with a wider covariance spectrum should require *longer* transients.

**Where the mechanism is not so simple.** Wiechert, Judkewitz, Riecke & Friedrich (2010) tested pattern decorrelation in the zebrafish bulb and found that simple broad lateral inhibition does *not* reproduce it; what works is sparse, specific recurrent connectivity that reshuffles activity between coactive channels. So the linear whitening derivation is the right *normative* statement and the wrong *mechanistic* one. That gap is where the interesting biology is.

### 2.3 Reading (c): winnerless competition and stable heteroclinic channels

**The claim.** The AL's connectivity implements a chain of saddle points connected by heteroclinic orbits. The trajectory is not "on its way" to a representation — the sequence of visited saddles, in order, with characteristic dwell times, *is* the representation. Different odours select different initial conditions and hence different channels through the same skeleton (Rabinovich et al., 2001; Afraimovich, Zhigulin & Rabinovich, 2004; Rabinovich, Huerta & Laurent, 2008).

**The mathematics.** Take generalised Lotka–Volterra dynamics for $m$ competing PN assemblies with activities $a_i \ge 0$:

$$
\dot a_i \;=\; a_i\Big(\sigma_i(S) - \sum_{j=1}^{m}\rho_{ij}a_j\Big) + \eta_i(t), \qquad \rho_{ii}=1,
$$

where $\sigma_i(S)$ is the odour-dependent drive to assembly $i$, $\rho_{ij} \ge 0$ is the inhibition of $i$ by $j$, and $\eta_i$ is noise of amplitude $\varepsilon$.

The axial fixed points are $A^{(k)}$: $a_k = \sigma_k$, all others zero. Linearise. Because $a_i = 0$ for $i\ne k$, the Jacobian is triangular and its diagonal gives the eigenvalues directly:

$$
\lambda^{(k)}_k = -\sigma_k \quad(\text{always stable, along the axis}), \qquad
\lambda^{(k)}_i = \sigma_i - \rho_{ik}\sigma_k \quad (i \ne k).
$$

For $A^{(k)}$ to be a saddle with a one-dimensional unstable manifold pointing at $A^{(k+1)}$ we need

$$
\boxed{\ \sigma_{k+1} > \rho_{k+1,k}\,\sigma_k \quad\text{and}\quad \sigma_i < \rho_{ik}\,\sigma_k \ \ \forall i \ne k, k+1 .\ }
$$

Read that in words: *the next assembly in the sequence must be the one the current assembly inhibits least; every other assembly must be inhibited enough to stay down.* The inhibition matrix $\rho$ must be strongly **asymmetric and specifically structured**. This is the crux of the whole unit. A random $\rho$ does not do this. Reading (c) makes a hard claim about connectivity that reading (a) explicitly denies.

**Stability of the channel.** Define the saddle value at $A^{(k)}$,

$$
\nu_k \;=\; \frac{\big|\lambda^{(k)}_{s,\text{lead}}\big|}{\lambda^{(k)}_{u}},
$$

the ratio of the least-negative stable eigenvalue to the unstable one. If $\nu_k > 1$ for every $k$ (all saddles dissipative), the heteroclinic sequence is an attracting set: trajectories starting in a finite-width tube around it are compressed toward it, and the object is a **stable heteroclinic channel**. If some $\nu_k < 1$, the trajectory leaks out and the sequence is not reproducible.

**Dwell times, and the prediction that discriminates.** Near saddle $k$, noise holds the unstable coordinate at a floor $\sim\varepsilon$; it must then grow to $O(1)$ at rate $\lambda_u^{(k)}$:

$$
T_k \;\approx\; \frac{1}{\lambda_u^{(k)}}\,\ln\frac{1}{\varepsilon} .
$$

Three consequences, each testable:

1. **Log sensitivity.** Increasing the noise by a factor of 10 shortens each dwell by only $\ln 10/\lambda_u$ — a fixed additive amount, not a proportional one. Sequence timing is robust but not rigid.
2. **Jitter accumulates.** Each escape time is a random variable with standard deviation $O(1/\lambda_u)$ that is *independent of $\varepsilon$ at leading order*. Escapes at successive saddles are approximately independent, so the time of the $k$-th transition has standard deviation growing like $\sqrt{k}/\lambda_u$. **Transition-time variance should grow linearly with transition index.**
3. **Order is rigid, geometry is not.** Perturbing the state transversally is corrected (that is what $\nu_k>1$ means); perturbing it along the channel is not corrected but only shifts phase.

Prediction 2 is the discriminator. In a *driven* system — a reservoir, or any circuit whose trajectory is clocked by the continuing stimulus — the timing of events is locked to the input and jitter does **not** accumulate: it stays roughly constant across the response. In an autonomous heteroclinic sequence, jitter compounds. Measuring the variance of PN transition times as a function of position in the response, with a stimulus whose own timing is precisely controlled, separates (a) from (c) cleanly. To my knowledge nobody has published this analysis on locust PN data, and the data to do it already exist.

**A caution about what has actually been shown.** The heteroclinic reading is supported by (i) beautiful mathematics, (ii) simulations of AL models that reproduce PN data (Bazhenov et al., 2001), and (iii) the qualitative observation of reproducible sequences. It is *not* supported by any direct measurement of $\rho_{ij}$ with the required asymmetry. That measurement is the missing link, and `structures/S-03-combinatorial-threshold-linear-networks.md` gives you the modern machinery for making it: CTLNs let you go from an inhibitory connectivity *graph* to a proof about which sequential attractors exist, which converts "measure $\rho$" into "measure the graph."

### 2.4 Reading (d): synchrony as gating for downstream coincidence detection

**The claim.** The oscillation is a clock, and transient synchrony is a *selection* signal: within each cycle, the subset of PNs that fire coincidentally is the subset that can drive a Kenyon cell. The AL is not computing the representation; it is formatting it for a coincidence detector.

**The mathematics.** Model a KC as counting inputs in a window $w$ (set by oscillation-locked feedforward inhibition, and by KC membrane and active properties — Perez-Orive et al., 2004). Let $K$ PNs converging on this KC be "synchronous," each firing once per cycle with timing jitter $\sigma$, and let $B$ other PNs fire as a Poisson process with rate $\nu$ over a cycle of period $T$. In a window of width $w$ centred on the packet,

$$
\mu_{\text{sig}} = K\,\mathrm{erf}\!\left(\frac{w}{2\sqrt2\,\sigma}\right), \qquad
\mu_{\text{bg}} = B\nu w, \qquad
\mathrm{Var} \approx B\nu w + K\,\mathcal{P}(1-\mathcal{P}),\ \ \mathcal{P}=\mathrm{erf}(\cdot).
$$

The discriminability of "packet present" against background is $d' = \mu_{\text{sig}}/\sqrt{\mathrm{Var}}$. For $w \gg \sigma$, $\mu_{\text{sig}}\to K$ while the noise grows as $\sqrt{w}$: shrinking $w$ helps. For $w \ll \sigma$, $\mu_{\text{sig}} \propto w/\sigma$ while noise still grows as $\sqrt w$: $d' \propto \sqrt{w}$, so shrinking $w$ hurts. The optimum is $w \sim \sigma$, giving

$$
d'_{\max} \;\approx\; \frac{K}{\sqrt{B\nu\sigma}} .
$$

Now desynchronise without touching rates — precisely what picrotoxin does (MacLeod & Laurent, 1996). The packet spreads over the whole cycle, $\sigma \to T$, and

$$
\frac{d'_{\text{sync}}}{d'_{\text{desync}}} \;=\; \sqrt{\frac{T}{\sigma}} .
$$

With $T = 50$ ms and $\sigma = 5$ ms this is a factor of $\sqrt{10} \approx 3.2$ loss in $d'$ — the difference between reliable and unreliable discrimination — obtained with **no change in any firing rate**. And it predicts the *selectivity* of the deficit: odour pairs that differ in which PNs fire are still discriminable from rates; odour pairs that differ only in which PNs are coincident are not. That is exactly the fine-versus-coarse pattern Stopfer et al. (1997) reported.

This reading closes a loop with C1 and C3. The KC's job (C1) is $\mathrm{WTA}(Ax)$; hypothesis (d) says the effective input $x$ is not the PN rate vector but the *per-cycle synchronous subset*, so the projection is re-sampled every 50 ms and the KC tag is a sequence of tags, not one tag. And the mechanism is the same coincidence-detection primitive that C3 develops in the barn owl, minus the delay lines.

Cassenaer & Laurent (2007) then complete the causal chain: STDP at KC→β-lobe synapses uses the oscillation as its temporal reference, so synchrony gates not only transmission but plasticity.

**My assessment.** (d) is the best-supported reading in the whole unit, because it is the only one with a complete causal chain — a manipulation that removes exactly one variable, a measured downstream consequence, and a behavioural consequence with the predicted selectivity. It is also the *narrowest*: it explains the oscillation and the synchrony, and says nothing about why the slow trajectory exists.

### 2.5 Adjudicating: how to turn each into a falsifiable level-2 claim

The discipline: for each hypothesis, write down (i) the quantity that must be measured, (ii) the value the hypothesis forbids, (iii) a manipulation that changes it.

**(a) Reservoir.** *Sharp form:* the AL is generic. *Measurement:* build a surrogate random network matched on $N_P$, $N_L$, mean connection probability, synaptic and adaptation time constants, and firing statistics. Train linear readouts on (i) odour identity, (ii) odour-concentration, (iii) a bank of random temporal functionals of the stimulus history (e.g. delayed XOR of two odour pulses). *Forbidden:* the real AL beating the surrogate on olfactory tasks by more than it beats it on generic tasks. *If the AL is specialised, (a) is demoted from hypothesis to null model* — which is still a useful role, and the one it should have had all along.

**(b) Whitening.** *Sharp form:* the covariance spectrum flattens monotonically over the transient, in the order predicted by $\sqrt{\lambda_k}$. *Measurement:* time-resolved eigenspectrum of the across-odour PN covariance; test monotonicity of $\lambda_1/\lambda_D$ and the ordering of decay rates. *Forbidden:* non-monotone flattening, or low-variance modes decaying first. *Manipulation:* change the input statistics (rear on a restricted odour environment; or use odour panels with engineered covariance) and ask whether the decorrelation *pattern* follows within the plasticity timescale. Whitening is a claim about matching input statistics, so if the transform does not track the statistics it is not whitening, it is just a fixed nonlinearity that happens to decorrelate.

**(c) Heteroclinic sequence.** *Sharp form:* metastability plus accumulating jitter plus specific asymmetric $\rho$. *Measurements:* (i) fit a hidden Markov model to simultaneously recorded PNs and test it against a smooth-flow (e.g. latent linear dynamical system) alternative by cross-validated likelihood — the approach Jones et al. (2007) took in gustatory cortex; (ii) transition-time variance versus transition index, which must be linear; (iii) dwell time versus injected noise amplitude, which must be $\propto \ln(1/\varepsilon)$ and not $\propto 1/\varepsilon$ or constant; (iv) reconstruct the LN→PN and LN→LN graph and test the saddle inequalities. *Forbidden:* constant jitter across the response; smooth-flow models winning on held-out likelihood.

**(d) Gating.** *Sharp form:* downstream selectivity depends on PN coincidence, not PN rate. *Measurement:* the modern version of Stopfer 1997 — optogenetically inject timing jitter into PNs while holding rates fixed (or drive PNs with rate-matched but temporally scrambled patterns) and read out KC responses and behaviour. *Forbidden:* KC sparseness and odour selectivity surviving rate-matched desynchronisation.

**Are they mutually exclusive?** No — and pretending otherwise is the standard error. A circuit can whiten its inputs *and* do so via a heteroclinic skeleton *and* format the result for a coincidence detector. What is at stake is which property is load-bearing: which one, if you removed it, would break the animal. My ranking, stated so you can disagree with it precisely:

1. **(d) is true and demonstrated** for the fast timescale, with the strongest causal chain in the literature.
2. **(b) is the best normative account** of the slow timescale and makes the sharpest untested prediction (spectrum flattening); its mechanistic implementation is probably not lateral inhibition in the textbook sense (Wiechert et al., 2010).
3. **(c) is the most specific dynamical account** and the most valuable *because* it is the most falsifiable; it stands or falls on the jitter-accumulation and connectivity measurements, neither of which has been done properly in locust.
4. **(a) is not a hypothesis in its usual form.** Promote it to a null model, match it carefully, and use it to measure how much of the AL's behaviour is generic. That is a real and useful job.

## 3. Worked example / model to build

A rate-based PN–LN network with slow GABA and spike-frequency adaptation, driven by eight overlapping odours. It gives you an odour-specific transient trajectory, progressive decorrelation, and a linear readout whose performance peaks in the transient and falls at the fixed point — the Mazor–Laurent signature — in ~35 lines.

```python
import numpy as np
NP, NL, dt = 60, 20, 1.0                  # PNs, LNs, timestep (ms)
tP, tL, tA = 20., 50., 400.               # PN, GABA_slow, adaptation (ms)
relu = lambda z: np.maximum(z, 0.)
rng  = np.random.default_rng(0)
spw  = lambda a,b,g,p=.5: (rng.random((a,b))<p)*rng.gamma(3.,1/3.,(a,b))*g/(b*p)
WLP, WPL = spw(NL,NP,3.0), spw(NP,NL,12.0)          # PN->LN , LN-|PN

def odor_bank(n=8, k=18, shared=14, seed=3):        # 8 overlapping odours
    r = np.random.default_rng(seed); S = np.zeros((n,NP))
    base = r.choice(NP,shared,replace=False); rest = np.setdiff1d(np.arange(NP),base)
    for i in range(n):
        idx = np.concatenate([base, r.choice(rest,k-shared,replace=False)])
        S[i,idx] = r.gamma(3.,1/3.,k)*2.0
    return S

def run(s, T=2000, on=200, off=1400, alpha=3.0):
    P=np.zeros(NP); L=np.zeros(NL); A=np.zeros(NP); R=np.zeros((T,NP))
    for t in range(T):
        d = s if on < t < off else 0.
        P += dt*(-P + relu(d - WPL@L - alpha*A))/tP
        L += dt*(-L + relu(WLP@P))/tL
        A += dt*(-A + P)/tA
        R[t] = P
    return R

S = odor_bank(); M = np.array([run(s) for s in S])   # (8, T, NP)
cc = lambda a,b: float((a-a.mean())@(b-b.mean())/
                       (np.linalg.norm(a-a.mean())*np.linalg.norm(b-b.mean())+1e-12))
pairs = [(i,j) for i in range(8) for j in range(i)]
noise, R_ = 0.6, np.random.default_rng(7)
print("input  <r> = %.3f" % np.mean([cc(S[i],S[j]) for i,j in pairs]))
print("  t(ms)   <r>_PN   |r|    8-way acc")
for t in [205,225,260,320,450,700,1000,1350,1450]:
    C = M[:,t]
    X = C[:,None,:] + noise*R_.standard_normal((8,300,NP))          # 300 noisy trials
    pred = np.argmin(((X[:,:,None,:]-C[None,None])**2).sum(-1),-1)  # nearest centroid
    print("%6d   %5.2f   %5.2f   %5.2f" % (t, np.mean([cc(C[i],C[j]) for i,j in pairs]),
          np.linalg.norm(C,axis=1).mean(), (pred==np.arange(8)[:,None]).mean()))
Z = M[:,200:1400].reshape(-1,NP); Z = Z - Z.mean(0)
ev = np.linalg.svd(Z, compute_uv=False)**2
print("PC var frac:", np.round((ev/ev.sum())[:5],3), " participation ratio %.1f" %
      (ev.sum()**2/(ev**2).sum()))
```

**What to look for.** Three things, and you should get them almost exactly.

*Progressive decorrelation.* Input odours have mean pairwise correlation $0.51$. The PN representation starts at $0.49$ at $t=205$ ms (5 ms after onset — no decorrelation yet, exactly as §2.2 predicts for $t \ll \tau$), falls to $0.30$ by 225 ms, and bottoms out around $0.16$ by 320 ms. It then creeps back up to $0.21$ at the fixed point. Correlations are cut by a factor of three, and — this is the point — *the amount of decorrelation is a function of when you look*.

*Transient beats fixed point.* Eight-way classification accuracy under fixed readout noise is $0.76$ at 205 ms, peaks at $0.97$ at 225 ms, and decays to $0.55$ at the fixed point (700–1350 ms). Compare the response norm column: the transient is both larger and less correlated, and both matter. This is a working model of Mazor & Laurent's central result.

*Dimensionality.* The participation ratio of the trajectory ensemble is $\approx 6$, with the top five PCs holding 26/21/16/12/10% — the trajectory is genuinely multidimensional, not a one-dimensional ramp.

**Things to do to it.**
1. **Make it oscillate.** As written it does not, and that is instructive: see E1. Add a 5 ms conduction delay in the LN→PN path and re-run.
2. **Test hypothesis (b) properly.** Compute the eigenspectrum of the across-odour covariance at each time point and check whether $\lambda_1/\lambda_D$ decreases monotonically and whether the top mode decays first. Then set $W_{PL}W_{LP} \propto C_s^{1/2}$ by construction and see how much better it gets.
3. **Test hypothesis (a).** Freeze the network and train a ridge readout to compute a *delayed* function of the odour sequence (was odour A presented two pulses ago?). Measure how memory depth scales with $\tau_A$. This is where the fading-memory theory earns its keep.
4. **Build (c) alongside it.** Implement the GLV system of §2.3 with $m=6$, $\sigma_i=1$, $\rho_{k+1,k}=0.5$, $\rho_{ik}=2$ otherwise, $\varepsilon=10^{-4}$, and verify the dwell-time law $T \propto \ln(1/\varepsilon)$ by sweeping $\varepsilon$ over four decades. Then ask what would distinguish its output from the PN–LN model's.

## 4. Exercises

**E1 (★★).** The model in §3 does not oscillate at 20–30 Hz. Prove that a two-population rate model with only reciprocal PN↔LN coupling and first-order dynamics *cannot* produce a growing oscillation, and compute the delay needed to place the instability at 25 Hz for $\tau_P=6$ ms, $\tau_L=12$ ms.

<details markdown="1"><summary>Solution</summary>

Linearise about an operating point where both populations are active. The loop transfer function for the negative feedback PN→LN→PN is
$$\mathcal{L}(i\omega) = \frac{-G}{(1+i\omega\tau_P)(1+i\omega\tau_L)},$$
with loop gain $G = g_{PL}g_{LP} > 0$. Instability (Nyquist) requires $\mathcal{L}(i\omega)=1$ for some real $\omega$, i.e. total phase lag of the *negative* loop equal to $180°$:
$$\arctan(\omega\tau_P) + \arctan(\omega\tau_L) = \pi .$$
Each arctangent is strictly less than $\pi/2$, so the sum is strictly less than $\pi$ for every finite $\omega$. The phase condition is never met: **two first-order lags cannot oscillate at any gain.** (Increasing $G$ only moves the two real poles further left along the negative real axis or makes them complex with negative real part; the system is unconditionally stable for a purely inhibitory loop with two poles.)

You need a third pole or a delay. With a pure delay $\Delta$,
$$\omega\Delta + \arctan(\omega\tau_P) + \arctan(\omega\tau_L) = \pi .$$
At $f = 25$ Hz, $\omega = 2\pi f = 0.157\ \mathrm{rad\,ms^{-1}}$. Then $\arctan(0.157\times 6)=\arctan(0.943)=0.756$ rad, $\arctan(0.157\times 12)=\arctan(1.885)=1.083$ rad. So
$$\omega\Delta = \pi - 0.756 - 1.083 = 1.303\ \mathrm{rad} \;\Longrightarrow\; \Delta = 1.303/0.157 = 8.3\ \mathrm{ms}.$$
Required gain at that frequency: $G \ge \sqrt{1+0.943^2}\,\sqrt{1+1.885^2} = 1.374\times 2.134 = 2.93$. So a ~8 ms effective delay and a loop gain above ~3 puts the instability at 25 Hz. Biologically that "delay" is not literal axonal conduction — it is the finite rise time of the GABA_A conductance plus synaptic and dendritic delay, which is exactly a third pole. Replacing the delay with a third lag $\tau_S$ and solving $\arctan(\omega\tau_P)+\arctan(\omega\tau_L)+\arctan(\omega\tau_S)=\pi$ at $\omega=0.157$ needs $\arctan(\omega\tau_S)=1.303$ rad, i.e. $\omega\tau_S = 3.67$, $\tau_S = 23$ ms — a slow GABA component, which is what the AL has.
</details>

**E2 (★★).** Derive the whitening result of §2.2 in the form actually needed for data analysis: show that if $W = c^{-1}C_s^{1/2} - I$, the time at which mode $k$ has completed a fraction $1-e^{-1}$ of its suppression is $t_k = c\tau/\sqrt{\lambda_k}$, and use it to predict the shape of the "decorrelation time course" that Friedrich & Laurent measured.

<details markdown="1"><summary>Solution</summary>

From $m_k(t) = (1-e^{-(1+w_k)t/\tau})/(1+w_k)$ with $1+w_k=\sqrt{\lambda_k}/c$, the exponential rate is $(1+w_k)/\tau = \sqrt{\lambda_k}/(c\tau)$, hence $t_k = c\tau/\sqrt{\lambda_k}$. High-variance modes have small $t_k$.

Now the observable. Friedrich & Laurent measured the mean pairwise correlation between odour representations, $\bar\rho(t)$. For zero-mean patterns, $\bar\rho(t)$ is a monotone decreasing function of the *spectral flatness* of $C_r(t)$: if $C_r$ has eigenvalues $\{\lambda_k m_k^2\}$ and the odour patterns are generic in that basis, then
$$\bar\rho(t) \;\approx\; \frac{\sum_k \lambda_k m_k(t)^2\,\beta_k}{\sum_k \lambda_k m_k(t)^2},$$
where $\beta_k \in [-1,1]$ is the mean pairwise alignment of the odour set in mode $k$; the dominant mode of a natural odour ensemble is a common "total activity / broadly-tuned" direction with $\beta_1 \approx 1$, and higher modes are near zero.

So to leading order $\bar\rho(t) \approx \lambda_1 m_1(t)^2 / \sum_k \lambda_k m_k(t)^2$. At $t\to 0$, all $m_k \propto t$, so $\bar\rho(0) = \lambda_1/\sum\lambda_k$ — the input correlation. As $t$ grows, $m_1$ saturates first (rate $\sqrt{\lambda_1}$ is largest) at the small value $c/\sqrt{\lambda_1}$, while the small-$\lambda$ modes are still growing linearly. Hence
$$\bar\rho(t) \approx \frac{c^2}{c^2 + \sum_{k\ge2}\lambda_k m_k(t)^2}$$
once the leading mode has saturated — a monotone sigmoidal decay from $\lambda_1/\sum\lambda_k$ to $1/D$, with the *fall time set by the second eigenvalue*, $t_2 = c\tau/\sqrt{\lambda_2}$, and the *completion time set by the smallest*, $t_D = c\tau/\sqrt{\lambda_D}$.

Prediction, matching the published shape: a sigmoidal decorrelation curve, not exponential; a plateau at $\bar\rho \approx 1/D$; and a duration that grows as the *dynamic range* of the input spectrum grows. The last of these is testable by comparing decorrelation time courses for odour panels of different spectral flatness, and I do not know of anyone who has done it.
</details>

**E3 (★★).** For the GLV system of §2.3 with $m$ assemblies, $\sigma_i \equiv \sigma$, $\rho_{k+1,k} = \rho_-$, and $\rho_{ik} = \rho_+$ otherwise (all indices mod $m$), find the conditions on $(\rho_-, \rho_+)$ for a stable heteroclinic cycle, and compute the period.

<details markdown="1"><summary>Solution</summary>

Eigenvalues at $A^{(k)}$: along the axis, $-\sigma$; for $i = k+1$, $\lambda_u = \sigma(1-\rho_-)$; for all other $i$, $\lambda_s = \sigma(1-\rho_+)$.

**Saddle condition.** Need $\lambda_u > 0 \Rightarrow \rho_- < 1$, and $\lambda_s < 0 \Rightarrow \rho_+ > 1$. So
$$\rho_- < 1 < \rho_+ .$$
Interpretation: assembly $k$ suppresses its successor *less* than it suppresses itself ($\rho_{kk}=1$), and suppresses everything else *more*. Note the successor also needs to be able to grow, which is the same inequality.

**Channel stability.** The leading stable eigenvalue is the least negative of $\{-\sigma, \sigma(1-\rho_+)\}$. Assume $\rho_+ < 2$ so that $\sigma(1-\rho_+) > -\sigma$ and the transverse direction is the leading stable one. Then
$$\nu = \frac{|\lambda_s|}{\lambda_u} = \frac{\rho_+ - 1}{1-\rho_-} \;>\; 1 \quad\Longleftrightarrow\quad \rho_+ + \rho_- > 2 .$$
Combined: $\rho_- < 1 < \rho_+$ **and** $\rho_+ > 2 - \rho_-$. E.g. $\rho_-=0.5$ requires $\rho_+>1.5$. For $\rho_+ \ge 2$ the axial direction becomes the leading stable one and $\nu = \sigma/[\sigma(1-\rho_-)] = 1/(1-\rho_-) > 1$ automatically, so the cycle is stable for any $\rho_-<1$.

**Period.** Each dwell is $T_k = \lambda_u^{-1}\ln(1/\varepsilon) = \big[\sigma(1-\rho_-)\big]^{-1}\ln(1/\varepsilon)$, and there are $m$ of them:
$$T_{\text{cycle}} \;=\; \frac{m}{\sigma(1-\rho_-)}\,\ln\frac{1}{\varepsilon}.$$
Sanity check with locust numbers: if a "state" lasts ~100 ms and $\ln(1/\varepsilon)\approx 9$ (for $\varepsilon\sim10^{-4}$), then $\sigma(1-\rho_-) \approx 90\ \mathrm{s^{-1}}$, i.e. an unstable rate with time constant ~11 ms — comfortably in the range of AL synaptic dynamics. The model is not asking for anything unphysiological.

Note the dependence structure that makes this falsifiable: the period is proportional to $\ln(1/\varepsilon)$ and inversely proportional to $(1-\rho_-)$. Pharmacologically reducing inhibition raises $\rho_-$? No — it *lowers* all $\rho$, which raises $\lambda_u$ and *shortens* dwells while simultaneously threatening $\nu>1$. So partial disinhibition should speed the sequence up and then abruptly destabilise its order. That is a distinctive, testable signature that neither (a) nor (b) predicts.
</details>

**E4 (★★).** Work out the jitter-accumulation prediction quantitatively. If the escape time from each saddle is $T_k = \lambda_u^{-1}\ln(\xi_k^{-1})$ with $\xi_k$ the noise-set initial condition drawn as $|\mathcal{N}(0,\varepsilon^2)|$, compute the mean and variance of $T_k$ and of the time of the $n$-th transition.

<details markdown="1"><summary>Solution</summary>

Let $\xi = \varepsilon|Z|$, $Z\sim\mathcal N(0,1)$. Then $T = \lambda_u^{-1}[\ln(1/\varepsilon) - \ln|Z|]$.

Mean: $\mathbb{E}[\ln|Z|] = \tfrac12(\psi(1/2)+\ln 2) = \tfrac12(-\gamma - 2\ln 2 + \ln 2) = -\tfrac12(\gamma+\ln 2) \approx -0.635$. So
$$\mathbb{E}[T] = \frac{1}{\lambda_u}\left[\ln\frac{1}{\varepsilon} + 0.635\right].$$

Variance: $\mathrm{Var}(\ln|Z|) = \tfrac14\psi'(1/2) = \tfrac14\cdot\frac{\pi^2}{2} = \frac{\pi^2}{8} \approx 1.234$. So
$$\mathrm{Var}(T) = \frac{\pi^2}{8\lambda_u^2} \approx \frac{1.234}{\lambda_u^2},\qquad \mathrm{SD}(T) \approx \frac{1.11}{\lambda_u},$$
**independent of $\varepsilon$.** The mean depends on the noise only logarithmically; the jitter does not depend on it at all.

For $n$ successive independent escapes, the time of the $n$-th transition has
$$\mathbb{E}[t_n] = \frac{n}{\lambda_u}\left[\ln\frac1\varepsilon + 0.635\right], \qquad \mathrm{Var}(t_n) = n\,\frac{\pi^2}{8\lambda_u^2},$$
so $\mathrm{SD}(t_n) = 1.11\sqrt{n}/\lambda_u$ and the **coefficient of variation falls as $1/\sqrt n$** while the absolute jitter grows as $\sqrt n$.

The experiment: align many trials to odour onset, identify state transitions (HMM), and plot $\mathrm{Var}(t_n)$ against $n$. Heteroclinic: strictly linear through the origin with slope $\pi^2/(8\lambda_u^2)$, and the slope must be *invariant* to changing the noise level. Stimulus-clocked (reservoir, or any input-driven trajectory): $\mathrm{Var}(t_n)$ roughly constant in $n$, set by stimulus and sensory jitter. Limit cycle with a phase-resetting input: also roughly constant. This single plot separates autonomous from driven sequence generation, and with locust PN data recorded during long, precisely-timed odour pulses it is a few hours of analysis.
</details>

**E5 (★★, computational).** Take the §3 model, add a 5 ms delay in the LN→PN pathway to make it oscillate, and then implement the hypothesis-(d) test: measure, cycle by cycle, which PNs spike within a $\pm5$ ms window of the population peak, and ask whether an 8-way classifier on the *synchronous subset* outperforms one on the full rate vector.

<details markdown="1"><summary>Solution</summary>

Implementation. Add `buf = np.zeros((5, NL))` and replace `WPL@L` with `WPL@buf[t % 5]`, storing `buf[t % 5] = L` after use; raise the loop gain until the population mean shows a >8 Hz peak in its power spectrum (E1 says you need $G \gtrsim 3$ at the crossover). Define the LFP as $-\sum_j L_j(t)$, find its peaks, and for each cycle form the binary vector $c_i = \mathbb 1[P_i(t) > \theta$ within $\pm 5$ ms of the peak$]$ with $\theta$ chosen so ~20% of PNs qualify.

Expected result and why. Two classifiers: (i) nearest-centroid on the continuous rate vector averaged over the cycle; (ii) nearest-centroid on the binary coincidence vector $c$. With *no readout noise*, (i) wins — it has strictly more information ($c$ is a thresholded function of the rates). With readout noise applied downstream (i.e. at the KC), (ii) wins above a noise threshold, because the binarisation is a matched nonlinearity for a coincidence detector: it discards amplitude, which carries little odour-identity information here, and keeps the ordinal pattern, which carries most of it.

The instructive part is the crossover. Sweep the readout noise $\sigma$ and find where the two curves cross. That crossover is the quantitative statement of hypothesis (d): *synchrony-based readout is superior precisely in the regime where the downstream reader is noisy and integrates over one cycle*. If you make the KC integrate over 200 ms instead of 20 ms, the advantage of $c$ disappears entirely — which is why Perez-Orive, Bazhenov & Laurent (2004) had to show that KC intrinsic properties *prevent* long integration. That paper is the hidden load-bearing member of the whole synchrony argument, and this simulation shows you why.
</details>

**E6 (★★★).** Design the experiment that would demote reading (a) from hypothesis to null model. Specify the surrogate network, the task bank, and the statistic, and state what result would *support* (a) rather than merely fail to reject it.

<details markdown="1"><summary>Solution</summary>

**Surrogate.** Construct a randomised control that preserves everything except putative task-specific structure. Three nested surrogates, in increasing strictness:
- *S1 (degree-preserving rewire):* shuffle the LN↔PN graph preserving in/out degrees and weight distributions. Kills any specific $\rho$ structure (§2.3) but preserves gross statistics.
- *S2 (spectrum-matched):* additionally match the eigenvalue spectrum of the linearised Jacobian, so the surrogate has the same intrinsic timescales and the same distance from the edge of chaos.
- *S3 (response-matched):* additionally match single-PN response statistics (rates, adaptation, oscillation phase-locking) by re-tuning intrinsic parameters.

**Task bank.** Two families, matched in difficulty:
- *Olfactory:* identity of 16 monomolecular odorants; identity in the presence of a background (Saha et al.-style); concentration-invariant identity; identity of a binary mixture.
- *Generic:* delayed match-to-sample on arbitrary pulse patterns; $n$-back on a pulse train; parity of the number of pulses in the last 500 ms; reconstruction of a random band-limited input from $t-\Delta$.

**Statistic.** For each network $\mathcal{N}\in\{\text{real}, S1, S2, S3\}$ and task $T$, train a ridge readout on the population state and record cross-validated performance $A(\mathcal N, T)$. Compute the *specialisation index*
$$\Sigma = \big[\bar A(\text{real}, \text{olf}) - \bar A(S, \text{olf})\big] - \big[\bar A(\text{real}, \text{gen}) - \bar A(S, \text{gen})\big],$$
i.e. the interaction term. $\Sigma \approx 0$ means whatever advantage the real network has is task-general — it is a better reservoir, not a specialised one. $\Sigma > 0$ means the AL is tuned to olfaction, and (a) in its strong form is false.

**What would support (a).** Not $\Sigma \approx 0$ alone — that is a null result. Positive support requires two further things. (i) A *dose-response*: performance on the whole task bank should track a single scalar summary of the network's dynamics (e.g. the effective rank of the state covariance, or $\lambda_{\max}$), across the real network and all surrogates, with the real network sitting on the same curve. That is the signature of a generic reservoir: one knob, many tasks. (ii) *Edge-of-chaos placement*: the real network should sit near the maximum of that curve, and manipulations that move $\lambda_{\max}$ (e.g. graded disinhibition) should move performance along it in the predicted direction. If both hold, "the AL is a well-tuned generic reservoir" is a substantive, supported claim rather than a restatement of the data. That is the standard I would hold it to, and nobody has met it.
</details>

**E7 (★★★).** The four hypotheses are usually treated as rivals. Construct a single model in which all four are simultaneously true, and identify the one measurement that would still distinguish which is load-bearing.

<details markdown="1"><summary>Solution</summary>

**The unified construction.** Take the GLV/heteroclinic skeleton of §2.3 as the slow dynamics on assembly activities $a_k$, and add a fast oscillatory instability on top of it (E1's delayed inhibition), so the network state is
$$x(t) = \Big(\textstyle\sum_k a_k(t)\,v_k\Big)\cdot\big(1 + \eta\cos(\omega t + \varphi(t))\big),$$
with $v_k$ the assembly's PN pattern. Choose the competition matrix $\rho$ by an anti-Hebbian rule on the odour ensemble, $\Delta\rho_{ij}\propto \langle a_ia_j\rangle$, symmetrised only partially so the required asymmetry survives.

Then: (c) holds by construction — the slow dynamics are a stable heteroclinic channel. (b) holds because anti-Hebbian $\rho$ makes the sequence of visited saddles progressively orthogonalise the odour set: successive states are chosen to be maximally competitive, hence maximally decorrelated, so the ensemble covariance flattens over the transient exactly as §2.2 requires. (d) holds because the multiplicative oscillation makes the assembly currently near its saddle the one whose PNs are coincident within a cycle, so a downstream coincidence detector reads out $\arg\max_k a_k$ — the current state — and nothing else. (a) holds trivially: the whole thing is a fixed nonlinear filter with fading memory read out linearly, and the LSM theorem applies.

So all four are true *of this model*. The rivalry was never about truth values; it was about which property, if abolished, breaks the computation.

**The distinguishing measurement.** Perturb along each hypothesis's own axis, holding the others fixed, and measure downstream (KC) selectivity:
- Abolish only the oscillation (set $\eta=0$), leaving $a_k(t)$ untouched. If KC selectivity survives, (d) is not load-bearing.
- Abolish only the sequence order (make $\rho$ symmetric, so the system goes to a single winner-take-all fixed point) while preserving the *set* of active assemblies and their total decorrelation. If odour discrimination survives, (c) is not load-bearing: the code is the visited set, not the sequence.
- Abolish only the decorrelation (scale $\rho \to \rho/\kappa$ with $\kappa$ large, preserving asymmetry and hence sequence order, but flattening the competition so states overlap). If discrimination of *correlated* odour pairs specifically collapses while uncorrelated pairs survive, (b) is load-bearing.
- Randomise the graph (E6's S2) while preserving dimension and timescales. If everything survives, (a) is the whole story.

The point of the exercise: **a level-2 claim is only meaningful relative to a perturbation that isolates it.** "Which algorithm is the circuit running" is not well-posed; "which property of the circuit is necessary for which behaviour" is. Every one of the four perturbations above is achievable in a model in an afternoon, and at least the first is achievable in the animal. That asymmetry — easy in silico, hard in vivo — is why this debate has lasted twenty-five years.
</details>

## 5. Reading path

- **Laurent & Davidowitz (1994)**, *Encoding of olfactory information with oscillating neural assemblies* (Science) — read it for: the founding observation; note how much of the eventual theory is already visible in the framing.
- **Laurent, Wehr & Davidowitz (1996)**, *Temporal representations of odors in an olfactory network* (J. Neurosci.) and **Wehr & Laurent (1996)**, *Odour encoding by temporal sequences of firing in oscillating neural assemblies* (Nature) — read them for: the demonstration that assembly membership changes cycle by cycle, which is the fact all four hypotheses must explain.
- **MacLeod & Laurent (1996)**, *Distinct mechanisms for synchronization and temporal patterning of odor-encoding neural assemblies* (Science) — read it for: the dissociation that makes the whole field experimentally tractable. Read it *before* Stopfer 1997, because it tells you what picrotoxin does and does not do.
- **Stopfer, Bhagavan, Smith & Laurent (1997)**, *Impaired odour discrimination on desynchronization of odour-encoding neural assemblies* (Nature) — read it for: the causal link from synchrony to behaviour, and for the selectivity of the deficit, which is the part usually mis-remembered.
- **Perez-Orive, Mazor, Turner, Cassenaer, Wilson & Laurent (2002)**, *Oscillations and sparsening of odor representations in the mushroom body* (Science) — read it for: the downstream consequence, and the bridge to C1.
- **Perez-Orive, Bazhenov & Laurent (2004)**, *Intrinsic and circuit properties favor coincidence detection for decoding oscillatory input* (J. Neurosci.) — read it for: why the KC cannot simply integrate, which is the hidden premise of hypothesis (d).
- **Mazor & Laurent (2005)**, *Transient dynamics versus fixed points in odor representations by locust antennal lobe projection neurons* (Neuron) — read it for: the trajectory framing and the transient-beats-fixed-point result. This is the paper the unit is built around.
- **Friedrich & Laurent (2001)**, *Dynamic optimization of odor representations by slow temporal patterning of mitral cell activity* (Science) — read it for: the cleanest measurement of progressive decorrelation anywhere, and the normative framing of hypothesis (b).
- **Wiechert, Judkewitz, Riecke & Friedrich (2010)**, *Mechanisms of pattern decorrelation by recurrent neuronal circuits* (Nat. Neurosci.) — read it for: the demonstration that the obvious mechanism (broad lateral inhibition) is the wrong one. This is a model of how to test a mechanism against a normative claim.
- **Rabinovich, Volkovskii, Lecanda, Huerta, Abarbanel & Laurent (2001)**, *Dynamical encoding by networks of competing neuron groups: winnerless competition* (Phys. Rev. Lett.) — read it for: the original WLC proposal in the AL context.
- **Afraimovich, Zhigulin & Rabinovich (2004)**, *On the origin of reproducible sequential activity in neural circuits* (Chaos) — read it for: the saddle-value condition and the actual theorems behind "stable heteroclinic channel."
- **Rabinovich, Huerta & Laurent (2008)**, *Transient dynamics for neural processing* (Science) — read it for: the mature statement of the transient-as-code position, and as the target for E4.
- **Bazhenov, Stopfer, Rabinovich, Huerta, Abarbanel, Sejnowski & Laurent (2001)**, the pair of *Neuron* papers on modelling transient oscillatory synchronisation and odour-evoked temporal patterning in the locust AL — read them for: the biophysical model that connects the abstractions above to conductances, and for how the two timescales are separately generated.
- **Maass, Natschläger & Markram (2002)**, *Real-time computing without stable states* (Neural Computation) — read it for: the LSM theorem, stated precisely. Then reread §2.1 and decide whether the theorem constrains anything about the locust.
- **Buonomano & Maass (2009)**, *State-dependent computations: spatiotemporal processing in cortical networks* (Nat. Rev. Neurosci.) — read it for: the broad case for reading (a) across systems, argued well by its strongest advocates.
- **Legenstein & Maass (2007)**, *Edge of chaos and prediction of computational performance for neural circuit models* (Neural Networks) — read it for: the one quantitative handle reservoir theory offers, which E6 leans on.
- **Stopfer, Jayaraman & Laurent (2003)**, *Intensity versus identity coding in an olfactory system* (Neuron) and **Broome, Jayaraman & Laurent (2006)**, *Encoding and decoding of overlapping odor sequences* (Neuron) — read them for: how the trajectory picture handles concentration and stimulus history, and for the strongest evidence that the AL's state is genuinely dynamic rather than a delayed rate code.
- **Cassenaer & Laurent (2007)**, *Hebbian STDP in mushroom bodies facilitates the synchronous flow of olfactory information* (Nature) — read it for: the plasticity rule that uses the oscillation as its clock, closing the loop on hypothesis (d).
- **Jones, Fontanini, Sadacca, Miller & Katz (2007)**, *Natural stimuli evoke dynamic sequences of states in sensory cortical ensembles* (PNAS) — read it for: the HMM methodology you need to run the metastability test in §2.5, developed on a different system.
- **Földiák (1990)**, *Forming sparse representations by local anti-Hebbian learning* (Biol. Cybern.) — read it for: the local learning rule that makes §2.2's whitening biologically plausible.

## 6. Open problems and what would settle them

**Does timing jitter accumulate?** This is the single highest-value unanswered question in the unit, and the cheapest to answer. §2.3/E4 give the prediction: for an autonomous heteroclinic sequence, $\mathrm{Var}(t_n) = n\pi^2/(8\lambda_u^2)$, linear in transition index with a noise-independent slope; for a stimulus-clocked trajectory, roughly constant. *Settling it:* HMM-segment simultaneously recorded locust PN populations during long, precisely-timed odour pulses; regress transition-time variance on index; repeat with elevated input noise. Data of the required type already exist in several labs.

**Is the LN→PN graph structured as WLC requires?** The saddle inequalities of §2.3 are a *hard, quantitative* prediction about connectivity: for each state $k$, exactly one assembly must be weakly inhibited and all others strongly. *Settling it:* insect AL connectomics plus paired recordings to sign and weight the graph, then check whether the measured $\rho$ admits heteroclinic cycles at all. `structures/S-03-combinatorial-threshold-linear-networks.md` gives the machinery to answer this from the graph alone, without simulating. A null result here would be the most informative outcome in the field, because it would kill the most beautiful hypothesis.

**Does the covariance spectrum actually flatten?** §2.2 predicts monotone flattening with the top mode decaying first, and a transient duration scaling as $1/\sqrt{\lambda_{\min}}$. *Settling it:* time-resolved eigenspectra on existing locust and zebrafish datasets; then the causal version — rear animals in restricted odour environments and ask whether the transform tracks the new statistics. If it does not adapt, it is not whitening.

**Is decorrelation the objective, or a side effect of sparsening?** C1 showed that expansion-plus-WTA decorrelates all by itself, with no learning and no lateral inhibition. If the KC layer decorrelates anyway, why should the AL bother? Two possibilities: the AL's decorrelation is what makes the KC layer's fixed threshold viable (a C8-style normalisation argument), or the AL's decorrelation is redundant and the normative story is wrong. *Settling it:* quantify how much of the total decorrelation from ORN to KC happens at each stage, in the same preparation, on the same odour panel. This is a straightforward experiment that would discipline a large amount of loose theory, including some of the above.

**What is the fixed point for?** Mazor & Laurent showed the fixed point is *less* informative about identity. It is nonetheless reliably reached and reliably odour-specific. Candidate answers: it is the state from which the off-transient is launched, and the off-transient is itself informative; it is a low-cost maintenance state for prolonged stimuli; it is an artefact of unnaturally long stimulus presentations that the animal never encounters in a turbulent plume. The third is my bet, and it has a clear test: deliver naturalistic intermittent plume statistics and ask whether the AL ever leaves the transient regime at all. If it does not, the entire fixed-point-versus-transient framing is a laboratory artefact — and the correct level-2 statement is about a system that is *always* transient, which sharply favours readings (b) and (c) over any account that needs a settled state.
