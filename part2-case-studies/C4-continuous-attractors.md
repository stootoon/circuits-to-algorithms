---
title: C4 · Continuous attractors
nav_order: 44
---

# C4 — Ring and line attractors → integration and analog memory

> **Circuit.** The *Drosophila* central complex (E-PG "compass" neurons in the ellipsoid body, P-EN neurons closing an offset loop through the protocerebral bridge), plus the vertebrate oculomotor velocity-to-position integrator and the mammalian head-direction system.
> **Primitive extracted.** A continuous symmetry in the connectivity produces a *manifold* of marginally stable fixed points; asymmetric drive slides the state along that manifold at a rate proportional to the asymmetry. The algorithm is **integration**, and the memory is the position on the manifold.
> **Status.** The most settled circuit→algorithm conversion in modern systems neuroscience. In the fly it is essentially closed: the predicted connectivity motif was found in the connectome, the predicted null direction was found by perturbation, and the predicted topology was found in the population activity. What is *not* settled is how the required fine-tuning is achieved biologically.
> **Structures thread.** [`../structures/S-05-toroidal-topology-of-grid-cells.md`](../structures/S-05-toroidal-topology-of-grid-cells.md) (the ring here is the $d=1$ case of the same empirical-topology programme); [`../structures/S-08-low-rank-connectivity.md`](../structures/S-08-low-rank-connectivity.md) (a cosine ring kernel is exactly a rank-2 perturbation of the identity); [`../structures/README.md`](../structures/README.md), situation 6.

---

## 1. The phenomenon

A tethered fly walks on a ball in the dark. Somewhere in its ellipsoid body — a torus of neuropil about 30 µm across — a single localized bump of calcium activity sits among ~16 wedges of E-PG neurons. When the fly turns right, the bump moves left; when it turns left, the bump moves right; when it stops, the bump stops and stays. In complete darkness, with no landmark, the bump tracks accumulated heading for tens of seconds before it noticeably slips (Seelig & Jayaraman 2015). Turn on a landmark and the bump snaps to a landmark-locked offset that is arbitrary but stable within a session.

That is an *angular path integrator*, running visibly, in a brain you can wire-trace. Three things make it the cleanest case study in this course.

**One state variable** — not a latent recovered by dimensionality reduction with error bars, but a single number readable off the raw imaging frame by eye.

**The predicted anatomy exists.** The theory (Skaggs et al. 1995; Zhang 1996) demanded, a decade before anyone looked in a fly, a recurrent excitatory kernel that is a function of angular difference, global inhibition, and a pair of *offset* pathways gated by angular velocity. E-PG neurons tile the ellipsoid body and the protocerebral bridge; P-EN neurons run bridge→body with a systematic one-wedge offset and receive left/right asymmetric turning signals; Δ7 neurons supply long-range inhibition. The motif is not analogous to the model. It *is* the model (Green et al. 2017; Turner-Evans et al. 2017, 2020; Hulse et al. 2021).

**A structural prediction was tested and passed.** Kim et al. (2017) optogenetically created two bumps and showed they annihilate down to one, and that a jump-started bump relaxes to the stereotyped profile — a test of attractor *structure*, not of tuning.

There is an older, harder cousin: the vertebrate oculomotor integrator, which converts eye-*velocity* commands into the tonic eye-*position* signal that holds your gaze still. Seung (1996) identified it as a **line** attractor and, in the same breath, identified the problem that has occupied the field since: it requires impossible-looking parameter precision.

---

## 2. The conversion

### 2.1 The continuous ring model

Index neurons by preferred heading $\theta \in [-\pi,\pi)$, with synaptic drive $u(\theta,t)$ and firing rate $\phi(u)$. Take

$$\tau\,\partial_t u(\theta,t) \;=\; -u(\theta,t) \;+\; \int_{-\pi}^{\pi}\frac{d\theta'}{2\pi}\,W(\theta-\theta')\,\phi\big(u(\theta',t)\big) \;+\; I(\theta,t).$$

The single structural assumption is that $W$ depends only on $\theta-\theta'$. Everything below is a consequence of it.

**Linear stability of the uniform state.** Put $I \equiv I_0$, linearize $\phi$ to gain $g$ about a uniform solution, and expand $W(\theta)=\sum_n \hat W_n e^{in\theta}$. Convolution is diagonal in the Fourier basis, so the perturbation $\delta u \propto e^{in\theta}$ evolves with rate $(-1+g\hat W_n)/\tau$. The uniform state loses stability in mode $n$ when $g\hat W_n > 1$. For the canonical Mexican-hat kernel

$$W(\theta) = W_0 + W_1\cos\theta, \qquad W_0<0,\;W_1>0,$$

we have $\hat W_0 = W_0$, $\hat W_{\pm1}=W_1/2$, so with $g=1$ the $n=1$ mode goes unstable at $W_1 = 2$ while $n=0$ stays damped because $W_0<0$. A *single* Fourier mode goes unstable, and it is the one that has a free phase. That phase is the whole story.

**The bump, in closed form.** With $\phi(u)=[u]_+$, seek $u^*(\theta)=u_0+u_1\cos\theta$, $u_1>0$, $u_0<0$, active on $|\theta|<\theta_c$ with $\cos\theta_c = -u_0/u_1$. Writing $\cos(\theta-\theta')=\cos\theta\cos\theta'+\sin\theta\sin\theta'$, the fixed-point equation separates:

$$u_0 = I_0 + W_0 A_0,\qquad u_1 = W_1 A_c,\qquad A_0=\!\!\int\!\frac{d\theta}{2\pi}\phi(u^*),\quad A_c=\!\!\int\!\frac{d\theta}{2\pi}\cos\theta\,\phi(u^*).$$

Doing the two integrals over $(-\theta_c,\theta_c)$ and using $u_0=-u_1\cos\theta_c$,

$$A_0 = \frac{u_1}{\pi}\big(\sin\theta_c-\theta_c\cos\theta_c\big),\qquad A_c=\frac{u_1}{2\pi}\big(\theta_c-\sin\theta_c\cos\theta_c\big).$$

The second self-consistency condition therefore reads

$$\boxed{\;W_1\big(\theta_c - \sin\theta_c\cos\theta_c\big) = 2\pi\;}$$

which determines the **bump width from $W_1$ alone** — independent of $I_0$, independent of amplitude. As $\theta_c\to\pi$ this gives $W_1\to 2$, exactly the linear instability threshold; as $\theta_c\to 0$, $W_1 \approx 3\pi/\theta_c^3 \to\infty$. Narrow bumps are expensive. The amplitude then follows:

$$u_1 = \frac{I_0}{-\cos\theta_c - \dfrac{W_0}{\pi}\big(\sin\theta_c - \theta_c\cos\theta_c\big)}.$$

Width fixed by recurrence, gain fixed by input. That is a real, falsifiable prediction and it is roughly what is seen: the fly's bump narrows and sharpens with development of the recurrent circuit but does not fatten when the visual drive is strengthened.

### 2.2 The Goldstone argument: why there is a manifold, not a point

Write the dynamics as $\tau\partial_t u = F[u]$ with $F[u] = -u + W\!*\!\phi(u) + I_0$. Let $T_\psi$ be rotation, $(T_\psi u)(\theta)=u(\theta-\psi)$. Translation invariance of $W$ and constancy of $I_0$ mean exactly that $F$ **commutes** with $T_\psi$:

$$F[T_\psi u] = T_\psi F[u].$$

So if $u^*$ is a fixed point, so is $u^*_\psi \equiv T_\psi u^*$ for every $\psi$: the fixed-point set contains a circle. Differentiate $F[T_\psi u^*]=0$ with respect to $\psi$ at $\psi=0$, using $\partial_\psi T_\psi u^*|_0=-\partial_\theta u^*$:

$$DF[u^*]\cdot \partial_\theta u^* = 0 .$$

**The Jacobian has $\partial_\theta u^*$ in its kernel.** Zero eigenvalue, one per continuous symmetry — the Goldstone mode, arriving here by exactly the argument that gives you massless bosons and phonons. The transverse directions are strictly stable if the bump is stable at all, so the state space near the manifold is a stable circle bundle: fast relaxation onto the ring, zero restoring force *along* it. That is analog memory. No parameter of the model stores $\psi$; $\psi$ is stored because nothing removes it.

The corresponding **left** null vector is worth computing, because we need it in a moment. With $L=DF[u^*]$, $Lv = -v+W\!*\!(\phi'(u^*)v)$, and $W$ even so $L^\dagger v = -v + \phi'(u^*)\,(W\!*\!v)$. Try $v=\phi'(u^*)\,\partial_\theta u^*$:

$$L^\dagger v = \phi'(u^*)\Big[-\partial_\theta u^* + W\!*\!\big(\phi'(u^*)\partial_\theta u^*\big)\Big] = \phi'(u^*)\cdot 0 = 0 .$$

So $v = \phi'(u^*)\partial_\theta u^*$, which for the threshold-linear bump is $-u_1\sin\theta$ restricted to the active region: an **odd** function centred on the bump.

### 2.3 Drift: perturbation theory on the manifold

Perturb the input, $I_0 \to I_0 + \varepsilon\,\delta I(\theta)$, with $\varepsilon\ll1$. The transverse directions are stiff, so the state stays near the manifold and the slow dynamics is entirely in $\psi$. Put $u = u^*_{\psi(t)} + \varepsilon u_1$ with $\dot\psi = O(\varepsilon)$. At order $\varepsilon$:

$$-\tau\dot\psi\;\partial_\theta u^*_\psi \;=\; L_\psi u_1 + \delta I .$$

Project onto the left null vector $v_\psi$; the $L_\psi u_1$ term dies by construction (Fredholm alternative), leaving the **drift law**

$$\boxed{\;\dot\psi \;=\; -\,\frac{\big\langle v_\psi,\;\delta I\big\rangle}{\tau\,\big\langle v_\psi,\;\partial_\theta u^*_\psi\big\rangle}\;},\qquad \big\langle v_\psi,\partial_\theta u^*_\psi\big\rangle = \int\!\frac{d\theta}{2\pi}\,\phi'(u^*)\,(\partial_\theta u^*)^2 > 0 .$$

Read the numerator. Because $v$ is odd about the bump centre, **only the component of the input that is antisymmetric about the current bump position moves the bump.** The symmetric component changes amplitude and is transverse — it is forgotten. This is the algorithmic content of the circuit in one line, and it is why a landmark input positioned *at* the bump does nothing while the same input offset by a few degrees pulls it.

### 2.4 The P-EN offset loop is a velocity term, exactly

Now build the fly's mechanism. E-PG activity is relayed through two populations that project back with anatomical offsets $\pm\delta$, and whose gains are pushed apart by left/right turning signals: $g_\pm = \tfrac12(1\pm\kappa\omega)$ for angular velocity $\omega$. The effective kernel is

$$W_{\rm eff}(\Delta) = \tfrac12(1+\kappa\omega)W(\Delta-\delta) + \tfrac12(1-\kappa\omega)W(\Delta+\delta) = W(\Delta) + \tfrac{\delta^2}{2}W''(\Delta) - \kappa\omega\delta\,W'(\Delta) + O(\delta^3).$$

The $W''$ term renormalizes the width; absorb it. The new term is $-\kappa\omega\delta\,W'$, and since $\partial_\Delta W * \phi = \partial_\theta(W*\phi)$, at the fixed point where $W*\phi(u^*) = u^*-I_0$ it acts as the effective input

$$\delta I = -\kappa\omega\delta\,\partial_\theta u^* .$$

Substituting into the drift law, the shape-dependent factors cancel *identically*:

$$\boxed{\;\dot\psi = \frac{\kappa\delta}{\tau}\,\omega\;}$$

This is the punchline of the whole unit. The bump slides at a rate strictly proportional to angular velocity, with gain $\kappa\delta/\tau$ set by an *anatomical offset* divided by a *membrane time constant* — and independent of bump amplitude, bump width, input strength, and the details of $\phi$. Integrating,

$$\psi(t) = \psi(0) + \frac{\kappa\delta}{\tau}\int_0^t \omega(t')\,dt' ,$$

which is angular path integration. The circuit does not "compute" an integral; the integral is what the manifold coordinate *is*.

For the cosine kernel this is even prettier. Exactly (no small-$\delta$ expansion),

$$W_{\rm eff}(\Delta) = W_0 + W_1\big[\cos\delta\,\cos\Delta + \kappa\omega\sin\delta\,\sin\Delta\big] = W_0 + \tilde W_1\cos(\Delta-\varphi),\quad \varphi=\arctan(\kappa\omega\tan\delta),$$

so the offset loop *rigidly rotates the kernel*, and a rotated kernel drags the bump at $\dot\psi \simeq \varphi/\tau$. The velocity-to-bump-speed transfer function therefore **saturates** as $\arctan$: a prediction that a fast enough turn under-rotates the bump. This is measurable and, in the fly, gain does fall off at high angular velocity.

### 2.5 Structural fragility: heterogeneity destroys the manifold

Everything above rests on exact translation invariance, which no biological circuit has. Perturb $W(\theta-\theta') \to W(\theta-\theta') + \eta\,V(\theta,\theta')$ with $V$ not translation-invariant and $\eta\ll1$. Then the effective input in the drift law depends on where the bump is:

$$\dot\psi = \frac{\eta}{\tau}\,h(\psi),\qquad h(\psi) = -\frac{1}{\mathcal N}\iint\frac{d\theta\,d\theta'}{(2\pi)^2}\,v(\theta-\psi)\,V(\theta,\theta')\,\phi\big(u^*(\theta'-\psi)\big).$$

$h$ is a smooth $2\pi$-periodic function, so the reduced dynamics is a vector field on $S^1$. Two regimes, and no others:

- **$h$ has zeros.** Generically they are simple and, by the Poincaré–Hopf index theorem on $S^1$, come in an **even number**, alternating stable and unstable. A weak random perturbation is dominated by its lowest Fourier component, so the typical outcome is exactly **two** fixed points: the continuum collapses to a bistable system with one preferred heading and one repelling heading. The bump now drifts systematically toward the preferred heading at rate $\eta|h|/\tau$.
- **$h$ is nowhere zero.** Then there are no fixed points at all and the bump rotates forever at a heading-dependent speed — a systematic, monotone heading error even when the animal is still.

Both are catastrophic for a compass, and both are *diagnosable*: the drift is reproducible within an animal and idiosyncratic across animals, with a drift-versus-heading curve you can measure directly. That is the experiment that distinguishes "true continuous attractor plus noise" (diffusive, unbiased, $\langle\Delta\psi^2\rangle=2Dt$) from "broken attractor" (ballistic, biased, $\Delta\psi \propto t$ with sign depending on $\psi$).

**How much precision is needed?** For i.i.d. relative synaptic errors of size $\eta$ across $N_a$ active neurons, the $\sim N_a$ independent contributions to $h$ self-average, giving $|\dot\psi| \sim c\,\eta/(\tau\sqrt{N_a})$ with $c=O(1)$. For the fly, $N_a\approx 8$ wedges and $\tau\approx 50$ ms, so tolerating a drift below $10°/\text{s} = 0.17$ rad/s requires

$$\eta \lesssim 0.17 \cdot \tau\sqrt{N_a} \approx 0.02,$$

i.e. synaptic weights matched to about **2%**. A mammalian head-direction ring with $N_a\sim10^3$ buys a factor $\sqrt{10^3/8}\approx11$ and tolerates ~20% heterogeneity. This $\sqrt{N}$ scaling is, I think, the single most important quantitative fact in this literature: it explains why small-brained systems must actively calibrate (visual anchoring, plasticity at the E-PG↔P-EN synapses) while large ones can coast, and it predicts that lesioning the ring should degrade heading accuracy as $1/\sqrt{N_a}$ before it degrades anything else.

Noise adds an independent problem: with no restoring force along the manifold, activity noise produces pure diffusion, $\langle\Delta\psi^2\rangle = 2Dt$ with $D \propto 1/N$ (Burak & Fiete 2012). Memory decays as $\sqrt{t}$ at best, and only population size fixes it.

### 2.6 The line attractor and the fine-tuning problem

Flatten the ring. For a linear recurrent network $\tau\dot{\mathbf x} = -\mathbf x + W\mathbf x + \mathbf I$, decompose along eigenvectors of $W$ with eigenvalue $\lambda$: each mode relaxes with

$$\tau_{\rm eff} = \frac{\tau}{1-\lambda}.$$

Persistent memory means $\lambda = 1$ *exactly*: a line attractor along that eigenvector. Seung's (1996) observation about the goldfish/primate oculomotor integrator is that this is not an idealization, it is a demand for absurd precision. The integrator holds eye position for $\tau_{\rm eff} \gtrsim 20$ s while neuronal $\tau\approx 100$ ms, so

$$1-\lambda \;=\; \frac{\tau}{\tau_{\rm eff}} \;\approx\; 0.005,$$

**half a percent** tuning of a feedback gain, in a system with turnover, neuromodulation, and development. Aksay and colleagues then did the decisive experiments: recording persistent firing in the goldfish area I integrator, showing graded persistent activity in single units, and pharmacologically detuning the loop to produce leaky ($\lambda<1$) and unstable ($\lambda>1$) integrators with the predicted exponential signatures (Aksay et al. 2001, 2007; Major et al. 2004). The fine-tuning is real, and it is maintained.

Three families of answers, all live. **Learn it:** corrective feedback from retinal slip continuously recalibrates the gain (Seung 1998; Major et al. 2004). **Make the tuning structural rather than parametric:** homeostatic synaptic scaling restores symmetry automatically because the perturbation itself drives the compensation (Renart, Song & Wang 2003), and negative-derivative feedback via balanced E/I makes $\tau_{\rm eff}$ robust to gain error (Lim & Goldman 2013). **Give up exact marginality:** many discrete bistable units — hysteretic dendrites, a staircase of saddle-node ghosts — coarse-grain into an integrator with no parameter needing to equal 1 (Koulakov et al. 2002; Goldman et al. 2003).

I find the third the most honest reading, and it changes the level-2 statement: the algorithm is not "integrate" but "integrate with quantization error $\Delta$, plus a mechanism keeping $\Delta$ below behavioural threshold." Different algorithm, different failure modes — and §2.5's drift-versus-state measurement is how you tell it from the idealization.

### 2.7 The general lesson

Strip the biology and the statement is:

> **A continuous symmetry group $G$ acting on the state space and commuting with the dynamics produces a $\dim G$-dimensional manifold of fixed points, with one zero eigenvalue per generator. Terms in the dynamics that break the symmetry antisymmetrically slide the state along the manifold at a rate given by a projection onto the left null vectors. The algorithm implemented is integration of those terms; the memory is the group coordinate.**

$G=SO(2)$ gives heading. $G=\mathbb{R}$ gives eye position. $G=\mathbb{R}^2/\Lambda$ gives grid-cell phase — see [C6](C6-grid-cells.md), where the *identical* offset-kernel trick appears in two dimensions. The recipe for spotting one of these in your own data: look for a stimulus dimension along which the population response is a rigid translation, then check whether spontaneous activity in the dark visits the same one-parameter family.

---

## 3. Worked example / model to build

A ring attractor with an offset-loop velocity input, plus a switch for connectivity heterogeneity.

```python
import numpy as np

N, tau, dt, T = 128, 0.05, 0.001, 30.0
theta = 2*np.pi*np.arange(N)/N
W0, W1, I0 = -12.0, 4.0, 1.0          # W1 > 2 => bump exists; W0 < 0 => global inhibition
delta, kappa = 2*np.pi/N*4, 0.8       # anatomical offset (4 wedges) and velocity gain
rng = np.random.default_rng(0)

eta = 0.0                             # <-- set to 0.02, 0.05, 0.10 for the fragility runs
het = eta*rng.standard_normal((N, N))

def kernel(shift):
    D = theta[:, None] - theta[None, :] - shift
    return (W0 + W1*np.cos(D))/N * (1.0 + het)

Kp, Km = kernel(+delta), kernel(-delta)   # W(Delta - delta), W(Delta + delta)
phi = lambda u: np.maximum(u, 0.0)

t     = np.arange(0, T, dt)
omega = 0.6*np.sin(2*np.pi*t/10.0)        # angular velocity, rad/s
u     = np.exp(3.0*np.cos(theta)); u /= u.max()
psi_b, psi_t, ang = np.zeros_like(t), np.zeros_like(t), 0.0

for n in range(len(t)):
    g = kappa*omega[n]
    W = 0.5*(1+g)*Kp + 0.5*(1-g)*Km
    u = u + (dt/tau)*(-u + W @ phi(u) + I0)
    ang += omega[n]*dt
    psi_t[n] = ang
    psi_b[n] = np.angle(phi(u) @ np.exp(1j*theta))   # population vector average

psi_b = np.unwrap(psi_b)
print("predicted gain kappa*delta/tau =", kappa*delta/tau)
print("measured gain                  =", np.polyfit(psi_t, psi_b, 1)[0])
print("bump width theta_c (deg)       =", np.degrees(np.sum(u > 0)*2*np.pi/N/2))
print("residual drift over 30 s (deg) =", np.degrees((psi_b - psi_t)[-1] - (psi_b - psi_t)[0]))
```

**What to look for.**
1. With `eta = 0` the measured gain is $3.18$ against a predicted $\kappa\delta/\tau=3.14$, and residual drift over 30 s is $0.4°$. You have built an integrator whose gain is an anatomical offset divided by a time constant.
2. Bump width is insensitive to `I0`: double it and the bump gets taller, not wider. Change `W1` from 4 to 6 and watch $\theta_c$ shrink as $W_1(\theta_c-\sin\theta_c\cos\theta_c)=2\pi$ predicts (at $W_1=4$ this gives $\theta_c\approx90°$; the code reports $91.4°$).
3. Turn on `eta`, set `omega` to zero, start the bump at 16 evenly spaced phases and run 10 s. At $\eta=0$ all 16 stay put — that is the marginal manifold. At $\eta=0.05$ they clump; at $\eta=0.10$ they collapse onto **two** roughly antipodal clusters, exactly the $n=1$-dominated two-fixed-point prediction of §2.5. Plot initial phase against drift rate: that is $h(\psi)$, measured.
4. Scan $\eta\in\{0.01,\dots,0.2\}$ and confirm drift $\propto\eta$; then double $N$ and confirm the $1/\sqrt{N_a}$ improvement. The fine-tuning problem, quantified in twenty lines.
5. Raise $\kappa\omega$ until $\kappa\omega\tan\delta\sim1$ and watch the gain saturate as $\arctan$.

---

## 4. Exercises

**★ E4.1 — Instability threshold and mode selection.** For $W(\theta)=W_0+W_1\cos\theta$ and gain $g$, find the condition for the uniform state to be unstable *only* to the $n=1$ mode, and explain why $n=2$ instability would be a disaster for a compass.

<details markdown="1"><summary>Solution</summary>

Fourier coefficients: $\hat W_0=W_0$, $\hat W_{\pm1}=W_1/2$, $\hat W_n=0$ for $|n|\ge2$. Growth rate of mode $n$ is $(-1+g\hat W_n)/\tau$. So:
- $n=0$ stable requires $gW_0<1$: automatic for $W_0<0$.
- $n=1$ unstable requires $gW_1/2>1$, i.e. $W_1>2/g$.
- $n\ge2$: rate $=-1/\tau<0$ always, since $\hat W_n=0$.

So the pure cosine kernel selects $n=1$ by construction. If instead $\hat W_2$ were also supercritical, the unstable manifold would contain a two-bump family. Two bumps on a ring is not a compass: the population vector average is degenerate (two antipodal bumps give zero net vector), and the symmetry group of the pattern is $\mathbb{Z}_2\ltimes SO(2)$, so heading is only recoverable modulo $\pi$. Mechanistically this is exactly what Kim et al. (2017) tested by inducing two bumps: a genuine ring attractor must *annihilate* the second one, which requires the $n=1$ mode to be the only unstable one and the inhibition ($W_0$) to be strong enough to prevent coexistence.
</details>

**★★ E4.2 — The zero mode, from scratch.** Without assuming any particular $\phi$ or $W$, prove that the Jacobian at a bump fixed point of a translation-invariant ring network has $\partial_\theta u^*$ in its kernel, and identify the left null vector. Then show that if $W$ is *not* even ($W(\theta)\ne W(-\theta)$) the right null vector is unchanged but the left one is not.

<details markdown="1"><summary>Solution</summary>

**Right null vector.** $F[u]=-u+W*\phi(u)+I_0$ commutes with rotation: $(W*\phi(T_\psi u))(\theta) = \int \frac{d\theta'}{2\pi}W(\theta-\theta')\phi(u(\theta'-\psi))$; substitute $\theta''=\theta'-\psi$ to get $\int\frac{d\theta''}{2\pi}W(\theta-\psi-\theta'')\phi(u(\theta''))=(W*\phi(u))(\theta-\psi)$. Hence $F[T_\psi u]=T_\psi F[u]$. With $F[u^*]=0$ this gives $F[T_\psi u^*]=0$ for all $\psi$. Differentiate at $\psi=0$: $DF[u^*]\,\partial_\psi(T_\psi u^*)|_0=0$, and $\partial_\psi u^*(\theta-\psi)|_0=-\partial_\theta u^*$, so $L\,\partial_\theta u^*=0$ with $L=DF[u^*]$. Note this used *no* property of $\phi$ or of $W$ beyond translation invariance.

**Left null vector.** $Lw = -w + W*(\phi'(u^*)w)$. The adjoint under $\langle a,b\rangle=\int \frac{d\theta}{2\pi}ab$: $\langle v, W*(\phi' w)\rangle = \langle \check W * v, \phi' w\rangle = \langle \phi'\,(\check W*v), w\rangle$ where $\check W(\theta)=W(-\theta)$. So $L^\dagger v = -v + \phi'(u^*)(\check W * v)$.

If $W$ is even, $\check W = W$, and $v=\phi'(u^*)\partial_\theta u^*$ works: $L^\dagger v = \phi'\big[-\partial_\theta u^* + W*(\phi'\partial_\theta u^*)\big]=\phi'\cdot L(\partial_\theta u^*)=0$.

If $W$ is not even, write $W = W_e + W_o$. The right null vector is still $\partial_\theta u^*$ (the argument above never used evenness), but the trial $v$ now gives $L^\dagger v = \phi'\big[-\partial_\theta u^* + \check W*(\phi'\partial_\theta u^*)\big] = \phi'\big[(\check W - W)*(\phi'\partial_\theta u^*)\big] = -2\phi'\big[W_o*(\phi'\partial_\theta u^*)\big] \neq 0$. So $v$ must be found from $v=\phi'(u^*)(\check W*v)$ directly. Consequence: with an odd component in $W$ the fixed points are still a circle, but the *drift law's* numerator projects onto a different function — and indeed an odd $W_o$ is precisely the $\kappa\omega\delta W'$ term of §2.4, which is why it drives motion.
</details>

**★★ E4.3 — Exact rotation for the cosine kernel.** Show that the offset-loop construction with a pure cosine kernel produces an exactly rotated cosine kernel, derive $\varphi=\arctan(\kappa\omega\tan\delta)$, and show that the resulting bump speed is $\varphi/\tau$ to leading order in $\varphi$. What is the maximum achievable bump speed?

<details markdown="1"><summary>Solution</summary>

$W(\Delta)=W_0+W_1\cos\Delta$. Then
$$\tfrac12(1+g)W(\Delta-\delta)+\tfrac12(1-g)W(\Delta+\delta) = W_0 + \tfrac{W_1}{2}\big[(1+g)\cos(\Delta-\delta)+(1-g)\cos(\Delta+\delta)\big].$$
Expand: $\cos(\Delta\mp\delta)=\cos\Delta\cos\delta\pm\sin\Delta\sin\delta$. The $\cos\Delta$ terms give $\tfrac{W_1}{2}\cdot 2\cos\delta\cos\Delta$; the $\sin\Delta$ terms give $\tfrac{W_1}{2}\cdot 2g\sin\delta\sin\Delta$. So
$$W_{\rm eff}(\Delta)=W_0+W_1\big[\cos\delta\cos\Delta + g\sin\delta\sin\Delta\big] = W_0 + \tilde W_1\cos(\Delta-\varphi),$$
$$\tilde W_1 = W_1\sqrt{\cos^2\delta + g^2\sin^2\delta},\qquad \tan\varphi = g\tan\delta,\quad g=\kappa\omega.$$

**Speed.** Look for a travelling solution $u(\theta,t)=U(\theta-ct)$. Substituting, $-\tau c\,U' = -U + (W_\varphi * \phi(U))$ where $W_\varphi(\Delta)=W(\Delta-\varphi)$, and $(W_\varphi*\phi(U))(\xi)=(W*\phi(U))(\xi-\varphi)$. Take $U\approx u^*$, the unrotated fixed point (valid to $O(\varphi)$), and use $W*\phi(u^*)=u^*-I_0$:
$$-\tau c\,\partial_\theta u^* = -u^*(\xi) + u^*(\xi-\varphi) - I_0 + I_0 = -\varphi\,\partial_\theta u^*(\xi) + O(\varphi^2).$$
Hence $c=\varphi/\tau + O(\varphi^2)$.

**Maximum speed.** $\varphi=\arctan(\kappa\omega\tan\delta) \to \pi/2$ as $\omega\to\infty$, so $c_{\max}\to \pi/(2\tau)$: with $\tau=50$ ms that is $\approx 31$ rad/s. But note also that $\tilde W_1$ *grows* with $g$ when $\delta<\pi/4$... no: $\tilde W_1 = W_1\sqrt{\cos^2\delta+g^2\sin^2\delta}$ increases with $g$, so the bump also narrows (larger $W_1$, smaller $\theta_c$) at high turning speed. Both are testable: gain saturation and speed-dependent bump sharpening.
</details>

**★★ E4.4 — Tuning precision for the oculomotor integrator.** A leaky integrator with membrane $\tau=100$ ms must hold eye position with a time constant of at least 20 s. (a) What relative precision on the loop gain is required? (b) If the gain drifts by 1% per day, how often must recalibration occur? (c) Repeat for a network of 100 units where the required eigenvalue is the top eigenvalue of $W$; if each synapse has independent 10% error, does averaging save you?

<details markdown="1"><summary>Solution</summary>

(a) $\tau_{\rm eff}=\tau/(1-\lambda) \ge 20$ s $\Rightarrow 1-\lambda \le 0.005$. Precision $0.5\%$.

(b) A 1%/day drift exceeds the 0.5% budget in about half a day, so recalibration must run on a timescale of hours or faster. Behaviourally this is what is seen: goldfish integrator time constants are re-trained by imposed retinal slip within tens of minutes (Major et al. 2004). The integrator is not a fixed circuit; it is a continuously-calibrated one.

(c) Write $W = W_{\rm ideal}+\eta E$ with $E$ i.i.d. zero-mean, $\mathrm{sd}(\eta E_{ij})=\eta \bar w$. If the target mode is $\mathbf u$ (unit norm), first-order perturbation theory gives $\delta\lambda = \mathbf u^\top(\eta E)\mathbf u = \eta\sum_{ij}u_iE_{ij}u_j$, a sum of $N^2$ zero-mean terms with weights $u_iu_j$ satisfying $\sum_{ij}u_i^2u_j^2=1$. So $\mathrm{sd}(\delta\lambda)=\eta\bar w$ — **independent of $N$**. Averaging does *not* save you at first order for a fixed mode: the eigenvalue error is the same as a single-synapse error.

But the relevant quantity for a distributed integrator is whether *some* eigenvalue is near 1 and the readout aligns with it. Second-order effects and the freedom to choose the readout direction do help; more importantly, homeostatic mechanisms that act on total drive rather than individual synapses can null $\delta\lambda$ directly. The lesson: unlike the ring attractor's drift (which self-averages as $1/\sqrt{N_a}$ because the errors enter through a *projection onto an odd function*), the line attractor's leak does not self-average, because it enters through a quadratic form aligned with the mode itself. This is a genuine structural difference between the two cases and is why the ring is more forgiving.
</details>

**★★ E4.5 — Diffusion of the bump.** Add independent white noise $\sigma\,\xi_i(t)$ to each of $N$ neurons in a ring attractor. Show that the bump position undergoes free diffusion and that $D\propto \sigma^2/N$ for a bump whose shape is held fixed.

<details markdown="1"><summary>Solution</summary>

Discretize: $\tau\dot u_i = -u_i + \sum_j W_{ij}\phi(u_j) + I_0 + \sigma\sqrt{\tau}\,\xi_i(t)$, $\langle\xi_i(t)\xi_j(t')\rangle=\delta_{ij}\delta(t-t')$. Apply the drift law with $\delta I_i = \sigma\sqrt{\tau}\xi_i$:
$$\dot\psi = -\frac{1}{\tau\mathcal N}\sum_i v_i(\psi)\,\sigma\sqrt{\tau}\,\xi_i(t),\qquad \mathcal N=\sum_i v_i \partial_\theta u^*_i .$$
This is a Wiener process with $\langle\dot\psi(t)\dot\psi(t')\rangle = \frac{\sigma^2}{\tau\mathcal N^2}\|v\|^2\,\delta(t-t')$, so $\langle\Delta\psi^2\rangle=2Dt$ with
$$D = \frac{\sigma^2\|v\|^2}{2\tau\,\mathcal N^2}.$$
Now scale. Let the bump profile be fixed and let $N$ neurons tile it, so $v_i = \bar v(\theta_i)$ with $\bar v$ an $N$-independent function. Then $\|v\|^2=\sum_i \bar v(\theta_i)^2 \propto N$ and $\mathcal N = \sum_i \bar v(\theta_i)\partial_\theta \bar u(\theta_i) \propto N$. Hence
$$D \propto \frac{\sigma^2 N}{N^2} = \frac{\sigma^2}{N}.$$
Two consequences. (i) Angular memory in the dark degrades as $\sqrt{2Dt}$, so halving the error requires quadrupling the population — a brutal scaling that Burak & Fiete (2012) turn into a fundamental bound. (ii) Diffusion is *unbiased*: $\langle\Delta\psi\rangle=0$. Heterogeneity, by contrast, gives $\langle\Delta\psi\rangle\propto t$ with a heading-dependent sign. Measuring the first two moments of drift versus heading therefore separates the two mechanisms cleanly. Do this experiment.
</details>

**★★★ E4.6 — How many fixed points does heterogeneity leave?** Show that a generic small perturbation leaves an even number of fixed points on the ring, that the expected number for a perturbation with Fourier content up to order $n$ is bounded by $2n$, and estimate the probability that *no* fixed point survives.

<details markdown="1"><summary>Solution</summary>

The reduced dynamics is $\dot\psi=\frac{\eta}{\tau}h(\psi)$ with $h:S^1\to\mathbb R$ smooth and $2\pi$-periodic. Fixed points are zeros of $h$.

**Evenness.** $h$ is continuous and periodic, so between consecutive zeros it has constant sign, and going once around the circle the sign must return to itself. Hence sign changes come in pairs: the number of simple zeros is even. Equivalently, by Poincaré–Hopf on $S^1$ (Euler characteristic $0$), $\sum_{\text{zeros}}\mathrm{index}=0$, and simple zeros have index $\pm1$, so there are equally many stable ($h'<0$) and unstable ($h'>0$) ones.

**Bound.** If $h(\psi)=\sum_{|n'|\le n}h_{n'}e^{in'\psi}$ is a real trigonometric polynomial of degree $n$, then $h(\psi)e^{in\psi}$ is a polynomial of degree $2n$ in $z=e^{i\psi}$, which has at most $2n$ roots. So at most $2n$ zeros on the circle. A perturbation whose kernel-error has smooth spatial statistics has rapidly decaying $h_{n'}$, so in practice $n=1$ dominates and you get **two** fixed points — one preferred heading, one repeller.

**Probability of no fixed point.** With $n=1$, $h(\psi)=h_0 + A\cos(\psi-\psi_0)$. Zeros exist iff $|h_0|\le A$. If $h_0$ and $A$ come from independent Gaussian contributions with variances $s_0^2$ and $s_1^2$ (so $A$ is Rayleigh with parameter $s_1$),
$$P(\text{no fixed point}) = P(|h_0|>A) = \int_0^\infty \frac{a}{s_1^2}e^{-a^2/2s_1^2}\,\mathrm{erfc}\!\Big(\frac{a}{s_0\sqrt2}\Big)\,da = 1-\frac{s_1}{\sqrt{s_0^2+s_1^2}} .$$
(The integral is standard: $P(|h_0|>A)=\mathbb E_A[\mathrm{erfc}(A/s_0\sqrt2)]$, and evaluating gives the stated form.) So if the rotationally-symmetric part of the error ($h_0$, a net bias) is comparable to the modulated part, persistent rotation is common. Concretely: an animal with a systematic left-turning bias in the dark is not necessarily behaving oddly — its attractor may simply have $|h_0|>A$. That is a testable and, as far as I know, untested prediction.
</details>

**★★ E4.7 — Computational: measure the drift law.** Using the Section 3 code, verify $\dot\psi = \kappa\delta\omega/\tau$ over a range of $\omega$, then verify that a *symmetric* input perturbation centred on the bump produces zero drift while an antisymmetric one produces drift proportional to its amplitude.

<details markdown="1"><summary>Solution</summary>

Set `eta=0`. Replace the velocity drive with a static extra input `dI = eps*f(theta - psi_now)` applied each step, and measure $d\psi/dt$ over 2 s.

- **Symmetric probe**, $f(\chi)=\cos\chi$ (peaked at the bump centre): measured drift $\approx 0$ to numerical precision, while the bump amplitude changes by $O(\varepsilon)$. Reason: $\langle v,\delta I\rangle = \int \phi'(u^*)(-u_1\sin\chi)\cos\chi\,d\chi/2\pi = 0$ by parity on the symmetric active interval $(-\theta_c,\theta_c)$.
- **Antisymmetric probe**, $f(\chi)=\sin\chi$: drift $\dot\psi = -\varepsilon\langle v,\sin\rangle/(\tau\mathcal N)$. With $v=-u_1\sin\chi\,\mathbb 1_{|\chi|<\theta_c}$ and $\partial_\theta u^*=-u_1\sin\chi$:
$$\langle v,\sin\rangle = -\frac{u_1}{2\pi}\!\!\int_{-\theta_c}^{\theta_c}\!\!\sin^2\!\chi\,d\chi = -\frac{u_1}{4\pi}(2\theta_c-\sin2\theta_c),\quad \mathcal N=\frac{u_1^2}{4\pi}(2\theta_c-\sin2\theta_c),$$
so $\dot\psi = +\varepsilon/(\tau u_1)$: linear in $\varepsilon$, and inversely proportional to bump amplitude. That last dependence is a sharp, non-obvious prediction — a *stronger* bump is harder to push — and it is why the velocity mechanism of §2.4, whose $u_1$-dependence cancels, is the right way to build an integrator. An input-based velocity signal would have a gain that drifts with arousal; a connectivity-based one does not.
- **Velocity sweep**: plot measured $\dot\psi$ against $\omega$ for $\omega\in[0,5]$ rad/s. Linear at low $\omega$ with slope $\kappa\delta/\tau$, bending over as $\arctan(\kappa\omega\tan\delta)/(\tau\omega)$ predicts.
</details>

---

## 5. Reading path

- **Amari (1977)**, *Dynamics of pattern formation in lateral-inhibition type neural fields* — read it for: the original existence-and-stability analysis of localized bumps in a continuous field; everything in §2.1 is in here.
- **Ben-Yishai, Bar-Or & Sompolinsky (1995)**, *Theory of orientation tuning in visual cortex* (PNAS) — read it for: the cosine-kernel ring model solved in closed form, and the marginal phase as a computational resource rather than a bug.
- **Skaggs, Knierim, Kudrimoti & McNaughton (1995)**, *A model of the neural basis of the rat's sense of direction* (NIPS) — read it for: the first explicit offset-loop proposal for angular velocity integration. The fly circuit is this diagram.
- **Zhang (1996)**, *Representation of spatial orientation by the intrinsic dynamics of the head-direction cell ensemble: a theory* (J. Neurosci.) — read it for: the clean derivation of the drift law and of what asymmetric input must look like; the mathematical core of this unit.
- **Seung (1996)**, *How the brain keeps the eyes still* (PNAS) — read it for: the line-attractor hypothesis and, more importantly, the fine-tuning problem stated with no hedging.
- **Aksay, Gamkrelidze, Seung, Baker & Tank (2001)**, on persistent activity in the goldfish oculomotor integrator (Nat. Neurosci.) — read it for: what graded persistent firing actually looks like in single neurons.
- **Goldman, Levine, Major, Tank & Seung (2003)**, on robust persistent activity from hysteretic dendrites (Cereb. Cortex) — read it for: the "ladder instead of a line" alternative. Take it seriously.
- **Renart, Song & Wang (2003)**, *Robust spatial working memory through homeostatic synaptic scaling* (Neuron) — read it for: symmetry restored by a mechanism driven by the very asymmetry it removes. The most elegant answer to §2.5.
- **Seelig & Jayaraman (2015)**, *Neural dynamics for landmark orientation and angular path integration* (Nature) — read it for: the bump, in a fly, in the dark.
- **Kim, Rouault, Druckmann & Jayaraman (2017)**, *Ring attractor dynamics in the Drosophila central brain* (Science) — read it for: perturbation experiments that test attractor structure rather than describing tuning.
- **Green et al. (2017)** (Nature) and **Turner-Evans et al. (2017)** (eLife) — read together for: two independent identifications of P-EN neurons as the velocity-gated offset pathway.
- **Peyrache, Lacroix, Petersen & Buzsáki (2015)** (Nat. Neurosci.) and **Chaudhuri, Gerçek, Pandey, Peyrache & Fiete (2019)** (Nat. Neurosci.) — read for: the ring persists in sleep, with no sensory input. Attractor, not stimulus tuning.
- **Burak & Fiete (2012)**, on fundamental limits on persistent activity in noisy networks (PNAS) — read it for: the $D\propto1/N$ bound.
- **Hulse et al. (2021)**, the central-complex connectome analysis (eLife) — read it for: everything the simple ring model leaves out.

---

## 6. Open problems and what would settle them

**1. How is the fine-tuning maintained?** Section 2.5 says the fly needs ~2% weight matching; nothing in development guarantees that. Candidates: plasticity at E-PG↔P-EN synapses driven by visual heading error; homeostatic scaling; or the possibility that the fly simply *does* have a broken attractor with a few discrete headings and compensates by anchoring to visual landmarks whenever available. **What would settle it:** measure drift-versus-heading, in the dark, in individual flies, across days. A true continuous attractor gives unbiased diffusion with $\langle\Delta\psi\rangle=0$ at all $\psi$. A broken one gives a reproducible $h(\psi)$ curve — and if that curve *changes* after visual experience with a specific heading, plasticity is the answer.

**2. Is the mammalian head-direction system one ring or many?** Ring topology has been recovered in ADn (Chaudhuri et al. 2019), but whether the postsubiculum/ADn/lateral mammillary loop is one attractor or several coupled ones with different gains is open. **What would settle it:** simultaneous recording across the loop, with topological analysis of the *joint* state space — one ring or a product of rings? The [`S-05`](../structures/S-05-toroidal-topology-of-grid-cells.md) methodology, one dimension down.

**3. Where does the velocity gain get calibrated?** $\kappa\delta/\tau$ must match real angular velocity in the animal's own units, under changing body mechanics. **What would settle it:** chronically alter the mapping from motor command to actual rotation (closed-loop gain manipulation in the fly) and ask whether bump velocity gain re-calibrates, on what timescale, and whether it survives darkness.

**4. Is marginality the right idealization?** If the true state space is a fine ladder rather than a continuum, the level-2 description must include the quantization. **What would settle it:** the *stationary distribution* of remembered values after long delays. A continuum gives a spreading Gaussian; a ladder gives multimodality at a characteristic spacing. The fly, with 16 wedges, is the place to try.

**5. What is the actual objective?** This unit assumes the level-1 problem is "estimate heading." But an animal integrating noisy angular velocity should be doing *Bayesian* integration, discounting old evidence and weighting landmarks by reliability — see [M4](../part1-foundations/04-probabilistic-computation.md) and [M6](../part1-foundations/06-control-and-filtering.md). A pure attractor is a flat-prior, no-forgetting integrator, which is the wrong estimator. **What would settle it:** manipulate landmark reliability and ask whether snap-to-landmark strength scales as the Kalman gain predicts. If it does, the algorithm is "run a circular Kalman filter," and the attractor is only its prediction step.
