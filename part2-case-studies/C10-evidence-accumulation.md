---
title: "C10 · Evidence accumulation"
parent: Case studies
nav_order: 10
---

# C10 — Evidence accumulation → the sequential probability ratio test

> **Circuit.** MT → LIP (and FEF, SC, caudate) in the macaque random-dot motion task; the homologous PPC/FOF circuit in rodents.
> **Primitive extracted.** Accumulate the log-likelihood ratio of the competing hypotheses and stop at a threshold: Wald's sequential probability ratio test, whose continuous limit is the drift-diffusion model. Implemented, approximately, by a two-population mutual-inhibition network operating near a saddle point.
> **Status.** The *algorithm* is on very firm ground — behaviourally, DDM fits are among the best quantitative models in psychology. The *circuit assignment* is much shakier than the textbooks suggest: causal tests of LIP have repeatedly come back negative.
> **Structures thread.** Optimal stopping and dynamic programming; first-passage problems for diffusions; saddle-point normal forms and slow manifolds. The two-population mutual-inhibition motif of §2.4 is the smallest case treated in `../structures/S-03-combinatorial-threshold-linear-networks.md`; index in `../structures/README.md`.

---

## 1. The phenomenon

Show a monkey a patch of dots, a fraction $c$ of which move coherently left or right while the rest are replotted randomly. The animal reports the direction with a saccade. Accuracy is a sigmoid in $c$; reaction time falls with $c$ and is longer for errors than correct choices at the same $c$. MT neurons carry a noisy, roughly *instantaneous* motion signal whose mean scales with $c$ and whose response does not grow over the trial.

LIP does something different. In the reaction-time version (Roitman & Shadlen, 2002), LIP neurons with the target in their response field *ramp*: the slope scales with coherence and with the sign of the eventual choice, and — this is the striking part — the firing rate reaches a common level just before the saccade, essentially independent of coherence and of reaction time. A ramp whose slope encodes stimulus strength and whose terminal value is fixed is the signature of an integrator hitting a bound. Gold & Shadlen (2001, 2007) made the identification explicit: the ramping variable is a running log-likelihood ratio, and the bound is Wald's.

Note what makes this case different from `./C7-dopamine-and-td-learning.md` and `./C9-adaptive-filters.md`. There, the algorithm was extracted from a learning rule or a teaching signal. Here it is fixed by an *optimality theorem* that predates the neuroscience by fifty years. That is a stronger starting point and, as we will see in §6, a more dangerous one: a theorem this elegant makes it very easy to stop asking whether the circuit you are recording from is actually doing the job.

## 2. The conversion

### 2.1 SPRT and its optimality

Observations $x_1,x_2,\dots$ arrive i.i.d. from either $p_1$ ($H_1$: rightward) or $p_2$ ($H_2$: leftward). Define the log-likelihood ratio

$$L_n = \sum_{i=1}^{n}\ell_i,\qquad \ell_i = \log\frac{p_1(x_i)}{p_2(x_i)} .$$

The SPRT continues sampling while $-b < L_n < a$, and stops the moment $L_n \ge a$ (declare $H_1$) or $L_n \le -b$ (declare $H_2$). Let $\alpha = P(\text{choose }H_1\mid H_2)$ and $\beta = P(\text{choose } H_2 \mid H_1)$. Neglecting overshoot,

$$a \approx \log\frac{1-\beta}{\alpha}, \qquad b\approx\log\frac{1-\alpha}{\beta}.$$

**Wald–Wolfowitz (1948).** Among *all* sequential tests with error probabilities no worse than $(\alpha,\beta)$, the SPRT minimises $\mathbb{E}[N]$ under both hypotheses simultaneously. Not "asymptotically", not "on average over priors" — exactly and uniformly. This is a rare thing in optimality theory, and it is why the identification carries weight: if a nervous system evolved to trade speed against accuracy on a two-alternative task with stationary evidence, it should look like an SPRT, and any deviation is a fact requiring explanation rather than a modelling nuisance.

By Wald's identity, $\mathbb{E}[L_N] = \mathbb{E}[N]\,\mathbb{E}[\ell]$, and $\mathbb{E}[\ell\mid H_1] = D_{\mathrm{KL}}(p_1\|p_2)$, so

$$\mathbb{E}[N\mid H_1] \approx \frac{(1-\beta)a-\beta b}{D_{\mathrm{KL}}(p_1\|p_2)} .$$

Decision time is threshold divided by information rate. Everything about the speed–accuracy trade-off is in that one line.

### 2.2 Why accumulating firing rates is accumulating evidence

Suppose the momentary evidence is the difference in pooled activity of two MT populations, and that under $H_{1,2}$ this quantity is Gaussian with means $\pm\mu$ and common variance $\sigma^2$. Then

$$\ell = \log\frac{\mathcal N(x;\mu,\sigma^2)}{\mathcal N(x;-\mu,\sigma^2)} = \frac{-(x-\mu)^2+(x+\mu)^2}{2\sigma^2} = \frac{2\mu}{\sigma^2}\,x .$$

The log-likelihood ratio is *linear in the raw evidence*. So a neuron that simply integrates the MT difference signal, with no knowledge of any likelihood function, computes a scaled LLR — and a fixed threshold on it is a fixed threshold on the LLR. This is the single fact that makes the whole story biologically plausible; it is also exactly where the story is fragile, because the equal-variance Gaussian assumption is doing enormous work (Exercise 1).

### 2.3 The continuous limit: DDM

Take small time steps: $dx = A\,dt + \sigma\,dW$, with absorbing bounds at $\pm z$ and $x(0)=x_0$. Let $p(x)$ be the probability of hitting $+z$ first. The backward Kolmogorov equation is

$$A\,p'(x) + \tfrac{\sigma^2}{2}\,p''(x) = 0,\qquad p(-z)=0,\; p(z)=1 .$$

General solution $p = C_1 + C_2 e^{-kx}$ with $k = 2A/\sigma^2$. Applying the boundary conditions,

$$p(x) = \frac{e^{kz}-e^{-kx}}{e^{kz}-e^{-kz}} \;\;\Longrightarrow\;\; P_{\text{error}} = 1-p(0) = \frac{1-e^{-kz}}{e^{kz}-e^{-kz}} = \frac{1}{1+e^{2Az/\sigma^2}} .$$

$$\boxed{\;P_c = \frac{1}{1+e^{-2Az/\sigma^2}}\;}$$

The psychometric function is a **logistic in $Az/\sigma^2$**, with no free shape parameters. For mean decision time, solve $A\,T'+\tfrac{\sigma^2}{2}T'' = -1$ with $T(\pm z)=0$, giving at $x_0=0$

$$\boxed{\;\bar T = \frac{z}{A}\tanh\!\Big(\frac{Az}{\sigma^2}\Big)\;}$$

Now eliminate the common argument. Since $2P_c-1=\tanh(Az/\sigma^2)$,

$$\bar T = \frac{z}{A}\,(2P_c-1).$$

This is a parameter-linking prediction of the strongest kind: chronometric and psychometric functions are not two curves to be fitted independently: given $z/A$, one determines the other. Fit accuracy across coherences, predict reaction times with a single scale factor. That the DDM survives this test across species and tasks is the real reason to believe it.

Two further exact consequences worth carrying around:

- **Priors and reward biases shift the starting point, not the drift.** Posterior log-odds are $\log(\pi_1/\pi_2) + L_n$, so an unequal prior is implemented by $L_0=\log(\pi_1/\pi_2)$, i.e. $x_0 = \frac{\sigma^2}{2A}\log(\pi_1/\pi_2)$ in evidence units (Exercise 2). This is a discriminating prediction: a starting-point shift produces a bias that *decays* with reaction time, whereas a drift bias produces one that grows. Both have been reported; which you see depends on what you manipulated.
- **Error RTs equal correct RTs** for the pure symmetric DDM. Real data show slow errors (and sometimes fast errors), which is why Ratcliff's (1978) full model adds across-trial variability in drift (→ slow errors) and in starting point (→ fast errors). These are not fudge factors: they are separately identifiable from the RT distribution shapes, and they are the reason the DDM fits rather than merely resembles data.

### 2.4 The circuit: mutual inhibition near a saddle

Wong & Wang (2006) reduced Wang's (2002) spiking attractor network to two variables, the NMDA gating variables $S_1,S_2$ of the two selective populations:

$$\frac{dS_i}{dt} = -\frac{S_i}{\tau_S} + (1-S_i)\,\gamma\,H(x_i), \qquad x_i = J_{11}S_i - J_{12}S_j + I_0 + I_i + I_i^{\text{noise}},$$

with the input–output function $H(x) = (ax-b)/(1-e^{-d(ax-b)})$. The mechanism is: strong recurrent NMDA excitation within each population (slow, $\tau_S\approx100$ ms) plus cross-inhibition between them.

Linearise about the symmetric state $S_1=S_2=\bar S$ at zero coherence and change coordinates to $u=S_1-S_2$, $v = S_1+S_2$. The antisymmetric mode decouples to first order and obeys

$$\dot u = \lambda u + \beta c + \varsigma\,\xi(t), \qquad \lambda = -\frac{1}{\tau_S} - \gamma H(\bar x) + (1-\bar S)\gamma H'(\bar x)\,(J_{11}+J_{12}),$$

with input gain $\beta \propto (1-\bar S)\gamma H'(\bar x)\,J_{A,\text{ext}}\mu_0$. So the network is an **Ornstein–Uhlenbeck process in the choice variable**, and it approximates a perfect integrator only when $\lambda\approx 0$ — i.e. when recurrent excitation plus effective inhibitory gain exactly cancels the leak. Solving,

$$\langle u(t)\rangle = \frac{\beta c}{\lambda}\big(e^{\lambda t}-1\big), \qquad \operatorname{Var}[u(t)] = \frac{\varsigma^2}{2\lambda}\big(e^{2\lambda t}-1\big),$$

which as $\lambda\to0$ reduces to $\beta c t$ and $\varsigma^2 t$: the DDM. Crucially, the *discriminability* saturates for $\lambda \ne 0$ of either sign,

$$\frac{\langle u\rangle}{\sqrt{\operatorname{Var}[u]}} \;\xrightarrow[t\to\infty]{}\; \frac{\beta c}{\varsigma}\sqrt{\frac{2}{|\lambda|}},$$

whereas the perfect integrator's grows as $\sqrt t$ without limit (Exercise 4). Leaky *or* unstable, you stop learning from the stimulus; only the critical case accumulates indefinitely. Fine-tuning is not an aesthetic complaint about the model — it is the model's central quantitative claim, and it is why measurements of accumulator leak in behaving animals matter. (Brunton, Botvinick & Brody, 2013, fitted per-trial accumulation in rats and found leak indistinguishable from zero, which is a genuinely surprising confirmation.)

The attractor implementation buys three things the bare DDM does not have: (i) the bound is not imposed by hand but emerges as the escape from the saddle into one of two stable attractors, so "commitment" is a dynamical state, not a rule; (ii) the instability provides an intrinsic, growing urgency; (iii) $S_i\in[0,1]$ gives saturation, which produces the observed compression of LIP activity at high coherence.

### 2.5 Complications you must not skip

**Collapsing bounds and urgency.** A constant bound is optimal only when the difficulty is known and fixed. If coherence varies unpredictably across trials, optimal stopping under a deadline or an opportunity cost requires a *time-varying* threshold: dynamic programming on the belief state yields bounds that collapse (Drugowitsch, Moreno-Bote, Churchland, Shadlen & Pouget, 2012). Physiologically, Hanks, Kiani & Shadlen (2014) found an evidence-independent, time-growing signal in LIP under speed pressure. Cisek, Puskas & El-Murr (2009) argue for an even stronger version — urgency gating, in which the evidence is *not* integrated at all but low-pass filtered and multiplied by a growing urgency. But Hawkins, Forstmann, Wagenmakers, Ratcliff & Brown (2015) reanalysed a large corpus and concluded that human data mostly favour fixed bounds while monkey data more often favour collapse — a species/training difference, which should make you suspicious that collapse is partly a signature of overtrained animals optimising reward rate.

**Reward rate as the objective.** Bogacz, Brown, Moehlis, Holmes & Cohen (2006) showed that if the animal maximises $\mathrm{RR}(z)=P_c(z)/(\bar T(z)+T_0+(1-P_c)D_p)$, there is a unique optimal $z$, and hence a one-parameter family of optimal (accuracy, RT) pairs — the "optimal performance curve". Testing whether subjects sit on it is a much sharper test than "does the DDM fit".

**More than two alternatives.** For $M$ hypotheses the exactly-optimal test is intractable; MSPRT (Baum & Veeravalli, 1994; Dragalin, Tartakovsky & Veeravalli) thresholds the posterior $\pi_i = e^{L_i}/\sum_j e^{L_j}$ and is only *asymptotically* optimal as error rates $\to0$. Note the algorithmic content: thresholding the posterior is a soft max-versus-*all* comparison, not the max-versus-next race that a naive mutual-inhibition network implements (Exercise 6). Churchland, Kiani & Shadlen (2008) found LIP activity starts lower with four targets than two — consistent with the normalisation term $\log\sum_j e^{L_j}$ having to be overcome, and with Hick's law.

---

## 3. Models to build

**(a) The DDM, analytics verified, plus the reward-rate optimum.**

```python
import numpy as np
rng = np.random.default_rng(0)

def ddm(A, sigma, z, dt=2e-4, n=20000, tmax=8.0):
    x = np.zeros(n); t = np.zeros(n); done = np.zeros(n, bool); ch = np.zeros(n)
    for _ in range(int(tmax/dt)):
        live = ~done
        x[live] += A*dt + sigma*np.sqrt(dt)*rng.standard_normal(live.sum())
        t[live] += dt
        hit = live & (np.abs(x) >= z)
        ch[hit] = np.sign(x[hit]); done |= hit
        if done.all(): break
    return ch, t

A, sigma, z = 0.6, 1.0, 0.8
ch, t = ddm(A, sigma, z); k = A*z/sigma**2          # note: kappa = A z / sigma^2
print("Pc   sim %.3f  theory %.3f" % ((ch > 0).mean(), 1/(1+np.exp(-2*k))))
print("E[T] sim %.3f  theory %.3f" % (t.mean(), (z/A)*np.tanh(k)))
print("RT(correct) %.3f  RT(error) %.3f" % (t[ch > 0].mean(), t[ch < 0].mean()))

T0, Dp = 0.3, 1.0                              # non-decision time, error penalty
zs = np.linspace(0.05, 4.0, 400); kk = A*zs/sigma**2
Pc = 1/(1+np.exp(-2*kk)); ET = (zs/A)*np.tanh(kk)
RR = Pc/(ET + T0 + (1-Pc)*Dp)
print("reward-rate-optimal z = %.3f, accuracy there = %.3f"
      % (zs[np.argmax(RR)], Pc[np.argmax(RR)]))
```

**What to look for.** Simulated $P_c$ and $\bar T$ match the closed forms to within Monte-Carlo error (they should — this is a check on your integrator, not on the theory). Correct and error RTs are equal: this is the *failure* of the pure DDM, and seeing it makes clear why across-trial drift variability is not optional. Finally, note how flat $\mathrm{RR}$ is near its maximum: every bound in $z\in[0.26,0.44]$ — accuracies from 57.7% to 62.9% — sits within 1% of the optimal reward rate. That is a standing warning about how little power "subjects behave optimally" tests actually have, and it is why the *shape* of the optimal performance curve (Exercise 5) is a better test than the location of any one point on it.

**(b) The Wong–Wang two-variable network.**

```python
a, b, d = 270.0, 108.0, 0.154                  # parameters from Wong & Wang (2006)
J11, J12, I0 = 0.2609, 0.0497, 0.3255
gam, tauS, tauA, sig = 0.641, 0.100, 0.002, 0.02
Jext, mu0 = 5.2e-4, 30.0

def H(x):
    y = a*x - b
    return np.where(np.abs(y) < 1e-9, 1.0/d, y/(1.0 - np.exp(-d*y)))

def trial(coh, dt=5e-4, tmax=2.0, thresh=15.0):
    S = np.array([0.1, 0.1]); In = np.zeros(2); t = 0.0
    Is = Jext*mu0*np.array([1 + coh/100.0, 1 - coh/100.0])
    while t < tmax:
        x = np.array([J11*S[0] - J12*S[1], J11*S[1] - J12*S[0]]) + I0 + Is + In
        r = H(x)
        S += dt*(-S/tauS + (1 - S)*gam*r); S = np.clip(S, 0, 1)
        In += -In*dt/tauA + sig*np.sqrt(dt/tauA)*rng.standard_normal(2)
        t += dt
        if r.max() > thresh: return int(np.argmax(r)), t
    return -1, t

for coh in [0, 3.2, 6.4, 12.8, 25.6, 51.2]:
    out = [trial(coh) for _ in range(300)]
    ok = [(c, tt) for c, tt in out if c >= 0]
    print("coh %5.1f  P(correct) %.3f  RT %.3f" %
          (coh, np.mean([c == 0 for c, _ in ok]), np.mean([tt for _, tt in ok])))
```

**What to look for.** A logistic psychometric function and a falling chronometric function, from a network that was never told about likelihood ratios. Then: (i) plot single-trial $S_1$ vs $S_2$ and see trajectories leave the diagonal along the unstable direction — that is the saddle; (ii) fit a DDM to the network's own behaviour and recover $A(c)$, check it is linear in $c$ at low coherence and saturates at high; (iii) scan $J_{11}$ by ±3% and watch performance collapse as $\lambda$ departs from zero — this is the fine-tuning cost made concrete, and it is the strongest argument that some homeostatic mechanism must set $\lambda$.

---

## 4. Exercises

**★ 1. When is integrating rates the same as integrating evidence?**
Derive $\ell = 2\mu x/\sigma^2$ for the equal-variance Gaussian case. Then redo it for $\mathcal N(\mu,\sigma_1^2)$ vs $\mathcal N(-\mu,\sigma_2^2)$ with $\sigma_1\ne\sigma_2$, and say what the circuit would have to do.

<details markdown="1"><summary>Solution</summary>

Equal variance: $\ell = \frac{-(x-\mu)^2 + (x+\mu)^2}{2\sigma^2} = \frac{4\mu x}{2\sigma^2}=\frac{2\mu x}{\sigma^2}$ — linear, so linear accumulation of $x$ suffices.

Unequal variance: $\ell = \log\frac{\sigma_2}{\sigma_1} - \frac{(x-\mu)^2}{2\sigma_1^2} + \frac{(x+\mu)^2}{2\sigma_2^2}$, which contains a term $\propto x^2\big(\tfrac{1}{2\sigma_2^2}-\tfrac{1}{2\sigma_1^2}\big)$. The optimal statistic is **quadratic** in the momentary evidence, so a linear integrator is no longer optimal and there is no bound on the linear accumulator that reproduces SPRT behaviour.

Why this matters: MT firing is approximately Poisson, so variance grows with mean and the two hypotheses genuinely have unequal variances. The saving grace is the *difference* of two opponent pools, which cancels the leading variance asymmetry — i.e. opponency is not merely a convenient way to get a signed signal, it is what makes linear accumulation near-optimal. That is a real algorithmic reason for a real anatomical feature.
</details>

**★★ 2. Priors, biases, and how to tell them apart.**
Show that an unequal prior is implemented exactly by a starting-point offset $x_0 = \frac{\sigma^2}{2A}\log(\pi_1/\pi_2)$. Derive the resulting choice probability with $x_0\ne0$, and show that the induced bias decays with decision time whereas a drift bias does not.

<details markdown="1"><summary>Solution</summary>

Posterior log-odds after $n$ samples are $\log\frac{\pi_1}{\pi_2}+L_n$. Since $L_n = k\,x_n$ with $k=2A/\sigma^2$, running the same symmetric test on $x$ with $x_0 = \frac{1}{k}\log\frac{\pi_1}{\pi_2}$ reproduces it exactly.

From §2.3, $P(\text{hit}+z) = p(x_0) = \dfrac{e^{kz}-e^{-kx_0}}{e^{kz}-e^{-kz}}$. Expanding for small $x_0$: $p(x_0)\approx p(0) + \dfrac{k e^{0}}{e^{kz}-e^{-kz}}x_0$, so the bias term is $\propto \dfrac{k x_0}{2\sinh(kz)}$ — it shrinks exponentially in the bound height.

Time dependence: condition on decision time $t$. For small $t$, only trajectories starting near a bound can have terminated, so the offset dominates and the bias is large; as $t$ grows, the accumulated evidence $At$ dwarfs $x_0$ and the offset's relative contribution $\to0$. Hence **starting-point bias produces choice bias concentrated in fast responses**. A drift offset $A\to A+\Delta$ instead contributes $\Delta t$, so its influence *grows* with $t$. Plot conditional accuracy (or choice bias) as a function of RT quantile: the two mechanisms have opposite slopes. This is the standard, and correct, way to dissociate prior effects from stimulus-gain effects, and it is free — no extra experiment, just a re-binning of data you already have.
</details>

**★★ 3. Wald's identity and the information rate.**
Derive $\mathbb{E}[N\mid H_1]$ in terms of the thresholds and $D_{\mathrm{KL}}(p_1\|p_2)$, and use it to estimate how many independent MT "samples" underlie a 700 ms decision at 6.4% coherence.

<details markdown="1"><summary>Solution</summary>

$N$ is a stopping time with $\mathbb{E}[N]<\infty$ and the $\ell_i$ are i.i.d. with finite mean, so Wald's identity gives $\mathbb{E}[L_N]=\mathbb{E}[N]\,\mathbb{E}[\ell]$. Under $H_1$, $\mathbb{E}[\ell]=D_{\mathrm{KL}}(p_1\|p_2)$. Ignoring overshoot, $L_N \in \{a,-b\}$ with probabilities $(1-\beta,\beta)$, so $\mathbb{E}[L_N\mid H_1]=(1-\beta)a-\beta b$ and

$$\mathbb{E}[N\mid H_1]=\frac{(1-\beta)a-\beta b}{D_{\mathrm{KL}}(p_1\|p_2)}.$$

Numbers: at 6.4% coherence a monkey is right about 80% of the time, so with symmetric errors $\alpha=\beta=0.2$, $a=b=\log 4\approx1.39$, and $\mathbb{E}[L_N]\approx 0.8(1.39)-0.2(1.39)=0.83$ nats. For the Gaussian case, $D_{\mathrm{KL}} = 2\mu^2/\sigma^2 = \tfrac12 d'^2$ per sample. So $\mathbb{E}[N] \approx 1.7/d'^2$ samples. If the decision takes 700 ms and MT provides roughly independent samples every ~50 ms (a fair guess given MT autocorrelation), $N\approx14$, giving $d'\approx0.35$ per sample — a plausible single-sample discriminability for pooled MT at threshold coherence.

The point of the exercise is not the number but the method: SPRT converts behaviour (accuracy + RT) into a *measurement* of the information rate of the sensory input, which you can then check against recordings. Disagreement by an order of magnitude would be a real problem for the theory. It has broadly checked out.
</details>

**★★ 4. Mutual inhibition, leak, and the saturation of discriminability.**
Linearise a two-unit leaky competing accumulator, show it becomes a DDM when leak equals inhibition, and derive the saturation of $\langle u\rangle/\mathrm{sd}(u)$ for $\lambda\ne0$.

<details markdown="1"><summary>Solution</summary>

Take $\dot y_i = -\kappa y_i - \omega y_j + I_i + \varsigma\xi_i$. In coordinates $u=y_1-y_2$, $v=y_1+y_2$: $\dot u = -(\kappa-\omega)u + (I_1-I_2) + \sqrt2\varsigma\xi$, $\dot v = -(\kappa+\omega)v + (I_1+I_2)+\sqrt2\varsigma\xi'$. The modes decouple exactly in the linear case. If $\kappa=\omega$ the $u$ equation loses its restoring term and becomes $\dot u = \Delta I + \sqrt2\varsigma\xi$ — a pure DDM (Bogacz et al., 2006). The $v$ mode remains stable and is discarded by the decision rule. Nonlinearity (rectification, saturation) is what turns $\kappa=\omega$ from a fine-tuned line into a robust regime, because the saddle-point structure survives parameter perturbations even when $\lambda$ does not vanish exactly.

For $\lambda = \omega-\kappa \ne 0$ write $\dot u = \lambda u + \beta c + \varsigma'\xi$. Then $\langle u\rangle = \frac{\beta c}{\lambda}(e^{\lambda t}-1)$ and $\operatorname{Var}(u)=\frac{\varsigma'^2}{2\lambda}(e^{2\lambda t}-1)$ (integrate the OU covariance; both formulas hold for either sign of $\lambda$, and tend to $\beta ct$, $\varsigma'^2t$ as $\lambda\to0$). Hence

$$\frac{\langle u\rangle}{\mathrm{sd}(u)} = \frac{\beta c}{\varsigma'}\cdot\frac{(e^{\lambda t}-1)/\lambda}{\sqrt{(e^{2\lambda t}-1)/(2\lambda)}} \;\xrightarrow[t\to\infty]{}\;\frac{\beta c}{\varsigma'}\sqrt{\frac{2}{|\lambda|}},$$

for both signs (for $\lambda<0$ both numerator and denominator saturate; for $\lambda>0$ the $e^{\lambda t}$ factors cancel). Meanwhile $\lambda=0$ gives $\frac{\beta c}{\varsigma'}\sqrt{t}$, unbounded.

Interpretation: a leaky accumulator forgets old evidence; an unstable one amplifies its own early noise as much as the signal. Both cap the achievable $d'$ at a value set by $|\lambda|$, so measuring the asymptotic accuracy on long-duration trials *measures the leak*. That is exactly what fixed-duration psychophysics with very long stimuli does, and why those experiments are more informative than they look.
</details>

**★★★ 5. The reward-rate-optimal bound.**
Derive the stationarity condition for $z^*$ maximising $\mathrm{RR}(z)=P_c/(\bar T + T_0 + (1-P_c)D_p)$ and solve it numerically for the parameters in §3. Then state the "optimal performance curve" prediction and explain why it is a stronger test than fitting.

<details markdown="1"><summary>Solution</summary>

Write $\kappa = Az/\sigma^2$, $P_c = \sigma(2\kappa)$ (logistic), $\bar T=(z/A)\tanh\kappa$. Then $\mathrm{RR}=P_c/Q$ with $Q=\bar T+T_0+(1-P_c)D_p$, and $\mathrm{RR}'(z)=0$ gives

$$P_c'\,Q = P_c\,Q' \quad\Longleftrightarrow\quad \frac{P_c'}{P_c} = \frac{\bar T' - P_c' D_p}{\bar T+T_0+(1-P_c)D_p},$$

with $P_c' = \frac{2A}{\sigma^2}P_c(1-P_c)$ and $\bar T' = \frac{1}{A}\tanh\kappa + \frac{z}{\sigma^2}\operatorname{sech}^2\kappa$. This is transcendental; solve it numerically (the scan in §3 does exactly that). With $A=0.6,\sigma=1,T_0=0.3,D_p=1$ the scan returns $z^*\approx0.35$, optimal accuracy $\approx0.60$, and $\bar T\approx0.12$ s; raising the error penalty to $D_p=2$ moves it to $z^*\approx0.57$, accuracy $\approx0.66$. Notice how far these are from near-perfect accuracy: with weak evidence and a modest error cost, the reward-rate-optimal agent is right only three times in five. **Optimal agents make errors on purpose**, and any claim that an animal is "suboptimally impulsive" has to be checked against this computation before it means anything.

The optimal performance curve: eliminating $z$ between $P_c(z^*)$ and $\bar T(z^*)$ leaves a relation between normalised decision time $\bar T/(T_0+D_p)$ and error rate that contains **no free parameters at all** once you know $T_0$ and $D_p$, which are set by the experimenter. Subjects either sit on that curve or they do not. Compare this with "the DDM fits my RT distributions": the DDM has enough parameters to fit almost any unimodal, right-skewed RT distribution, so fit quality is weak evidence. A parameter-free curve that the data can miss is strong evidence. This distinction — fitting versus predicting with parameters spent elsewhere — is the subject of `../part3-synthesis/S2-degeneracy-and-limits.md`.
</details>

**★★★ 6. Multi-alternative: the race is not the test.**
For $M$ alternatives, show that thresholding the posterior is not the same as a max-versus-next race, and quantify the difference by simulation for $M=4$.

<details markdown="1"><summary>Solution</summary>

Posterior: $\pi_i = e^{L_i}/\sum_j e^{L_j}$. Requiring $\pi_i>\theta$ is $L_i - \log\sum_j e^{L_j} > \log\theta$, i.e. $L_i - \operatorname{LSE}_j(L_j) > \log\theta$. Since $\operatorname{LSE}$ is a soft maximum over *all* accumulators, the stopping rule depends on the whole population, including the losers. A max-versus-next rule uses only $L_{(1)}-L_{(2)}$; a race-to-a-fixed-bound uses only $L_{(1)}$. These coincide only in the limit where one $L$ dominates, which is the asymptotic ($\alpha\to0$) regime in which MSPRT is provably optimal (Baum & Veeravalli, 1994). At the error rates animals actually work at (10–30%), they differ.

Simulation: draw $M=4$ accumulators with drifts $(\mu,0,0,0)$ and unit noise; implement (i) $\pi_i>\theta$, (ii) $L_{(1)}-L_{(2)}>\theta'$, (iii) $L_{(1)}>\theta''$; tune each threshold to a matched 15% error rate and compare mean decision times. You will find the posterior rule fastest, max-vs-next within a few percent of it, and the pure race clearly worst — because the race ignores the information that the other alternatives are *doing badly*, which is itself evidence. The practical conclusion for neuroscience: divisive normalisation across the choice-selective population is not decoration, it is the $\operatorname{LSE}$ term, and a circuit lacking it is implementing a measurably suboptimal test. Churchland, Kiani & Shadlen's (2008) observation that LIP baseline scales down with the number of targets is the kind of signature you would look for.
</details>

---

## 5. Reading path

- **Wald & Wolfowitz (1948)**, *Optimum character of the sequential probability ratio test* — read it for: the optimality theorem itself, which is what gives this case study its authority.
- **Gold & Shadlen (2001)**, *Neural computations that underlie decisions about sensory stimuli* — read it for: the explicit "LIP accumulates log-likelihood ratio" proposal, stated compactly.
- **Roitman & Shadlen (2002)**, *Response of neurons in the lateral intraparietal area during a combined visual discrimination reaction time task* — read it for: the ramping data, coherence-dependent slopes and coherence-independent terminal rate.
- **Ratcliff & McKoon (2008)**, *The diffusion decision model: theory and data for two-choice decision tasks* — read it for: what the model needs (variability parameters) to fit real RT distributions, and how those parameters are identified.
- **Bogacz, Brown, Moehlis, Holmes & Cohen (2006)**, *The physics of optimal decision making* — read it for: the reduction of LCA/race models to DDM and the reward-rate-optimal threshold.
- **Wong & Wang (2006)**, *A recurrent network mechanism of time integration in perceptual decisions* — read it for: the two-variable reduction and the phase-plane picture; then read Wang (2002) for the spiking original.
- **Gold & Shadlen (2007)**, *The neural basis of decision making* — read it for: the synthesis, and as the document whose claims the later causal work partly undermines.
- **Churchland, Kiani & Shadlen (2008)**, *Decision-making with multiple alternatives* — read it for: the multi-alternative case and the normalisation signature.
- **Drugowitsch, Moreno-Bote, Churchland, Shadlen & Pouget (2012)**, *The cost of accumulating evidence in perceptual decision making* — read it for: how variable difficulty makes collapsing bounds optimal.
- **Brunton, Botvinick & Brody (2013)**, *Rats and humans can optimally accumulate evidence for decision-making* — read it for: per-trial, per-click model fitting; the best available estimate of accumulator leak.
- **Latimer, Yates, Meister, Huk & Pillow (2015)**, *Single-trial spike trains in parietal cortex reveal discrete steps during decision-making* — read it, with the accompanying technical comment and reply, for: how hard it is to tell ramps from steps in averaged data.
- **Katz, Yates, Pillow & Huk (2016)**, *Dissociated functional significance of decision-related activity in the primate dorsal stream* — read it for: the inactivation result that should have changed the field's confidence more than it did.

---

## 6. Open problems and what would settle them

**Is LIP causally necessary?** Katz et al. (2016) inactivated LIP pharmacologically and found essentially intact motion-discrimination performance, while MT inactivation impaired it as expected — despite LIP's textbook decision-related activity. Erlich et al. (2015) reached a similar conclusion for rat PPC. Against this: microstimulation in LIP does bias choices and RTs in the predicted direction (Hanks, Ditterich & Shadlen, 2006), and Zhou & Freedman (2019) reported causal contributions in categorisation tasks. My reading: LIP activity is a *readout* of an accumulation process distributed across a loop (FEF, SC, caudate, and cortex broadly), sufficient to perturb but individually redundant. Settling move: simultaneous multi-area inactivation with dose-response, plus a within-trial causal manipulation timed to specific epochs; and, critically, the prediction that removing *all* candidate nodes should degrade accumulation in a way that a single node's removal does not.

**Ramp or step?** Latimer et al. (2015) argued single trials are better described by discrete jumps than by continuous ramps; the exchange that followed showed both models are hard to distinguish given spiking noise. Settling move: model comparison with matched generative complexity on simultaneously recorded populations, where a population step is far more identifiable than a single-neuron one. This is really a question about statistical power, and it should be treated as such rather than as a debate.

**What sets and maintains $\lambda\approx0$?** The Wong–Wang mechanism requires near-critical tuning; nothing in the model says how it is achieved or held against drift, development, or neuromodulation. Settling move: identify a homeostatic controller and show that perturbing it moves measured leak in behavioural fits. Until then, "the circuit implements a perfect integrator" is a description, not an explanation.

**Are the bounds fixed, collapsing, or urgency-gated, and is the answer species- or training-dependent?** Settling move: within-subject manipulation of difficulty variability with everything else held constant, since that is the variable the normative theory says should switch the optimal policy. A theory that predicts *when* bounds collapse is worth much more than a model class flexible enough to fit either.

*Next:* `../part3-synthesis/S1-rnn-as-model-organism.md` asks what happens when, instead of guessing the algorithm and testing it, you train a network on the animal's task and reverse-engineer what it invented.
