---
title: "C7 · Dopamine & TD learning"
parent: Case studies
nav_order: 7
---

# C7 — Dopamine → temporal-difference learning

> **Circuit.** Midbrain dopamine neurons (VTA/SNc) and their projections to striatum and frontal cortex, embedded in the basal-ganglia loop.
> **Primitive extracted.** A scalar *temporal-difference error* $\delta_t = r_t + \gamma V(s_{t+1}) - V(s_t)$, broadcast as a global teaching signal to a value learner (critic) and a policy learner (actor).
> **Status.** The best-validated circuit→algorithm conversion in neuroscience — and now visibly incomplete. The core TD reading survives thirty years of adversarial testing; the claim that dopamine *is* a scalar TD error does not.
> **Structures thread.** Stochastic approximation and fixed points of contraction maps; the projection/least-squares thread (TD as a *semi*-gradient); the three-factor plasticity thread. See `../structures/README.md`.

---

## 1. The phenomenon

Schultz's group recorded from putative dopamine neurons in behaving monkeys through the late 1980s and early 1990s, and found something that looked, at first, like a disappointment. The neurons fired a burst to unpredicted juice. Fine — reward neurons. But once a light reliably predicted the juice, the burst *moved*: it appeared at the light and vanished at the juice. Reward delivery, the thing the animal actually wanted, produced no response at all. And if you then withheld the juice, the neurons fell silent — a pause below baseline, precisely at the time the juice *would* have arrived, with sub-100-ms timing fidelity.

Three signatures, then:

1. **Unpredicted reward → burst.**
2. **Learned prediction → the burst transfers to the earliest reliable predictor.**
3. **Omission → a dip, time-locked to the expected delivery, not to any stimulus.**

Signature 3 is the killer. A "reward detector" cannot dip at a moment when nothing happens. Whatever these cells encode, it is defined relative to an internally maintained expectation with a clock. Ljungberg, Apicella & Schultz (1992) and Schultz, Apicella & Ljungberg (1993) laid out the transfer; Mirenowicz & Schultz (1994) showed the responses are appetitive-specific rather than generic salience (at least for the majority population); the synthesis appeared as Schultz, Dayan & Montague (1997).

Meanwhile Montague, Dayan & Sejnowski (1996) — building on Sutton & Barto's temporal-difference theory and on Barto's earlier basal-ganglia sketch — had the conversion in hand. The three signatures are not three facts. They are one equation.

## 2. The conversion

### 2.1 From the Bellman equation to TD(0)

Fix a Markov policy $\pi$ in a Markov reward process with states $s$, reward $r_t$, and discount $\gamma\in[0,1)$. Define the return and the value function

$$G_t = \sum_{k=0}^{\infty}\gamma^{k} r_{t+k+1}, \qquad V^\pi(s) = \mathbb{E}\!\left[G_t \mid s_t = s\right].$$

Split the first term off the sum:

$$V^\pi(s) = \mathbb{E}\!\left[r_{t+1} + \gamma\sum_{k=0}^\infty \gamma^k r_{t+k+2}\,\middle|\, s_t=s\right] = \mathbb{E}\!\left[r_{t+1} + \gamma V^\pi(s_{t+1}) \mid s_t = s\right],$$

using the tower property and the Markov property. This is the Bellman equation, $V = T^\pi V$, with $T^\pi V = r^\pi + \gamma P^\pi V$. Since $\|T^\pi V - T^\pi U\|_\infty \le \gamma\|V-U\|_\infty$, $T^\pi$ is a $\gamma$-contraction and $V^\pi$ is its unique fixed point.

Now define the **temporal-difference error**

$$\boxed{\;\delta_t = r_{t+1} + \gamma V(s_{t+1}) - V(s_t)\;}$$

Its conditional expectation is $\mathbb{E}[\delta_t \mid s_t=s] = (T^\pi V)(s) - V(s)$, the Bellman residual. So $\mathbb{E}[\delta_t\mid s]=0$ for all $s$ **iff** $V=V^\pi$. The TD error is a zero-mean noise process exactly when your predictions are right, and its sign tells you which way you are wrong. Signatures 1–3 fall out immediately: an unpredicted reward is a positive residual; a fully predicted one is zero; an omission at a time when $V$ was high and $V(s_{t+1})$ drops to zero is a residual of $-V(s_t)$.

Solving $\mathbb{E}[\delta]=0$ by Robbins–Monro stochastic approximation gives tabular **TD(0)**:

$$V(s_t) \leftarrow V(s_t) + \alpha\,\delta_t .$$

With linear function approximation $V_w(s) = w^\top x(s)$ the update is $\Delta w = \alpha\,\delta_t\, x(s_t)$. Note carefully what this is *not*: it is not gradient descent on $\tfrac12\mathbb{E}[\delta^2]$, because the true gradient would include a $-\gamma x(s_{t+1})$ term from differentiating the bootstrap target. TD is a **semi-gradient**: you differentiate the prediction and treat the target as a constant. This is not a technicality — it is why TD converges to a different (and usually better) solution than Bellman-residual minimisation, and why off-policy TD can diverge. On-policy linear TD does converge (Tsitsiklis & Van Roy, 1997) to the fixed point of $\Pi T^\pi$, where $\Pi$ projects onto the span of the features under the stationary-distribution norm $\|\cdot\|_D$, with the guarantee

$$\|V_{\mathrm{TD}} - V^\pi\|_D \;\le\; \frac{1}{1-\gamma}\,\|\Pi V^\pi - V^\pi\|_D .$$

Exercise 6 asks you to prove this in four lines. It matters biologically: it says the quality of the animal's value estimate is bounded by the expressiveness of its *state representation*, amplified by $1/(1-\gamma)$. Bad state representations are catastrophic at long horizons. Hold that thought for §2.5.

### 2.2 Eligibility traces

The one-step target $r_{t+1}+\gamma V(s_{t+1})$ is low-variance and high-bias; the Monte-Carlo target $G_t$ is the reverse. The $\lambda$-return interpolates geometrically,

$$G^\lambda_t = (1-\lambda)\sum_{n=1}^{\infty}\lambda^{n-1}G^{(n)}_t,\qquad G^{(n)}_t = \sum_{k=1}^{n}\gamma^{k-1}r_{t+k} + \gamma^n V(s_{t+n}),$$

and Sutton's (1988) equivalence theorem shows that, for offline updates, regressing on $G^\lambda_t$ is identical to the causal **backward view**

$$e_t = \gamma\lambda e_{t-1} + x(s_t), \qquad \Delta w = \alpha\,\delta_t\, e_t .$$

The eligibility trace $e$ is a decaying synapse-local memory of recent activation, multiplied by a global scalar that arrives later. That is a *three-factor* learning rule, and it is the strongest structural prediction TD makes about synapses: there must be a silent, decaying, presynaptically-tagged eligibility state with a time constant of order seconds. Yagishita et al. (2014) found it — dendritic spines in D1 striatal neurons enlarge only if dopamine arrives in a narrow window (~0.3–2 s) after glutamatergic input, and not before. The theory demanded a mechanism nobody had asked for, and it was there.

### 2.3 Actor–critic, and why the basal ganglia

The critic learns $V$; the actor learns a policy $\pi_\theta(a\mid s)$. The policy-gradient theorem gives

$$\nabla_\theta J(\theta) = \mathbb{E}_{s\sim d^\pi,\,a\sim\pi_\theta}\!\left[\nabla_\theta\log\pi_\theta(a\mid s)\,Q^\pi(s,a)\right].$$

Any state-dependent baseline $b(s)$ can be subtracted for free, because $\mathbb{E}_a[\nabla_\theta\log\pi_\theta(a|s)\,b(s)] = b(s)\nabla_\theta\!\sum_a \pi_\theta(a|s) = b(s)\nabla_\theta 1 = 0$. Take $b=V^\pi$, giving the advantage $A^\pi(s,a) = Q^\pi(s,a)-V^\pi(s)$. And now the punchline:

$$\mathbb{E}[\delta_t \mid s_t=s,a_t=a] = \mathbb{E}[r_{t+1}+\gamma V^\pi(s_{t+1})\mid s,a] - V^\pi(s) = Q^\pi(s,a)-V^\pi(s) = A^\pi(s,a).$$

**The same scalar $\delta$ that trains the critic is an unbiased sample of the advantage that trains the actor.** One broadcast signal, two learners:

$$\Delta w = \alpha\,\delta_t\, x(s_t) \quad\text{(critic)}, \qquad \Delta\theta = \beta\,\delta_t\,\nabla_\theta\log\pi_\theta(a_t\mid s_t)\quad\text{(actor)}.$$

This is why the anatomy is suggestive rather than merely compatible. Dopamine neurons project diffusely; a scalar broadcast is exactly what the algorithm wants. Houk, Adams & Barto (1995) and Joel, Niv & Ruppin (2002) map ventral striatum → critic, dorsolateral striatum → actor. The direct (D1, "Go") and indirect (D2, "NoGo") pathways give the actor an *opponent* parameterisation: Frank's (2005) reading is that D1 medium spiny neurons learn from positive $\delta$ and D2 from negative $\delta$, which is exactly what Shen, Flajolet, Greengard & Surmeier (2008) found at the synaptic level — D1 spines potentiate under a dopamine burst, D2 spines potentiate under a dopamine *dip*, via opposite-signed G-protein cascades. A rectifying nonlinearity in the messenger (firing rates cannot go far below a ~5 Hz baseline) is repaired by splitting the channel in two. Good engineering; good evidence.

### 2.4 Why this is the gold standard

Because it made *quantitative, falsifiable, non-obvious* predictions that were then tested, and could have failed.

- **Blocking.** TD learns from prediction error, not contiguity. So if cue A already predicts reward, adding cue B in compound should produce no learning about B — and dopamine should not respond to B. Waelti, Dickinson & Schultz (2001) confirmed both the behavioural and the dopaminergic blocking. A contiguity-based account is dead here.
- **Error magnitude is graded and signed.** Bayer & Glimcher (2005) regressed firing rates against a model-derived $\delta$ trial by trial: positive errors are linearly encoded; negative errors are compressed by the firing-rate floor (see §2.3).
- **Probability and risk.** Fiorillo, Tobler & Schultz (2003) showed cue responses scale with reward probability and found a slow ramp scaling with variance.
- **Sufficiency, not just correlation.** Optogenetics turned the correlation into a causal claim: Tsai et al. (2009) showed phasic dopamine bursts support conditioning; Steinberg et al. (2013) showed an artificial burst at the time of a fully predicted reward *unblocks* learning about a blocked cue — the theory's prediction, run in reverse; Chang et al. (2016) showed brief inhibition at reward time acts like a negative prediction error.
- **The backward shift itself.** Long treated as qualitative, it was quantified by Amo et al. (2022), who reported a gradual, trial-by-trial retreat of the response toward the cue in mice, of the kind TD with a suitable state representation produces.

That is what a mature level-2 claim looks like: a compact equation, derived from a normative problem, with parameters that transfer across tasks and species, generating predictions that survive attempts to kill them.

### 2.5 The honest complications

**Distributional coding.** Dabney et al. (2020) noticed that if different dopamine neurons scale positive and negative errors *asymmetrically*, the population does not learn the mean — it learns a set of expectiles of the reward distribution. Take channel $i$ with gains $\alpha_i^+,\alpha_i^-$ and update $\Delta V_i \propto \alpha_i^+[\delta_i]_+ - \alpha_i^-[\delta_i]_-$. At the fixed point,

$$\alpha_i^+\,\mathbb{E}\big[(R-V_i)_+\big] = \alpha_i^-\,\mathbb{E}\big[(V_i-R)_+\big],$$

which is precisely the defining condition for the $\tau_i$-expectile of the reward distribution with $\tau_i = \alpha_i^+/(\alpha_i^++\alpha_i^-)$. Two consequences are directly measurable: (i) each neuron's *reversal point* — the reward magnitude at which its response changes sign — should equal its own $V_i$, and (ii) reversal point should be monotonically related to asymmetry $\tau_i$ measured independently from response slopes. Both held in mouse VTA. So dopamine heterogeneity is not noise around a scalar; it is a code for a distribution. The algorithm got richer, not wrong.

**Non-value signals.** A substantial fraction of dopamine neurons respond to aversive and salient events (Matsumoto & Hikosaka, 2009), and dorsal-striatum-projecting dopamine axons carry movement-onset and acceleration signals (Howe & Dombeck, 2016; da Silva et al., 2018; Coddington & Dudman, 2018). Novelty responses require an add-on (Kakade & Dayan, 2002, argued novelty bonuses and shaping bonuses are the only two ways this can be made consistent). Mesolimbic *release* dynamics dissociate from spiking and track motivation on slow timescales (Mohebi et al., 2019; Berke, 2018).

**State inference.** Much apparent TD failure is actually the $\Pi$ in §2.1 — the animal's state representation. Starkweather et al. (2017) showed dopamine responses in a task with hidden trial structure match a TD agent that performs belief-state inference, not one with a naive clock. Gardner, Schoenbaum & Gershman (2018) push this further: dopamine may report prediction errors over a *generalised* successor-like representation, not scalar value. Kim et al. (2020) argue ramps and phasic responses can be unified once the state representation and discounting are handled properly.

**The alternative-account debate.** Jeong et al. (2022) reported dopamine release patterns they argue TD cannot produce, proposing ANCCR — a retrospective, contingency-based mechanism computing adjusted net contingency for causal relations rather than a prospective value error. The claim has been contested on both experimental and modelling grounds (much of the disputed evidence turns on which state representation and which measurement modality — spiking vs. release, dorsal vs. ventral — you grant the TD account). My read: ANCCR is a genuine alternative that has not yet produced the density of *risky* confirmed predictions TD has, and the debate is mostly a debate about state representations wearing a different hat. But you should hold that opinion loosely, and §6 says what would settle it.

---

## 3. Model to build

Two short programs. The first reproduces all three Schultz signatures.

```python
import numpy as np

T, N_TRIALS = 40, 400          # 40 bins/trial (25 ms bins), 400 trials
CUE, REW    = 5, 25            # cue at 125 ms, reward at 625 ms
GAMMA, ALPHA, LAM = 0.98, 0.05, 0.9

def csc(t_cue, T):             # complete serial compound: one unit per bin since cue
    X = np.zeros((T, T))
    for t in range(t_cue, T):
        X[t, t - t_cue] = 1.0
    return X

X = csc(CUE, T)
w = np.zeros(T)
delta_hist = np.zeros((N_TRIALS, T))

for trial in range(N_TRIALS):
    e = np.zeros(T)
    omit = trial >= 300                       # probe block: reward withheld
    for t in range(T):
        r    = 1.0 if (t == REW and not omit) else 0.0
        v_t  = w @ X[t]
        v_t1 = w @ X[t + 1] if t + 1 < T else 0.0
        delta = r + GAMMA * v_t1 - v_t
        e = GAMMA * LAM * e + X[t]
        if not omit:                          # freeze learning during probes
            w += ALPHA * delta * e
        delta_hist[trial, t] = delta
```

**What to look for.** Plot `delta_hist` as a trials × time image. Early trials: a single positive band at $t=25$. Late trials: the band has migrated to $t=4$ — index $t$ stores the transition $s_t\!\to\!s_{t+1}$, so the cue response appears one bin *before* the cue arrives, at value $\gamma^{21}\approx0.654$ — and reward-time $\delta\to 0$. In the omission block, a sharp negative deflection at exactly $t=25$ of magnitude $\to -1$ (Exercise 1). Then set `LAM = 0` and re-run: you will see the response creep backwards one bin per trial, a *visible* travelling bump. Real dopamine does not show that intermediate bump clearly, which is one of the standing embarrassments of the complete-serial-compound representation and motivated microstimulus representations (Ludvig, Sutton & Kehoe, 2008) and the belief-state accounts of §2.5.

Second, the distributional variant — twelve lines that reproduce the core of Dabney et al. (2020):

```python
rng = np.random.default_rng(0)
K = 21
alpha_p = np.linspace(0.05, 0.95, K); alpha_m = 1.0 - alpha_p   # tau_i = alpha_p[i]
V = np.zeros(K)
R = rng.choice([1.0, 5.0], size=200000, p=[0.5, 0.5])           # bimodal rewards
for r in R:
    d = r - V
    V += 0.005 * np.where(d > 0, alpha_p * d, alpha_m * d)
print(np.round(V, 2))
```

**What to look for.** `V` is monotone in $\tau$ and spans (1, 5) without ever visiting the mean-only solution; a channel with $\tau=0.5$ sits at 3.0, low-$\tau$ channels near 1, high-$\tau$ near 5. Each $V_i$ *is* the reward magnitude at which channel $i$'s $\delta$ reverses sign — the measured quantity in the paper. Now run with a *unimodal* reward distribution of the same mean and variance and confirm the expectile spread differs: that is the decodability claim.

---

## 4. Exercises

**★ 1. The omission dip has a predicted magnitude.**
In the deterministic CSC task above with cue at $t_c$ and reward $r=1$ at $t_R$, compute the converged $V(t)$ for $t_c\le t\le t_R$, and hence the value of $\delta$ at each time step on a probe trial where reward is omitted.

<details markdown="1"><summary>Solution</summary>

At the TD fixed point $\mathbb{E}[\delta_t]=0$ everywhere, so $V(t)=\gamma V(t+1)$ for $t<t_R$ and $V(t_R)=r+\gamma V(t_R+1)=r$ (nothing follows). Recursing, $V(t)=\gamma^{\,t_R-t}r$ for $t_c\le t\le t_R$, and $V(t)=0$ for $t>t_R$ and $t<t_c$.

Adopt the code's indexing convention: $\delta_t$ is the error on the transition $s_t\to s_{t+1}$, i.e. $\delta_t = r_{t+1}+\gamma V(t+1)-V(t)$.

On a probe trial, for $t_c\le t<t_R$: $\delta_t = 0+\gamma V(t+1)-V(t) = \gamma\cdot\gamma^{t_R-t-1}r - \gamma^{t_R-t}r = 0$. Nothing happens until the expected moment. At $t=t_R$: $\delta_{t_R} = 0 + \gamma V(t_R+1) - V(t_R) = 0 - r = -r$. And on the transition *into* the cue state, $t=t_c-1$: $\delta_{t_c-1} = 0 + \gamma V(t_c) - V(t_c-1) = \gamma^{\,t_R-t_c+1}r$ — the transferred cue response, discounted by the delay. (With $\gamma=0.98$, $t_c=5$, $t_R=25$: $0.98^{21}=0.654$, which is what the simulation returns.)

So the theory predicts three numbers, not three shapes: cue response $\gamma^{\Delta+1}r$ with $\Delta=t_R-t_c$, reward response $0$, omission dip $-r$ at exactly $t_R$. The ratio of cue response to unpredicted-reward response measures $\gamma^{\Delta+1}$ and therefore *identifies the discount factor* from the data. That is the kind of prediction that makes a theory testable — and note that it is sensitive to the assumed bin width, which is one concrete way the state representation (§2.1) leaks into every quantitative claim.
</details>

**★★ 2. $\delta$ is the advantage.**
Prove $\mathbb{E}[\delta_t\mid s_t=s,a_t=a]=A^\pi(s,a)$ when the critic is exact, and state precisely what bias arises when $V\ne V^\pi$.

<details markdown="1"><summary>Solution</summary>

$\mathbb{E}[\delta_t\mid s,a] = \mathbb{E}[r_{t+1}\mid s,a] + \gamma\mathbb{E}[V^\pi(s_{t+1})\mid s,a] - V^\pi(s)$. The first two terms are exactly the one-step Bellman expansion of $Q^\pi(s,a) = \mathbb{E}[r_{t+1}+\gamma V^\pi(s_{t+1})\mid s,a]$. Hence the expression equals $Q^\pi(s,a)-V^\pi(s)=A^\pi(s,a)$.

With an approximate critic $V = V^\pi + \varepsilon$, the expectation becomes $A^\pi(s,a) + \gamma\mathbb{E}[\varepsilon(s_{t+1})\mid s,a] - \varepsilon(s)$, i.e. a bias equal to the *negative Bellman residual of the error function*. Crucially this bias is not purely state-dependent — it depends on $a$ through the transition distribution — so it does **not** vanish as a baseline. Critic error therefore corrupts the policy gradient direction, not just its scale. Biologically: a ventral-striatal lesion should not merely slow actor learning but bias which actions are learned. That is testable and, broadly, observed.
</details>

**★★ 3. Forward = backward.**
Show, for offline (end-of-episode) updates in the tabular case, that the sum of backward-view TD($\lambda$) updates equals the sum of forward-view $\lambda$-return updates.

<details markdown="1"><summary>Solution</summary>

Backward view accumulates, for state $s$ visited at time $k$, the total $\alpha\sum_{t\ge k}(\gamma\lambda)^{t-k}\delta_t$. So it suffices to show $G^\lambda_k - V(s_k) = \sum_{t\ge k}(\gamma\lambda)^{t-k}\delta_t$.

Write $G^{(n)}_k - V(s_k) = \sum_{j=0}^{n-1}\gamma^{j}\delta_{k+j}$ (telescoping: each $\delta_{k+j}=r_{k+j+1}+\gamma V(s_{k+j+1})-V(s_{k+j})$, and consecutive $\gamma^{j}V$ terms cancel). Then

$$G^\lambda_k - V(s_k) = (1-\lambda)\sum_{n\ge1}\lambda^{n-1}\sum_{j=0}^{n-1}\gamma^j\delta_{k+j} = \sum_{j\ge0}\gamma^j\delta_{k+j}\,(1-\lambda)\!\!\sum_{n\ge j+1}\!\!\lambda^{n-1} = \sum_{j\ge0}(\gamma\lambda)^j\delta_{k+j},$$

since $(1-\lambda)\sum_{n\ge j+1}\lambda^{n-1} = \lambda^{j}$. Swapping the order of summation is justified by absolute convergence for $\gamma\lambda<1$ and bounded rewards. $\square$
</details>

**★★ 4. Expectiles from asymmetric learning rates.**
Show that the asymmetric update $\Delta V \propto \alpha^+[R-V]_+ - \alpha^-[V-R]_+$ has as its unique stable fixed point the $\tau$-expectile with $\tau=\alpha^+/(\alpha^++\alpha^-)$, and compute $V$ for the two-point distribution $R\in\{1,5\}$ with equal probability.

<details markdown="1"><summary>Solution</summary>

Setting the expected update to zero: $\alpha^+\mathbb{E}[(R-V)_+] = \alpha^-\mathbb{E}[(V-R)_+]$, which is the first-order condition for minimising the asymmetric squared loss $\mathbb{E}[\,|\tau - \mathbb{1}\{R<V\}|\,(R-V)^2\,]$ — the definition of the $\tau$-expectile. Uniqueness/stability: the derivative of the expected update w.r.t. $V$ is $-(\alpha^+ P(R>V) + \alpha^- P(R<V)) < 0$, so the map is strictly decreasing in $V$, giving one crossing and local stability.

For $R\in\{1,5\}$ with $p=1/2$ and $1<V<5$: $\alpha^+\cdot\tfrac12(5-V) = \alpha^-\cdot\tfrac12(V-1)$, so $V = \dfrac{5\alpha^+ + \alpha^-}{\alpha^++\alpha^-} = 1 + 4\tau$. Thus $\tau=0.05\Rightarrow V=1.2$; $\tau=0.5\Rightarrow V=3$; $\tau=0.95\Rightarrow V=4.8$. Note the code's printed values should match $1+4\tau$ to within the residual learning-rate noise — a two-line analytic check on a 200 000-step simulation.
</details>

**★★★ 5. Simulate blocking.**
Extend the model of §3 to two cues. Phase 1: cue A alone predicts reward, 200 trials. Phase 2: compound AB predicts the same reward, 200 trials. Phase 3: probe B alone. Predict and then verify the $\delta$ response to B. Then show what a contiguity ("Hebbian cue–reward pairing") rule predicts instead.

<details markdown="1"><summary>Solution</summary>

Represent each cue with its own CSC block and let $V(t) = w_A^\top x_A(t) + w_B^\top x_B(t)$ (linear, shared value). After phase 1, $w_A$ solves the task, so on every compound trial $\delta_t\approx0$ at all $t$: with $\Delta w_B \propto \delta\, e_B$ and $\delta\approx0$, $w_B$ never grows. Probe B alone yields $\delta_{t_B}\approx0$ — no cue response, no learning. This is Kamin blocking, and it is *forced* by the algorithm: the error term is the only teacher.

A contiguity rule $\Delta w_B \propto (\text{B active})\times(\text{reward})$ has no error term, so $w_B$ grows in phase 2 exactly as fast as $w_A$ did in phase 1, and B alone produces a full cue response. The two accounts differ by everything on this experiment. Waelti, Dickinson & Schultz (2001) ran it: dopamine blocked.

Numerically, expect residual non-zero $\delta$ for B of order the phase-1 asymptotic error; with $\alpha=0.05$ and 200 trials this is a few percent of the unblocked response. Report it — the interesting quantitative question is not "zero vs. non-zero" but whether the residual matches the residual TD error, which is a much sharper test.
</details>

**★★★ 6. The TD fixed point is not the least-squares projection.**
With on-policy linear function approximation, prove $\|V_{\mathrm{TD}}-V^\pi\|_D \le \frac{1}{1-\gamma}\|\Pi V^\pi - V^\pi\|_D$, and explain in one sentence why this bound is the mathematical statement of "state representation is destiny".

<details markdown="1"><summary>Solution</summary>

Two facts. (i) $T^\pi$ is a $\gamma$-contraction in $\|\cdot\|_D$ when $D$ is the stationary distribution of $P^\pi$ (Tsitsiklis & Van Roy, 1997), because $\|P^\pi u\|_D \le \|u\|_D$ by Jensen plus stationarity. (ii) $\Pi$ is an orthogonal projection in that norm, hence non-expansive. So $\Pi T^\pi$ is a $\gamma$-contraction with unique fixed point $V_{\mathrm{TD}}$.

Then
$$\|V_{\mathrm{TD}}-V^\pi\|_D \le \|V_{\mathrm{TD}}-\Pi V^\pi\|_D + \|\Pi V^\pi - V^\pi\|_D = \|\Pi T^\pi V_{\mathrm{TD}} - \Pi T^\pi V^\pi\|_D + \|\Pi V^\pi-V^\pi\|_D$$
$$\le \gamma\|V_{\mathrm{TD}}-V^\pi\|_D + \|\Pi V^\pi - V^\pi\|_D,$$
using $V^\pi = T^\pi V^\pi$ and $\Pi V_{\mathrm{TD}} = V_{\mathrm{TD}}$. Rearranging gives the bound. (A tighter constant $1/\sqrt{1-\gamma^2}$ is available.)

Why it matters: the approximation error of the *representation*, $\|\Pi V^\pi - V^\pi\|_D$, is amplified by $1/(1-\gamma)$, which at $\gamma=0.99$ is a factor of 100. Almost every apparent failure of TD as a dopamine theory (ramps, ANCCR-style anomalies, timing) is a candidate failure of $\Pi$, not of $\delta$. Any experiment claiming to refute TD must first pin down the state representation, or it refutes nothing.
</details>

---

## 5. Reading path

- **Sutton (1988)**, *Learning to Predict by the Methods of Temporal Differences* — read it for: the $\lambda$-return, the forward/backward equivalence, and the original convergence argument.
- **Montague, Dayan & Sejnowski (1996)**, *A framework for mesencephalic dopamine systems based on predictive Hebbian learning* — read it for: the conversion itself, stated before most of the confirming data existed.
- **Schultz, Dayan & Montague (1997)**, *A neural substrate of prediction and reward* — read it for: the canonical figure and the three signatures in one place.
- **Waelti, Dickinson & Schultz (2001)**, *Dopamine responses comply with basic assumptions of formal learning theory* — read it for: the blocking experiment, i.e. the moment the theory could have died.
- **Bayer & Glimcher (2005)**, *Midbrain dopamine neurons encode a quantitative reward prediction error signal* — read it for: the regression against a model-derived $\delta$ and the rectification of negative errors.
- **Joel, Niv & Ruppin (2002)**, *Actor–critic models of the basal ganglia: new anatomical and computational perspectives* — read it for: a hard-nosed audit of which actor–critic mappings the anatomy actually permits.
- **Steinberg et al. (2013)**, *A causal link between prediction errors, dopamine neurons and learning* — read it for: unblocking by optogenetic dopamine, the cleanest sufficiency test.
- **Dabney, Kurth-Nelson, Uchida, Starkweather, Hassabis, Munos & Botvinick (2020)**, *A distributional code for value in dopamine-based reinforcement learning* — read it for: how to upgrade a theory using its own heterogeneity as signal.
- **Starkweather, Babayan, Uchida & Gershman (2017)**, *Dopamine reward prediction errors reflect hidden-state inference across time* — read it for: the state-representation move, which you will need in S3.
- **Berke (2018)**, *What does dopamine mean?* — read it for: an honest inventory of everything the scalar-TD reading does not cover.
- **Jeong et al. (2022)**, *Mesolimbic dopamine release conveys causal associations* — read it for: the strongest current alternative account, and practice at evaluating a challenge to a paradigm.

---

## 6. Open problems and what would settle them

**Is the dopamine signal a value error or a general-purpose vector error?** Successor-representation and generalised-prediction-error accounts (Gardner, Schoenbaum & Gershman, 2018) predict dopamine errors about *sensory features*, not just reward. Settling move: measure dopamine responses to surprising, reward-neutral sensory transitions in a task where reward prediction is exactly satisfied, and where the successor representation is independently pinned down by behaviour. A null is informative; a graded, feature-specific response kills the scalar-value reading outright.

**TD vs. ANCCR.** These theories mostly disagree about what the animal's *state* is. Settling move: an experiment in which the two accounts require *opposite-signed* responses under a state representation that is independently constrained (e.g. by having the animal report its belief). Until someone constructs that, the debate will be relitigated on each new dataset.

**What is the eligibility trace, mechanistically, and does its time constant predict behaviour?** Yagishita et al. give a window; nobody has shown that manipulating the window's width shifts the credit-assignment horizon in behaviour as $\lambda$ predicts. Settling move: pharmacological or optogenetic control of the trace's decay, plus a task whose optimal $\lambda$ is known.

**Do the anatomically-distinct dopamine channels implement different algorithms, or the same one on different state spaces?** Settling move: simultaneous recording of ventral and dorsal projection populations with a task that dissociates value error from movement vigour and from novelty, analysed with the expectile decoding of §2.5 applied per channel.

**Is the distributional code read out?** The distributional finding is a claim about the *encoder*. Nothing yet shows a downstream population using the higher moments to make risk-sensitive choices. Settling move: decode the implied distribution from dopamine on single trials, then predict risk attitude trial-by-trial. If risk preference is unrelated to the decoded spread, distributional dopamine is an epiphenomenon of asymmetric plasticity, not an algorithm.

*Next:* `./C9-adaptive-filters.md` takes the opposite tack — a circuit where the *learning rule itself* is the extracted algorithm, and where the teaching signal is not broadcast but locally manufactured.
