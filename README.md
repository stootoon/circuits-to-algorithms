---
title: Home
nav_order: 1
permalink: /
---

# Circuits to Algorithms

A self-study course in extracting **algorithms** from **neural circuits**.

The premise: a circuit is an implementation, an algorithm is an equivalence class over
implementations, and the interesting scientific object is usually the second one. This course
is a set of techniques for climbing from the first to the second — and for knowing when you
have actually arrived rather than told yourself a story.

It is written for someone with a strong mathematical background who wants the derivations, not
the summaries. Every unit contains real mathematics, runnable models, and exercises with full
worked solutions.

---

## The two threads

The course runs **two interleaved threads**. Read them together.

**Thread A — The conversion.** How do you look at a circuit and extract the algorithm? Six
methodological modules, then ten case studies where the conversion has actually been carried
out (with wildly varying degrees of success), then a synthesis section on when the enterprise
fails and a capstone project.

**Thread B — The structures.** Neural circuits have repeatedly forced the invention of *new
mathematical objects* — not just new applications of old math, but genuinely new invariants,
new correspondences, new formal systems. Convex codes and the neural ideal. Clique topology.
Combinatorial threshold-linear networks. Manifold capacity. Assembly calculus. This thread
catalogues them, shows how each was forced into existence by a circuit that refused to fit the
math already on the shelf, and — most importantly — teaches you to *recognize the pattern* so
you can do it yourself.

Thread B lives in [`structures/`](structures/) and has its own
[index and methodology note](structures/README.md). If you are the kind of person who reads
math textbooks for pleasure, start there and let it pull you into Thread A.

---

## Map

### Phase 0 — Orientation

| Unit | Title |
|---|---|
| [00](00-orientation/README.md) | Marr's levels, multiple realizability, and the identifiability problem |

### Phase 1 — Methods

| Unit | Title | Primitive |
|---|---|---|
| [01](part1-foundations/01-dynamical-systems.md) | Dynamical systems as the lingua franca | Phase portraits; fixed-point reverse-engineering |
| [02](part1-foundations/02-population-geometry.md) | The linear algebra and geometry of populations | Latent manifolds; low-rank connectivity |
| [03](part1-foundations/03-efficient-coding.md) | Efficient coding: deriving circuits from principles | Running the derivation *forward* |
| [04](part1-foundations/04-probabilistic-computation.md) | Circuits that represent posteriors | Sampling vs. parametric codes |
| [05](part1-foundations/05-learning-rules.md) | Learning rules as the extracted algorithm | Rule → implied objective |
| [06](part1-foundations/06-control-and-filtering.md) | Control, filtering, and internal models | "The circuit subtracts a prediction" |

### Phase 2 — Case studies

Sensory coding is weighted heavily, but the point of the breadth is comparative: you learn the
*moves* by seeing the same move made in four different systems.

| Unit | Circuit | Algorithm extracted |
|---|---|---|
| [C1](part2-case-studies/C1-expansion-and-sparsening.md) | PN→KC, cerebellar granule layer | Sparse random projection, LSH, random-feature kernels |
| [C2](part2-case-studies/C2-transient-synchrony.md) ★ | Locust antennal lobe | Reservoir computation / dynamic decorrelation / heteroclinic channels |
| [C3](part2-case-studies/C3-coincidence-detection.md) | Barn owl nucleus laminaris | Cross-correlation |
| [C4](part2-case-studies/C4-continuous-attractors.md) | *Drosophila* ellipsoid body; oculomotor integrator | Integration on an attractor manifold |
| [C5](part2-case-studies/C5-autoassociative-memory.md) | CA3 recurrent collaterals | Content-addressable memory (→ transformer attention) |
| [C6](part2-case-studies/C6-grid-cells.md) | Entorhinal grid modules | Modular error-correcting code *vs.* successor representation |
| [C7](part2-case-studies/C7-dopamine-and-td-learning.md) | Midbrain dopamine, basal ganglia | Temporal-difference learning; actor–critic |
| [C8](part2-case-studies/C8-divisive-normalization.md) | Cortex, antennal lobe | Divisive normalization — *several* valid readings |
| [C9](part2-case-studies/C9-adaptive-filters.md) | Electric fish ELL, cerebellum | Adaptive filter / LMS; learned negative image |
| [C10](part2-case-studies/C10-evidence-accumulation.md) | LIP, parietal cortex | Sequential probability ratio test |

★ = flagship unit.

### Phase 3 — Synthesis

| Unit | Title |
|---|---|
| [S1](part3-synthesis/S1-rnn-as-model-organism.md) | Train an RNN, reverse-engineer it, compare to biology |
| [S2](part3-synthesis/S2-degeneracy-and-limits.md) | When there is no clean algorithm |
| [S3](part3-synthesis/S3-capstone.md) | Capstone: a level-2 specification for transient-synchrony olfactory coding |

### Thread B — Structures

See [`structures/README.md`](structures/README.md) for the full index and the methodology of
recognizing new mathematical objects in the wild.

### Companion code

[`notebooks/`](notebooks/README.md) — three runnable models (C1, C2, C4) with plots, parameter
sweeps, and null models. numpy + matplotlib only.

---

## A suggested 26-week pace

Roughly 6–10 hours/week. The structures notes are short and deliberately offset from the main
track so that you meet an object *after* you have met the circuit that produced it.

| Weeks | Thread A | Thread B |
|---|---|---|
| 1 | Orientation | Methodology note (in `structures/README.md`) |
| 2–3 | 01 Dynamical systems | S-07 Random matrices & chaos |
| 4–5 | 02 Population geometry | S-08 Low-rank connectivity |
| 6–7 | 03 Efficient coding | S-06 Hyperbolic odour space |
| 8 | 04 Probabilistic computation | — |
| 9 | 05 Learning rules | — |
| 10 | 06 Control & filtering | — |
| 11–12 | C1 Expansion & sparsening | S-11 Expanders & optimal degree |
| 13–15 | **C2 Transient synchrony** | S-10 Tropical & piecewise-linear geometry |
| 16 | C3 Coincidence detection · C8 Normalization | — |
| 17–18 | C4 Continuous attractors | S-05 Toroidal topology of grid cells |
| 19 | C5 Autoassociative memory | S-04 Manifold capacity |
| 20 | C6 Grid cells | S-01 Convex codes & the neural ideal |
| 21 | C7 Dopamine & TD · C9 Adaptive filters | S-02 Clique topology |
| 22 | C10 Evidence accumulation | S-03 Combinatorial threshold-linear networks |
| 23 | S1 RNN as model organism | S-09 Assembly calculus |
| 24 | S2 Degeneracy and limits | — |
| 25–26 | **S3 Capstone** | — |

If you want a faster pass: Orientation → 01 → 02 → C2 → C1 → S1 → S2 → S3, plus the whole
structures thread. That is the spine.

---

## How to use this

**Do the exercises.** Every unit has six to ten, tagged ★ (consolidation), ★★ (real work),
★★★ (research-adjacent, occasionally open). All solutions are fully worked, not sketched. The
★★★ problems are where the course actually happens — several of them are live research
questions that I believe are unpublished, and they are marked as such.

Solutions sit in collapsible blocks directly beneath each exercise, so the page reads as a
problem set until you choose otherwise. There are 58 exercises in Foundations alone, and roughly
150 across the course.

**Run the code.** Each case study contains a small numpy model, deliberately small enough to
read in one sitting and modify in an afternoon. A model you have perturbed is worth ten you
have read about. Three of them are expanded into [notebooks](notebooks/README.md) with plots
and null models.

**Keep a conversion log.** For every case study, write one page: (1) the level-1 problem
statement, (2) the algorithm as a mathematical object, (3) the mapping to circuit elements,
(4) the falsifiable prediction, (5) what you would measure to kill it. By C10 you will have ten
of these and the pattern-matching becomes automatic. This is the real deliverable of the course;
the capstone is just the eleventh one, done on your own problem.

**Be adversarial.** Section 6 of every unit ("Where it breaks" / "Open problems") is not
throat-clearing. A large fraction of published circuit→algorithm conversions are, on inspection,
unfalsifiable. [S2](part3-synthesis/S2-degeneracy-and-limits.md) gives you a checklist; apply it
retroactively to every unit, including the ones the course presents approvingly.

---

## Prerequisites

Assumed: linear algebra, probability, ODEs, and enough real analysis to be unbothered by an
argument. Helpful but developed as needed: dynamical systems (bifurcations, normal forms),
information theory, statistical mechanics of disordered systems (the replica method is used and
quoted, not derived from scratch), and elementary algebraic topology (simplicial complexes,
homology) for the structures thread.

Not assumed: any neuroscience. Every case study states the biology from scratch.

---

## Standing bibliography

Three things worth owning rather than borrowing:

- **Dayan & Abbott (2001)**, *Theoretical Neuroscience* — the standard reference; still the best
  single source for mechanics and notation.
- **Strogatz**, *Nonlinear Dynamics and Chaos* — for M1, if the dynamical systems material is
  not already reflexive.
- **Marr (1982)**, *Vision* — read chapter 1 in week 1 and again in week 26. It reads
  differently the second time.

---

## Contributing / errata

This is a living document. Corrections to the mathematics, and especially to the citations, are
welcome via issue or PR.

*Licensed CC BY 4.0 — see [LICENSE](LICENSE).*
