---
title: Case studies
nav_order: 4
has_children: true
---

# Phase 2 — Case studies

Ten circuits where someone has attempted the conversion. They are not equally successful, and
the variation is the point: C7 (dopamine → TD) is about as complete as neuroscience gets, while
C2 (transient synchrony) has four live competing readings and no consensus. Learning to tell
those situations apart is more valuable than learning either one.

Sensory coding is weighted heavily. The breadth is there for comparison — you learn a *move* by
watching it made in systems that share no biology.

| Unit | Circuit | Algorithm extracted | Status |
|---|---|---|---|
| [C1](C1-expansion-and-sparsening.md) | PN→KC; cerebellar granule layer | Sparse random projection, LSH, random-feature kernels | Well established |
| [C2](C2-transient-synchrony.md) ★ | Locust antennal lobe | Reservoir / decorrelation / heteroclinic channels | **Contested** |
| [C3](C3-coincidence-detection.md) | Barn owl nucleus laminaris | Cross-correlation | Established for owl; contested for mammals |
| [C4](C4-continuous-attractors.md) | *Drosophila* ellipsoid body; oculomotor integrator | Integration on an attractor manifold | Well established |
| [C5](C5-autoassociative-memory.md) | CA3 recurrent collaterals | Content-addressable memory | Established as theory; contested *in vivo* |
| [C6](C6-grid-cells.md) | Entorhinal grid modules | Modular error-correcting code *vs.* successor representation | **Two live readings** |
| [C7](C7-dopamine-and-td-learning.md) | Midbrain dopamine; basal ganglia | Temporal-difference learning; actor–critic | The field's crown jewel |
| [C8](C8-divisive-normalization.md) | Cortex; antennal lobe | Divisive normalization | Motif certain, *reading* non-unique |
| [C9](C9-adaptive-filters.md) | Electric fish ELL; cerebellum | Adaptive filter / LMS; learned negative image | Well established |
| [C10](C10-evidence-accumulation.md) | LIP; parietal cortex | Sequential probability ratio test | Established as behaviour; causality contested |

★ = flagship unit, and the one the [capstone](../part3-synthesis/S3-capstone.md) builds on.

---

**Keep a conversion log.** One page per case study: (1) the level-1 problem statement, (2) the
algorithm as a mathematical object, (3) the mapping to circuit elements, (4) the falsifiable
prediction, (5) what you would measure to kill it. Ten of these is the real output of this phase.

**Read the Status column adversarially.** Where it says "well established," ask what would have
had to happen for the field to notice it was wrong. Sometimes the answer is uncomfortable —
which is what [S2](../part3-synthesis/S2-degeneracy-and-limits.md) is for.
