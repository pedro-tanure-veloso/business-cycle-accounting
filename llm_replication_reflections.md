---
title: "Reflections on Vibe Coding the Replication"
topic: "methodology"
layer: "all"
status: "reference"
last_updated: "2026-05-06"
---

# Reflections on Vibe Coding the Replication 
This document synthesizes observations regarding the replication of quantitative research—specifically the conversion of Matlab code to Python—using large language models (LLMs).

## 1. The Research Intuition Gap
There is a significant divergence between an LLM's technical proficiency and a human researcher's procedural intuition. While the model exhibits high-level knowledge, it lacks the instinctive drive to treat the original source as a "source of truth."

* **Absence of Benchmarking:** A human researcher naturally compares outputs at every margin: data construction, individual function outputs, and final tables. The LLM tends to create a "surface-level interpretation" rather than a rigorous audit.
* **Source Neglect:** Despite having access to original `.m` files, the model frequently defaults to guessing logic. It fails to "stop and inspect" when results diverge, a step that is second nature to an experienced researcher.

## 2. Context Decay and Memory Management
LLMs demonstrate a "shorter memory" than their context windows suggest during long, complex workflows.

* **Regression Errors:** The model often reverts to previously discarded solutions or re-introduces bugs it had already solved. This occurs because the "win" or the specific solution eventually drifts out of its immediate focus.
* **The "Diary" Protocol:** To counter context decay, it is essential to force the model to maintain a **running log of findings**, failed attempts, and core logic. Frequently prompting the model to re-read this log prevents it from cycling back to inefficient or incorrect methods.

## 3. Workflow Constraints & Input Quality
The marginal utility of the tool is heavily dependent on the quality of the provided environment and the specificity of the guidance.

* **The "Messy Folder" Tax:** Providing a non-machine-readable or disorganized directory (e.g., a cluttered Dropbox folder) overwhelms the model, leading to degraded synthesis and increased errors.
* **The RA Analogy:** The LLM is best viewed as a **Research Assistant with a PhD-level understanding of theory but an undergraduate-level grasp of research methodology.** It knows *how* to write code but does not instinctively know *when* or *why* to validate it against the baseline.

## 4. Strategic Recommendations for Serious Replication
For a successful replication, the researcher must provide the procedural structure the model lacks:

* **Enforce Marginal Benchmarking:** Instruct the model to replicate numbers from the start. Verify that Function A in Python produces the exact output of Function A in Matlab before proceeding to the next step.
* **Documentation as a Context Anchor:** Require the model to document variable definitions and function logic one by one. This keeps the logic within the active context and minimizes "hallucinated" model versions.
* **Intentionality Over Automation:** Optimal results are achieved where marginal benefit equals marginal cost; the "cost" of providing heavy upfront guidance on process is rewarded by significantly higher accuracy in the final output.
