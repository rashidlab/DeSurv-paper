# Jen Jen Yeh — Written Comments on paper_jenjen_comments.pdf

## Page 1 (Abstract, left margin)
"isn't this true by definition since using survival to train in the 1st place" — questioning the abstract's claim that DeSurv yields "clearer survival separation" than unsupervised counterparts. This is a significant conceptual challenge: if the model is trained on survival, better survival separation in training/validation may be expected by construction rather than representing a genuine methodological advance.

## Page 1 (Abstract, right margin)
Notes inconsistent nomenclature for basal/classical terminology in the abstract.

---

## Page 2 (Introduction/Results)
- **Top right margin**: "important to note that these are treatment naive cohorts so at least homogenous in that aspect → also a future limitation" — flagging that the training cohorts are treatment-naive, which is both a strength (homogeneity) and a limitation to note explicitly.
- **Right column**: Asks if all the samples are nonmetastatic.
- **Highlights "more reliably"** in the abstract carry-over text — likely connecting to the same circular reasoning concern from page 1: is it surprising or noteworthy that a survival-supervised method more reliably recovers prognostic programs?
- **Left margin** (near "while the exocrine-compositional variation that dominates standard NMF"): Asks whether this is a citation or their own finding — needs clarification on whether this claim is supported by a reference or is a result being reported here.
- **Left margin** (near iCAF-associated stroma): Asks whether "iCAF-associated" is referring to Elyada iCAF or deCAF proCAF — these are distinct CAF subtypes and the terminology needs to be precise and consistent.

---

## Page 3 (Factor structure section)
- **Large left-margin annotation**: A major methodological critique — essentially asking *"isn't this handicapping NMF by limiting it to k=3? Would standard NMF perform better when not limited to k=3?"* This challenges the fairness of the head-to-head comparison at matched rank.
- **Left margin** (near iCAF/restCAF): Asks whether the paper is claiming that iCAF and restCAF are the same thing — a significant biological concern given that a lot of work has gone into proving they are *not* the same. The paper needs to be careful not to conflate these two distinct CAF subtypes.

---

## Page 4 (Validation + Discussion)
- **Top left**: "also treatment naive?" — asking whether the validation cohorts are also treatment naive (extending the concern from page 2).
- **Beginning of Discussion**: Asks to remove the phrase "rather than assessing them retrospectively" from the opening sentence ("Incorporating survival outcomes directly into NMF factorization rather than assessing them retrospectively reorganizes the transcriptional landscape...").
- **Left margin**: Wants a comparison to previously published hazard ratios before publication — does not specify which published HRs to compare against.

---

## Page 5 (Discussion/Methods)
- **Top margin** (spanning the page): "because DeSurv was trained on the same data as the D1-3 may have meaning. The interesting test would be to use DeSurv on later patients + see what programs it recovers if it can be validated." — asking us to apply DeSurv to newer datasets that include treated patients, since the training cohorts are treatment-naive. This connects to her repeated concern about the homogeneity/limitations of the current cohorts.
- **Left margin**: Extended note encouraging discussion of DeSurv's ability to learn biologically meaningful programs and make predictions; emphasizes that DeSurv is finding biological programs (not just statistical artifacts).

---

## Page 10 (Figure 3 — Heatmaps)
- **Top of page**: "Would use consistent nomenclature throughout — basal vs basal-like, Classical Tumor B ('T' or 't') vs Tumor Classical" — a clear request for uniform terminology across the paper.
- **Circled items on heatmap (Fig 3A)**: Added "- like?" next to "SCISSORS: iCAF" and "SCISSORS: myCAF" — questioning whether those labels are correct/precise for those gene sets.
- **Right margin**: "is this nomenclature in Puleo?" — asking whether the naming conventions used for gene programs actually come from the Puleo reference.

---

## Overall Themes
1. **Circular reasoning** — repeatedly questions whether DeSurv's advantages (clearer survival separation, more reliable recovery) are genuine or expected by construction given that survival is used during training; needs to be addressed head-on in the abstract and introduction.
2. **CAF subtype biology** — concerned about conflation of iCAF, restCAF, and deCAF proCAF, which are distinct subtypes; the paper needs precise, consistent usage throughout, especially given her expertise in this area.
3. **Nomenclature consistency** — wants uniform naming of basal/classical and CAF subtypes throughout the paper and figures.
4. **Cohort homogeneity and limitations** — repeatedly flags that training cohorts are treatment-naive and possibly all nonmetastatic; wants this acknowledged as both a strength and a limitation, and suggests applying DeSurv to newer datasets with treated patients as a future direction.
5. **Fair comparison to NMF** — concerned that restricting NMF to k=3 handicaps it; the comparison may not be fair without also showing NMF at its optimal rank.
6. **Missing benchmarks** — wants comparison to previously published hazard ratios (unspecified) before publication.
7. **Citations/sourcing** — questions whether the claim about exocrine-compositional variation dominating standard NMF is a citation or their own finding.
