# Cover letter for Nature Cancer (DeSurv), Technical Report

<!-- REMAINING AUTHOR ITEMS BEFORE SENDING:
     1. Referee contact emails are entered directly in the submission system.
     2. Confirm the editor name (Julieta Alfonso) against the presubmission thread.
     3. Confirm no suggested referee is a recent co-author/collaborator, shares a
        grant, or shares an institution.
     4. Submission form: article type = Technical Report (not Article). The
        category must match the first sentence of this letter.
     Render: pandoc -f markdown+hard_line_breaks -V geometry:margin=1in
     (hard_line_breaks keeps the address and signature blocks on separate lines;
     the 1in margin keeps the letter to two pages). Lab style: no em-dashes,
     direct voice. -->

7 September 2026

Dr. Julieta Alfonso
Senior Editor, Nature Cancer

Dear Dr. Alfonso,

We are submitting a Technical Report, "DeSurv identifies a replicable tumor-stroma prognostic architecture in pancreatic cancer through survival-supervised matrix factorization," for consideration at Nature Cancer. Thank you for inviting a full submission following our presubmission enquiry. The underlying study is unchanged; the revised title names the method directly, and we believe that a Technical Report best reflects the nature of its contribution: a new survival-supervised matrix-factorization framework for discovering and transferring clinically relevant gene programs, with pancreatic cancer providing the biological and clinical validation setting.

Molecular subtypes of pancreatic cancer are typically found by factorizing bulk expression to recover programs that best reconstruct the transcriptome, testing those programs against outcome, and validating retrospectively. DeSurv instead makes survival information part of discovery: the survival gradient acts on the gene weights that define each program, while the reconstruction objective keeps those programs representing the observed transcriptome. Because supervision acts on the gene weights, the trained basis is frozen and a new tumor is scored directly against the fixed gene-weight matrix, with no refitting and no use of validation outcomes, so cross-study transportability becomes a direct test of the same fixed representation. The related CoxNMF and SurvNMF formulations place the survival term on the sample-side latent representation; the manuscript tabulates that distinction.

In pancreatic cancer, DeSurv recovered a three-program tumor-stroma architecture whose dominant survival-aligned program captured coordinated classical malignant and tumor-restraining (restCAF-like) stromal features. That program was recovered de novo without subtype labels and retained prognostic information beyond the established PurIST and DeCAF classifiers, so it is not reducible to them, while an unsupervised control matched on rank, optimizer, consensus-initialization procedure and tuning budget did not recover it, supporting a role for survival supervision beyond tuning effort alone in producing the reorganization. In simulations with known truth, supervision improved recovery of the genes defining the prognostic program when those genes contributed relatively little transcriptomic variation. The frozen programs then transferred without refitting across five external datasets spanning four independent patient cohorts, and the principal survival-aligned direction was independently supported by alternative supervised methods, GATA6 RNA in situ hybridization, and exploratory analysis in treated metastatic disease. The DeSurv R package and the reproducibility repository that renders the manuscript from archived results are released on GitHub and archived on Zenodo.

The advance is not a claim of superior discrimination. Higher-rank unsupervised factorization and strong supervised predictors reach comparable external discrimination. Instead, DeSurv recovers the dominant survival-aligned signal within a compact three-program representation that separately preserves tumor and stromal programs and can be scored unchanged across cohorts. The contribution is therefore in the biological representation learned and transported by the method, rather than in predictive performance alone.

Because the central claims concern the factorization itself, we ask that at least one referee have expertise in matrix factorization or related latent-variable models. We confirm that this manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript and agree with its submission to Nature Cancer. Competing interests are disclosed in the manuscript: N.U.R. and J.J.Y. are inventors on two pancreatic-cancer subtype-classification patents; the PurIST and DeCAF classifiers are used only as published, independent comparators and are not components of DeSurv.

We suggest the following referees:

- Elana J. Fertig (University of Maryland School of Medicine), computational cancer biology and non-negative matrix factorization of tumor transcriptomes; a natural methodological reviewer for interpretable, outcome-aligned latent-program models.
- Andrew J. Aguirre (Dana-Farber Cancer Institute), pancreatic cancer functional genomics, molecular subtypes, and precision oncology.
- Anguraj Sadanandam (The Institute of Cancer Research, London), transcriptomic molecular subtyping and classification across pancreatic and gastrointestinal cancers.
- Peter Bailey (Botton-Champalimaud Pancreatic Cancer Centre, Lisbon), pancreatic cancer molecular subtypes and tumor-microenvironment transcriptomics.

We have no strong objections to particular referees and defer to the editors on exclusions; the suggested reviewers above are, to our knowledge, independent and unconflicted.

Thank you for considering our work.

Sincerely,

Naim U. Rashid, PhD
On behalf of all authors
Department of Biostatistics, University of North Carolina at Chapel Hill
nur2@email.unc.edu
