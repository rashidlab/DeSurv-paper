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

We are submitting a Technical Report, "DeSurv, a survival-supervised matrix factorization, identifies a replicable tumor-stroma prognostic architecture in pancreatic cancer," for consideration at Nature Cancer. Thank you for inviting a full submission following our presubmission enquiry. The underlying study is unchanged; the revised title names the method directly, and we believe that a Technical Report best reflects the nature of its contribution: a new survival-supervised matrix-factorization framework for discovering and transferring clinically relevant gene programs, with pancreatic cancer providing the biological and clinical validation setting.

Molecular subtypes of pancreatic cancer are typically found by factorizing bulk expression to recover programs that best reconstruct the transcriptome, testing those programs against outcome, and validating retrospectively. DeSurv instead makes survival information part of discovery: the survival gradient acts on the gene weights that define each program, while the reconstruction objective keeps those programs representing the observed transcriptome. Because supervision acts on the gene weights, the trained basis is frozen and a new tumor is scored by projection alone, with no refitting and no use of validation outcomes, so cross-study transportability becomes a direct test of the same fixed representation. The related CoxNMF and SurvNMF formulations place the survival term on the sample-side latent representation; the manuscript tabulates that distinction.

In pancreatic cancer DeSurv recovered a three-program tumor-stroma architecture whose dominant survival-aligned program captures covariation of classical malignant and tumor-restraining (restCAF-like) stromal features. This program was recovered de novo without subtype labels and is not reducible to the established PurIST and DeCAF classifiers. Together, those classifiers explain roughly a third of its variance, yet the program remains prognostic after adjustment for both. An unsupervised control matched on rank, optimizer, consensus-initialization procedure and tuning budget did not recover that program, supporting survival supervision, rather than computational effort alone, as the source of the reorganization. In simulations with known ground truth, supervision recovered the genes defining the prognostic program more accurately when the factor-specific prognostic signal contributed little transcriptomic variation relative to the shared outcome-neutral background. Independently seeded replicates of the full procedure closely reproduced the reported decomposition. The frozen programs transferred without refitting across five external validation datasets spanning four independent patient cohorts; the survival-aligned program was also recovered by three alternative survival-supervised methods; its classical tumor component was corroborated by GATA6 RNA in situ hybridization; and in an exploratory, arm-stratified analysis it remained associated with overall survival in treated metastatic disease. The DeSurv R package and the reproducibility repository that renders the manuscript from archived results are released on GitHub and archived on Zenodo.

We want to be explicit about what is not claimed. DeSurv does not predict survival better than the alternatives. Higher-rank unsupervised factorization and strong penalized models reach comparable external discrimination. In the PDAC training data, the unsupervised control required more than twice as many factors to approach the same cross-validated prognostic information; at the higher rank, it also broadened coverage of established programs. Conventional supervised predictors recovered the prognostic direction without separately resolving the stromal programs. The contribution is in the biological representation learned by the method and transferred unchanged across cohorts, not in superior discrimination.

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
