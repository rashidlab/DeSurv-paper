# Cover letter for Nature Cancer (DeSurv), Technical Report

<!-- REMAINING AUTHOR ITEMS BEFORE SENDING:
     1. [Date] below.
     2. Referee contact emails (line marked "contact details to be added").
     3. Confirm the editor name (Julieta Alfonso) against the presubmission thread.
     4. Confirm no suggested referee is a recent co-author/collaborator, shares a
        grant, or shares an institution.
     5. Submission form: article type = Technical Report (not Article). The
        category must match the first sentence of this letter.
     Render: pandoc -f markdown+hard_line_breaks (keeps the address and signature
     blocks on separate lines). Lab style: no em-dashes, direct voice. -->

[Date]

Dr. Julieta Alfonso
Senior Editor, Nature Cancer

Dear Dr. Alfonso,

We are submitting a Technical Report, "DeSurv, a survival-supervised matrix factorization, identifies a replicable tumor-stroma prognostic architecture in pancreatic cancer," for consideration at Nature Cancer. We submit this work as a Technical Report because its principal contribution is a new survival-supervised matrix-factorization framework for discovering and transferring clinically relevant gene programs, with pancreatic cancer providing the biological and clinical validation setting. Thank you for inviting a full submission following our presubmission enquiry. The underlying study is unchanged; the revised title names the method directly, and we believe Technical Report best reflects the nature of its contribution.

The report makes three claims.

First, DeSurv is a new methodological framework. Molecular subtypes of pancreatic cancer are typically found by factorizing bulk expression to recover the axes of greatest variation, testing those axes against outcome, and validating retrospectively. The programs that best reconstruct a transcriptome are not necessarily the ones that carry prognosis, and reconstruction-based factorization is non-identifiable, so which programs count as clinically relevant is settled only after discovery, by a second outcome-dependent selection step. DeSurv instead makes clinical relevance part of discovery: the survival gradient acts on the gene weights that define each program, while the reconstruction objective keeps those programs representing the observed transcriptome. The related CoxNMF and SurvNMF formulations place the survival term on the sample-side latent representation; DeSurv places it on the gene-program basis, and the manuscript tabulates that distinction against the related methods.

Second, that design has a practical consequence. Because supervision acts on the gene weights, the trained basis is frozen after training and a new tumor is scored by projection alone, with no refitting and no use of validation outcomes. The programs discovered in training are therefore the ones tested everywhere else, and cross-study transportability becomes a direct test of the same fixed object rather than of a re-estimated one.

Third, the result is not an artifact of one dataset, one cohort or one fit. In pancreatic cancer DeSurv recovered a three-program tumor-stroma architecture whose dominant survival-aligned program joins classical malignant identity to tumor-restraining (restCAF-like) stromal biology, recovered de novo without subtype labels and not reducible to the established PurIST and DeCAF classifiers, which jointly explain roughly a third of its variance while it remains prognostic after adjustment for both. An unsupervised model given the same rank, optimizer, consensus initialization and tuning budget did not recover that program, so the reorganization is attributable to supervision rather than to computational effort. In simulations with known ground truth, supervision recovered the genes defining the prognostic program more accurately, most when that program contributed little to overall variance. Independently seeded replicates of the full procedure recover the reported decomposition closely, whereas single runs of the underlying optimizer do not, which is what supports treating it as a fixed object. The frozen programs transferred without refitting across five external validation datasets spanning four independent patient cohorts; the survival-aligned program was recovered independently by three other survival-supervised methods; its classical tumor identity was corroborated by GATA6 RNA in situ hybridization; and in an exploratory, arm-stratified analysis it remained associated with overall survival in treated metastatic disease. The DeSurv R package and the reproducibility repository that renders the manuscript from archived results are released on GitHub and archived on Zenodo.

We want to be explicit about what is not claimed. DeSurv does not predict survival better than the alternatives. Higher-rank unsupervised factorization and strong penalized models reach comparable external discrimination; the prognostic signal is accessible to them. Unsupervised factorization needs more than twice as many factors to approach the same cross-validated prognostic information, and at higher rank it broadens coverage of established programs while approaching that prognostic information. Penalized models return a risk score without the decomposed tumor and stromal programs. The contribution is in the biological representation the method learns and can carry unchanged across cohorts, not in superior discrimination.

Because the central claims concern the factorization itself, we ask that at least one referee have expertise in matrix factorization or related latent-variable models. Among the suggestions below, Dr. Fertig has that expertise.

This manuscript is not under consideration elsewhere and all authors have approved the submission. Competing interests are disclosed in the manuscript: N.U.R. and J.J.Y. are inventors on two pancreatic-cancer subtype-classification patents; the PurIST and DeCAF classifiers are used only as published, independent comparators and are not components of DeSurv.

We suggest the following knowledgeable and, to our knowledge, unconflicted referees (contact details to be added):

- Elana J. Fertig (Johns Hopkins University), computational cancer biology and non-negative matrix factorization of tumor transcriptomes; a natural methodological reviewer for interpretable, outcome-aligned latent-program models.
- Andrew J. Aguirre (Dana-Farber Cancer Institute), pancreatic cancer functional genomics, molecular subtypes, and precision oncology.
- Anguraj Sadanandam (Institute of Cancer Research, London), transcriptomic molecular subtyping and classification across pancreatic and gastrointestinal cancers.
- Peter Bailey (University of Glasgow), pancreatic cancer molecular subtypes and tumor-microenvironment transcriptomics.

We have no strong objections to particular referees and defer to the editors on any exclusions, noting only that, per journal policy, current or recent collaborators and members of our own institutions should not be assigned.

Thank you for considering our work. We look forward to the reviewers' assessment.

Sincerely,

Naim U. Rashid, PhD
On behalf of all authors
Department of Biostatistics, University of North Carolina at Chapel Hill
nur2@email.unc.edu
