# Bioconductor Submission Notes

Current package status:

- Package version: `0.99.0`
- Default branch content: package source only
- Full repository tooling: `full-repo` branch
- Local `R CMD check --no-manual`: OK
- Local `BiocCheck`: 0 errors, 0 warnings, notes only

Official references:

- Submission guide: https://contributions.bioconductor.org/bioconductor-package-submissions.html
- DESCRIPTION guidance: https://contributions.bioconductor.org/description.html
- Version numbering: https://contributions.bioconductor.org/versionnum.html
- Branch rename guidance: https://contributions.bioconductor.org/branch-rename-faqs.html

Recommended final repository setup:

1. Ensure the GitHub default branch is `devel`.
2. Keep the `devel` branch package-only.
3. Keep Docker, GitHub Actions, and other auxiliary tooling on `full-repo`.

Manual GitHub step still required:

- Open `https://github.com/BioCompNet/hcocena/branches`
- Rename or switch the default branch to `devel`

Suggested submission issue text:

```text
Package submission

Confirm the following by editing each check box to '[x]':

- [x] I understand that Bioconductor package submissions are reviewed for
  technical quality and project fit.
- [x] The package is available at:
      https://github.com/BioCompNet/hcocena
- [x] The package builds and checks locally with R CMD build / R CMD check.
- [x] I have run BiocCheck and addressed all errors and warnings.
- [x] The DESCRIPTION file uses version 0.99.0 for this initial submission.
- [x] I am the maintainer listed in DESCRIPTION and I am registered on the
      Bioconductor support site.
- [x] I have added the package name to Watched Tags on the support site.

Package name: hcocena
Repository: https://github.com/BioCompNet/hcocena
Maintainer: Thomas Ulas <t.ulas@uni-bonn.de>

Short description:
hcocena provides a network-centric workflow for the horizontal integration and
analysis of transcriptomics datasets. The package combines a modern S4 API based
on MultiAssayExperiment with compatibility for the legacy hcobject workflow and
supports clustering, integrated network analysis, functional enrichment,
upstream inference, cell-type annotation, and longitudinal module analyses.

Bioconductor classes used:
- MultiAssayExperiment
- SummarizedExperiment
- S4Vectors

Availability:
- Software package
- Source package only on the default branch
- Extended development tooling is maintained separately from the submission
  branch
```
