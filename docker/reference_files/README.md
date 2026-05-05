Docker reference files
======================

The Docker image copies this directory to `/home/rstudio/reference_files/` so
hCocena workflows can use the standard supplementary-data path immediately.

Bundled files include MSigDB v2023.1 human and mouse gene-set GMT files for
Hallmark, Reactome, KEGG, and GO collections, plus `TFcat.txt` and
`immune_sig_m.csv` helper references used by hCocena workflows.

These files are bundled for reproducible Docker runs. Users remain responsible
for complying with the license and attribution terms of the upstream reference
databases used in their analyses.

Upstream reference pages:

- MSigDB collections: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
- MSigDB license terms: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp
- KEGG legal information: https://www.kegg.jp/kegg/legal.html
