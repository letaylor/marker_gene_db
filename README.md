
# Description

This repository contains list of various marker genes.

# Description

```bash
<marker_group_id>/README.md
<marker_group_id>/database.tsv
```

* Each marker gene set is grouped by a `marker_group_id` (e.g., `hypoxia`, `cellular_stress`, `islet_celltypes`).
* `<marker_group_id>/README.md`: description of this marker_group_id.
* `<marker_group_id>/database.tsv`: the actual database file.


NOTE: There may be duplicate hgnc_symbol in the case where a gene id mapped to multiple ensemble gene ids.
Each file contains at least following columns, in addition to any other columns one may add:


Column | Description
--- | ---
hgnc_symbol | HUGO gene nomenclature (HGNC) gene symbol.
grouping_id | Additional grouping of (example: 'immune' to group all immune cells together or `hypoxia` `beta_cell`)
study_id | Unique identifier of the reference study. `<first_author_last_name><first_author_initials>-<PubMed_ID>`.
description | Any additional notes.
ensembl_gene_id | Ensembl gene id.
ensembl_version | Ensembl database version. [Optional]
biopsy_tissue | Tissue biopsied for the single cell experiment. [Optional]
biopsy_cell_selection | Cell type selection applied to the tissue biopsy (example: 'lamina propria' if selected lamina propria from colon biopsies). [Optional]
notes | Any additional notes. [Optoinal]
