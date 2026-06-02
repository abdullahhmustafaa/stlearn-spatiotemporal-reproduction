# Datasets

All datasets used in this reproduction are **publicly available** and load automatically inside the Colab notebook via `scanpy.datasets.visium_sge()`.

## Datasets Used

| Dataset | Technology | Spots | Genes | Used for |
|---------|-----------|-------|-------|---------|
| Mouse Brain Sagittal Anterior | 10x Visium | 2702 | ~20,000 | stSME clustering + PSTS trajectory |
| Human Breast Cancer Block A Section 1 | 10x Visium | 3813 | ~33,000 | SCTP cell-cell interaction |

## Manual Download Links (if needed)

If the automatic download fails in the notebook, you can manually download from:

- **Mouse Brain:** https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Mouse_Brain_Sagittal_Anterior
- **Breast Cancer:** https://support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_1

## Other datasets mentioned in the paper (not reproduced here)

| Dataset | Access |
|---------|--------|
| Mouse TBI Visium (in-house) | GEO: GSE236171 |
| Mouse embryo sci-Space | GEO: GSE166692 |
| SeqFISH+ mouse cortex | https://github.com/CaiGroup/seqFISH-PLUS |
| Slide-seq mouse hippocampus | Dropbox link in paper |
