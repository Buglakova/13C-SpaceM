This directory contains data necessary to reproduce the analysis of Hypoxia/Normoxia experiment.

hypoxia_metadata.csv: list of datasets
IsoCorrectoR data.csv: bulk mass spectrometry data for samples with the same conditions
anndata - Scanpy objects which contain the final isotope tracing data:
    - hypoxia_adata_am_concat.h5ad: normalized and corrected for natural isotope abundance data per pixel of MALDI image
    - hypoxia_adata_concat_fit.h5ad: normalized and corrected for natural isotope abundance data per cell
slide2 and slide3: raw data in the form of SpaceM output. Each separate directory corresponds to one dataset (one MALDI image + one microscopy image) and has the same structure:

    W1
    |-- config.json: pipeline settings and info about the dataset
    |-- analysis
    |   |-- ablation_mark_analysis
    |   |   |-- spatiomolecular_adata.h5ad: Scanpy Anndata object with isotope tracing data per pixel of MALDI image
    |   |   |-- spatiomolecular_matrix.csv: same but as csv
    |   |-- cell_segmentation_external
    |   |   |-- cell.labels.tif: cell segmentation
    |   |-- ion_image_registration
    |   |   |-- ablation_mark.labels.tif: ablation mark segmentation registered with the microscopy
    |   |   |-- transformation.json: information about the transform needed to register AM segmentation image with microscopy
    |   |-- overlap_analysis2
    |   |   |-- ablation_mark.regions.csv: list of AM masks and their bounding boxes
    |   |   |-- cell.regions.csv: list of cell masks and their bounding boxes
    |   |   |-- overlap.labels.tif: overlap of AM segmentation masks with cell segmentation
    |   |   |-- overlap.regions.csv: list overlap areas of AMs with cell masks
    |   |-- single_cell_analysis
    |   |   |-- spatiomolecular_adata.h5ad: Scanpy Anndata object with isotope tracing data and morphological features per cell
    |   |   |-- spatiomolecular_matrix.csv: same but as csv
    |   |-- transformation
    |   |   |-- ablation_mark_mask_pre_maldi_cropped.tif: ablation mark segmentation (AM numbering may not be consistent with the microscopy orientation)
    |   |   |-- cropped_pre_maldi_channels.json: paths to premaldi microscopy files
    |   |   |-- pre_maldi_crop_frame.json: registration metadata
    |   |   |-- cropped_pre_maldi_channels
    |   |   |   |-- img_t1_z1_c0.tif: transmission
    |   |   |   |-- img_t1_z1_c1.tif: DAPI
    |   |   |   |-- img_t1_z1_c2.tif: GFP
    |   |   |   |-- img_t1_z1_c2_registered.tif: GFP before desiccation, was used for the analysis
    |   |   |-- cropped_post_maldi_channels
    |   |   |   |-- img_t1_z1_c0.tif: transmission
    |   |   |   |-- img_t1_z1_c1.tif: DAPI
