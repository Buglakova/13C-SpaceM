# Single-cell isotope tracing - 13C-SpaceM  

## About this repository

This repository contains the code necessary to reproduce the figures from the manuscript [13C-SpaceM: Spatial single-cell isotope tracing reveals heterogeneity of de novo fatty acid synthesis in cancer](https://www.biorxiv.org/content/10.1101/2023.08.18.553810v2) 

If you are interested in modeling your own bulk/MSI isotope tracing data with a binomial distribution, check out this easier example with more explanations: [13C-binomial](https://github.com/Buglakova/13C_binomial) 

## Installation

```
    git clone https://github.com/Buglakova/13C-SpaceM
    cd 13C-SpaceM
    conda env create -f environment.yml
```

Usage:

```
    conda activate sciso_env
```

### Possible issues
#### It takes too long to solve environment
It's a widespread issue with using `conda`. Consider using `mamba` solver and refer to the official `conda` website.

#### Scanpy errors
This project uses Scanpy's pickled `anndata` objects to store the single-cell data. Sometimes different combinations of versions of `pandas`, `anndata`, `pickle` and `scanpy` can cause errors in reading the files or having imports work correctly. For possible troubleshooting take a look at the dump of the environment used for processing `environment_dump.yml`

## Download data
Data for Hypoxia/Normoxia and ACLYkd experiments can be downloaded from Biostudies.

Data for mouse brain experiments can be downloaded from [Metaspace project](https://metaspace2020.eu/project/buglakova-2024).


## Single-cell workflow

### SpaceM processing
Raw data for single-cell lipidomics consists of microscopy images and MALDI images. `SpaceM` software was used for correlating the modalities and producing single-cell data by overlaying mass spectrometry images with microscopy and attributing ablation marks with segmentation cell masks. This method is described in the corresponding paper: [SpaceM reveals metabolic states of single cells](https://www.nature.com/articles/s41592-021-01198-0).   
`SpaceM` is a closed-sourced software, so for reproducibility we share the output of `SpaceM`.

Normally mass spec peak annotations downloaded from [Metaspace](https://metaspace2020.eu/) are used in `SpaceM` pipeline, but Metaspace only generates isotopologue distributions based on natural isotope abundance, therefore here after generating all possible isotopologues for a given molecule, [pyimzML](https://github.com/alexandrovteam/pyimzML) is used to extract ion images and store in the format expected by `SpaceM`.


Therefore in this project we performed image registration and segmentation with `SpaceM` in a usual way until Metaspace annotation, then used `extract_isotopologue_ion_images.py` script to parse each `imzML` file and produce ion images for all isotopologues of the ions from molecule list.

``` 
    python sciso/extract_isotopologue_ion_images.py analysis/ion_list_fa.csv 20220322_AB_DKFZHypoxia/W2 
```

Example of extracted ion images and resulting objects is in `analysis/00_extract_ion_images.ipynb`. You don't need to download raw `imzML` files and rerun this step, as we share the output of `SpaceM`.

### Single-cell analysis
Although `SpaceM` is used to register together pre-MALDI microscopy and MALDI ion images, signal normalization and single-cell analysis (aggregation of signal from each ablation mark overlapping with a given cell in a single cell's readout) steps are specific to this project, therefore they are implemented outside of the SpaceM pipeline.  
`SpaceM` outputs results of analysis as a directory with the following structure:  
```
    W1
    |-- config.json
    |-- analysis
    |   |-- ablation_mark_analysis
    |   |   |-- spatiomolecular_adata.h5ad
    |   |   |-- spatiomolecular_matrix.csv
    |   |-- cell_segmentation_external
    |   |   |-- cell.labels.tif
    |   |   |-- cell_segmentation.tif
    |   |-- ion_image_registration
    |   |   |-- ablation_mark.labels.tif
    |   |   |-- transformation.json
    |   |-- overlap_analysis2
    |   |   |-- ablation_mark.regions.csv
    |   |   |-- cell.regions.csv
    |   |   |-- overlap.labels.tif
    |   |   |-- overlap.regions.csv
    |   |-- single_cell_analysis
    |   |   |-- spatiomolecular_adata.h5ad
    |   |   |-- spatiomolecular_matrix.csv
    |   |-- transformation
    |   |   |-- ablation_mark_mask_pre_maldi_cropped.tif
    |   |   |-- cropped_pre_maldi_channels.json
    |   |   |-- pre_maldi_crop_frame.json
    |   |   |-- cropped_pre_maldi_channels
    |   |   |   |-- img_t1_z1_c0.tif
    |   |   |   |-- img_t1_z1_c1.tif
    |   |   |   |-- img_t1_z1_c2.tif
    |   |   |-- cropped_post_maldi_channels
    |   |   |   |-- img_t1_z1_c0.tif
    |   |   |   |-- img_t1_z1_c1.tif
    |
    |
    W2
    |--config.json
    |--analysis
    |  |--...
    ...
```
Run single-cell analysis for all wells on each slide:
```
    python sciso/sc_analysis.py $data_dir/20220322_AB_DKFZHypoxia/slide2/spacem_data $data_dir/20220322_AB_DKFZHypoxia/slide2/anndata
    python sciso/sc_analysis.py $data_dir/20220322_AB_DKFZHypoxia/slide3/spacem_data $data_dir/20220322_AB_DKFZHypoxia/slide3/anndata
    python sciso/sc_analysis.py $data_dir/20220411_AB_DKFZACLYac/slide1/spacem_data $data_dir/20220411_AB_DKFZACLYac/slide1/anndata
    python sciso/sc_analysis.py $data_dir/20220411_AB_DKFZACLYac/slide2/spacem_data $data_dir/20220411_AB_DKFZACLYac/slide2/anndata
    python sciso/sc_analysis.py $data_dir/20220411_AB_DKFZACLYac/slide3/spacem_data $data_dir/20220411_AB_DKFZACLYac/slide3/anndata
```
This script performs the following steps:
- Natural isotope abundance correction for each fatty acid in each ablation mark
- Averaging of signal in each cell by taking a median for each annotated m/z value
- Normalizing isotopologue distributions for each fatty acid for each cell 

The results are saved in an `anndata` object joined for all wells on the slide. This analysis is done per slide to keep the file sizes manageable. 

### Quality control and data aggregation
After single-cell analysis data from all slides is aggregated together and quality control is performed, filtering out cells by size and number and total intensity of detected ions. Subsequently the cells are classified using normalized GFP intensity. Running all notebooks in `analysis/hypoxia` and `analysis/ACLYkd` produces `anndata` objects that contain all information per cell that is further used to generate manuscript plots. We share these files as well.

### Reproducing manuscript figures
`manuscript_figures` contains notebooks that allow to reproduce plots for each figure panel based on the `anndata` objects produced in the previous step and SpaceM output (for registered microscopy images and segmentation). Figures 2-4 contain plots based on single-cell data.

## Mouse brain tissue image analysis
For the tissue imaging analysis the analysis was done per pixel of the MALDI image instead of per cell.
In each pixel we do the same analysis steps as for the single-cell data:
- Perform natural isotope abundance correction for each ion's isotopologue distribution
- Fit binomial distribution to each fatty acid's isotopologue distribution

To reproduce Figure 5 (panels c and d) and Figure 6 (panel a) do the following:
- Download `imzml` for the dataset `2023-09-25_ME_PattiLab_MouseBrain_4-4_010G_AIF_600-1600_30NCE_100-600_10umss_285x230_32at` from the [Metaspace project](https://metaspace2020.eu/project/buglakova-2024)
- Run the script   
  ```python analysis/brain/00_extract_ion_images.py $data_dir/raw/2023-09-25_ME_PattiLab_MouseBrain_4-4_AIF_600-1600_30NCE_100-600_10umss_285x230_32at.imzML  analysis/brain/ion_list_fa.csv output_dir```
`