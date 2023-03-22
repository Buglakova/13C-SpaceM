# Single cell isotope tracing - 13C-SpaceM  

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
## SpaceM processing  
Normally annotations downloaded via [Metaspace](https://metaspace2020.eu/) are used in SpaceM pipeline, but Metaspace only generates isotopologue distributions based on natural isotope abundance, therefore here after generating all possible isotopologues for a given molecule, [pyimzML](https://github.com/alexandrovteam/pyimzML) is used to extract ion images and store in the format expected by SpaceM (closed-source).


Process dataset with SpaceM in a usual way until Metaspace annotation, then use `extract_isotopologue_ion_images.py` to parse imzML file and produce ion images for all isotopologues of the ions from molecule list.

``` 
    python sciso/extract_isotopologue_ion_images.py analysis/ion_list_fa.csv 20220322_AB_DKFZHypoxia/W2 
```

Example of extracted ion images and resulting objects is in `analysis/00_extract_ion_images.ipynb`

## Single-cell analysis
SpaceM is used to register together pre-MALDI microscopy and MALDI ion images, but normalization and single-cell analysis steps are specific to this project, therefore reimplemented outside of the SpaceM pipeline.  
SpaceM outputs results of analysis as a directory with the following structure (only relevant files shown):  
```
    W1
    |-- analysis
    |   |-- ablation_mark_analysis
    |   |   |-- spatiomolecular_adata.h5ad
    |   |-- metaspace
    |   |   |-- ion_images.pickle
    |   |   |-- ion_images_metadata.csv
    |   |-- overlap_analysis2
    |   |   |-- overlap.regions.csv
    |   |-- single_cell_analysis
    |   |   |-- spatiomolecular_adata.h5ad
    |   |-- transformation
    |   |   |-- cropped_pre_maldi_channels
    |   |   |   |--
    |   |   |   |--
    |-- config.json
```
Run single-cell analysis for all wells on slide:
```
    python sciso/sc_analysis.py 20220322_AB_DKFZHypoxia/slide3/spacem_data 20220322_AB_DKFZHypoxia/slide3/anndata
```