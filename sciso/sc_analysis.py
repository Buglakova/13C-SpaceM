import numpy as np
import isocor
import scanpy as sc
import pandas as pd
import argparse
from pathlib import Path
import json
import time
from tqdm import tqdm


def apply_isocorr(ion: str, X: np.array):
    """
    Use isocor to correct isotopologue distributions of certain ion for isotopologue abundances.

    Args:
        ion: formula of the molecule, for which correction is done
        X: 
            matrix of raw isotopologue abundances,
            each row is a separate measurement,
            number of columns should be equal to number of tracer atoms + 1.

    Returns:
        Tuple of np.arrays with the same shape as X:
            X_corr_abs_ion - corrected distributions with the same sum of all peaks as before correction
            X_corr_norm_ion - normalized corrected distributions, sum of all peak intensities == 1
            X_isocor_res - residuals of isocor fitfor each peak
    """
    corrector_HR = isocor.mscorrectors.MetaboliteCorrectorFactory(
        ion,
        tracer="13C",
        resolution=140000,
        mz_of_resolution=200,
        charge=-1,
        correct_NA_tracer=True,
    )
    X_corr_abs_ion = np.zeros(X.shape)
    X_corr_norm_ion = np.zeros(X.shape)
    X_isocor_res = np.zeros(X.shape)
    for idx, row in enumerate(X):
        if idx % 1000 == 0:
            print("%d out of %d" % (idx, X.shape[0]))
            
        # Not do calculations for AMs only with 0s
        if np.all(row == 0):
            corrected_area = np.zeros_like(row)
            iso_fraction = np.zeros_like(row)
            res = np.zeros_like(row)
        else:
            corrected_area, iso_fraction, res, m_enr = corrector_HR.correct(row)
        X_corr_abs_ion[idx, :] = np.array(corrected_area)
        X_corr_norm_ion[idx, :] = np.array(iso_fraction)
        X_isocor_res[idx, :] = np.array(res)

    return X_corr_abs_ion, X_corr_norm_ion, X_isocor_res


def apply_isocorr_image(ion: str, X: np.array):
    """
    Use isocor to correct isotopologue distributions of certain ion for isotopologue abundances for an ion image.

    Args:
        ion: formula of the molecule, for which correction is done
        X: 
            matrix of raw isotopologue abundances NxHxW,
            where HxW is the size of the image,
            N should be equal to number of tracer atoms + 1.

    Returns:
        Tuple of np.arrays with the same shape as X:
            X_corr_abs_ion - corrected distributions with the same sum of all peaks as before correction
            X_corr_norm_ion - normalized corrected distributions, sum of all peak intensities == 1
            X_isocor_res - residuals of isocor fitfor each peak
    """
    corrector_HR = isocor.mscorrectors.MetaboliteCorrectorFactory(
        ion,
        tracer="13C",
        resolution=140000,
        mz_of_resolution=200,
        charge=-1,
        correct_NA_tracer=True,
    )
    X_corr_abs_ion = np.zeros(X.shape)
    X_corr_norm_ion = np.zeros(X.shape)
    X_isocor_res = np.zeros(X.shape)
    start = time.time()
    for iy, ix in tqdm(np.ndindex(X.shape[1:]), total = X.shape[1] * X.shape[2], desc="Isocor"):
        dist = X[:, iy, ix]
        if (iy % (X.shape[1] - 1) == 0):
            end = time.time()
            # print(f"{iy} {ix} out of {X.shape}, {end - start} sec")
            start = time.time()
            
        # Not do calculations for AMs only with 0s
        if np.all(dist == 0):
            corrected_area = np.zeros_like(dist)
            iso_fraction = np.zeros_like(dist)
            res = np.zeros_like(dist)
        else:
            corrected_area, iso_fraction, res, m_enr = corrector_HR.correct(dist)
        X_corr_abs_ion[:, iy, ix] = np.array(corrected_area)
        X_corr_norm_ion[:, iy, ix] = np.array(iso_fraction)
        X_isocor_res[:, iy, ix] = np.array(res)

    return X_corr_abs_ion, X_corr_norm_ion, X_isocor_res


def correct_isotopic_abundance(adata_raw: sc.AnnData, int_thresh: float = 200):
    """
    Perform natural isotopic abundance correction for single-cell measurements.

    Args:
        adata_raw: AnnData object with single-cell measurements of isotopologue distributions
        int_thresh: Mark AMs with raw sum intensity per fatty acid less than int_thr.

    Returns:
        AnnData with correction results added as new layers.
        X - raw intensities
        norm - raw intensities, divided by sum of distribution for each ion
        corr_abs - corrected intensities, sum of distribution for each ion == sum before correction
        corr_norm - corrected intensities, sum of distribution for each ion == 1
        corr_res - IsoCor residuals of fit for each peak
    """

    adata = adata_raw.copy()
    adata.var["x_idx"] = np.arange(adata.shape[1])
    X_norm = np.zeros(adata.X.shape)
    X_corr_abs = np.zeros(adata.X.shape)
    X_corr_norm = np.zeros(adata.X.shape)
    X_isocor_res = np.zeros(adata.X.shape)
    X_low_int = np.zeros(adata.X.shape)
    
    for idx, ion in adata.var[adata.var["M+"] == 0].iterrows():
        print(ion.ionFormula, flush=True)
        ion_idx = adata.var.loc[
            adata.var.unlabeled_annotation_id == ion.unlabeled_annotation_id, "x_idx"
        ]
        isotopologue_sum = adata.X[:, ion_idx].sum(axis=1)[:, np.newaxis]
        X_norm[:, ion_idx] = np.divide(
            adata.X[:, ion_idx], isotopologue_sum, where=(isotopologue_sum > 0)
        )

        X_low_int[:, ion_idx] = X_low_int[:, ion_idx] + (isotopologue_sum < int_thresh)

        X_corr_abs_ion, X_corr_norm_ion, X_isocor_res_ion = apply_isocorr(
            ion.ionFormula, adata.X[:, ion_idx]
        )
        X_corr_abs[:, ion_idx] = X_corr_abs_ion
        X_corr_norm[:, ion_idx] = X_corr_norm_ion
        X_isocor_res[:, ion_idx] = X_isocor_res_ion

    adata.layers["norm"] = X_norm
    adata.layers["corr_abs"] = X_corr_abs
    adata.layers["corr_norm"] = X_corr_norm
    adata.layers["corr_res"] = X_isocor_res
    adata.layers["low_int"] = X_low_int

    return adata


def overlapping_am(cell_id, memb_table):
    """ Get list of AMs overlapping with certain cell """
    return list(memb_table.columns[memb_table.loc[cell_id, :] > 0])


def median_norm(adata_sc: sc.AnnData, adata_am: sc.AnnData, overlap_df: pd.DataFrame):
    """
    Single-cell deconvolution method: for each cell assign median of normalized measurements
    of all AMs overlapping with a given cell, then divide isotopologue distribution for each molecule
    by the sum of all peaks to have total area = 1.

    Args:
        adata_sc: 
            AnnData with single-cell measurments as output by SpaceM
            (has fluorescence and cell shape readouts, but appropriate single-cell normalization method is not implemented)
        adata_am:
            AnnData with single-AM measurements normalized with IsoCor
        overlap_df:
            DataFrame output by SpaceM, lists which AMs overlap with which cell

    Returns:
        AnnData with correction results added as new layers.
        X - values output by SpaceM
        corr_norm - median of normalized measurements of all AMs overlapping with a given cell
            divided by sum over all isotopologue peaks for given ion
        min_val - minimum of distribution over all AMs overlapping with a given cell
        max_val - maximum of distribution over all AMs overlapping with a given cell
        variance - variance of distribution over all AMs overlapping with a given cell
    """
    memb_table = pd.pivot_table(
        overlap_df, index="cell_id", columns="am_id", values="area"
    ).fillna(0)
    am_ids = {id: id - 1 for id in memb_table.columns}
    memb_table = memb_table.rename(columns=am_ids)
    memb_table = memb_table[memb_table.columns.sort_values()]

    adata_sc_norm = adata_sc.copy()
    old_X = np.array(adata_am.layers["corr_norm"])
    low_int = np.array(adata_am.layers["low_int"])
    new_X = np.zeros(adata_sc_norm.X.shape)
    min_val = np.zeros(adata_sc_norm.X.shape)
    max_val = np.zeros(adata_sc_norm.X.shape)
    variance = np.zeros(adata_sc_norm.X.shape)

    for cell_idx, cell_id in enumerate(memb_table.index):
        ams = overlapping_am(cell_id, memb_table)
        low_int_filter = low_int[ams, :].astype(bool)
        overlap = old_X[ams, :]

        # filter AMs: remove those with too low sum intensity
        overlap_filtered = overlap
        overlap_filtered[low_int_filter] = np.nan
        
        norm_val = np.nanmedian(overlap_filtered, axis=0)

        min_val_i = np.nanmin(overlap_filtered, axis=0)
        max_val_i = np.nanmax(overlap_filtered, axis=0)
        variance_i = np.nanvar(overlap_filtered, axis=0)
        
        new_X[cell_idx, :] = norm_val
        min_val[cell_idx, :] = min_val_i
        max_val[cell_idx, :] = max_val_i
        variance[cell_idx, :] = variance_i

    adata_sc_norm.layers["corr_norm"] = new_X

    # Normalize per fatty acid again
    for fa in adata_sc_norm.var.unlabeled_annotation_id.unique():
        fa_sum = (
            adata_sc_norm[:, adata_sc_norm.var.unlabeled_annotation_id == fa]
            .layers["corr_norm"]
            .sum(axis=1)
        )
        adata_sc_norm[:, adata_sc_norm.var.unlabeled_annotation_id == fa].layers[
            "corr_norm"
        ] /= fa_sum[:, None]

    adata_sc_norm.layers["min_val"] = min_val
    adata_sc_norm.layers["max_val"] = max_val
    adata_sc_norm.layers["variance"] = variance

    return adata_sc_norm


def load_json(json_path) -> dict:
    with open(json_path, "r") as json_file:
        json_dict = json.load(json_file)
    return json_dict


def get_dataset_name(dataset_path: Path) -> str:
    config_path = Path(dataset_path) / "config.json"
    if config_path.exists():
        config = load_json(dataset_path)
        return config["input"]["dataset_name"]
    else:
        return dataset_path.name


def concatenate_adatas(spacem_dir, sc_analysis_dir):
    ds_dirs = [fname.parent for fname in spacem_dir.glob("**/config.json")]
    configs = [load_json(fname / "config.json") for fname in ds_dirs]
    ds_names = [config["input"]["dataset_name"] for config in configs]
    sm_paths = [
        ds_dir / "analysis" / sc_analysis_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]

    # Read adatas
    adatas = []
    for sm_path in sm_paths:
        if sm_path.is_file():
            adatas.append(sc.read_h5ad(sm_path))

    for adata in adatas:
        print(adata.obs.dataset.unique()[0], "cells: %d, ions: %d" % adata.shape)
        # Save cell id
        adata.obs["cellid"] = adata.obs.index

    adata_concat = sc.AnnData.concatenate(
        *adatas,
        batch_categories=ds_names,
        batch_key="batch",
        join="inner",
        fill_value=0.0,
    )

    print(f"concatenated:\t{len(adata_concat.var.index)} annotations")

    return adata_concat


if __name__ == "__main__":
    # Parse parameters
    parser = argparse.ArgumentParser(
        description="""Script to normalize and concatenate single-cell measurements within one slide.

        The script takes as input several SpaceM-generated directories, 
        performs natural isotope abundance correction for each ion's isotopologue distribution in each ablation mark,
        then summarizes normalized measurements from each ablation mark overlapping with each cell.
        """
    )
    parser.add_argument(
        "spacem_dir",
        type=str,
        help="path to directory, containing several SpaceM directories from the same slide.",
    )
    parser.add_argument(
        "anndata_path",
        type=str,
        help="path to directory, in which to store resulting anndata objects",
    )

    args = parser.parse_args()
    spacem_dir = Path(args.spacem_dir)

    print("Create directory for storing anndata objects if it doesn't exist")
    anndata_path = Path(args.anndata_path)
    anndata_path.mkdir(parents=True, exist_ok=True)

    print("Concatenate ablation mark AnnDatas")
    sc_analysis_dir = "ablation_mark_analysis"
    adata_concat = concatenate_adatas(spacem_dir, sc_analysis_dir)
    output_file = anndata_path / "adata_concat_am.h5ad"
    adata_concat.write(output_file.with_suffix(".h5ad"))

    print("Normalize and correct for natural isotopic abundance")
    adata_corrected = correct_isotopic_abundance(adata_concat)
    output_file = anndata_path / "adata_concat_am_isocor.h5ad"
    adata_corrected.write(output_file.with_suffix(".h5ad"))

    print("Correct and normalize AM AnnData separately for each well")
    sc_analysis_dir = "ablation_mark_analysis"
    sc_analysis_norm_dir = "ablation_mark_analysis_norm"

    ds_dirs = [fname.parent for fname in spacem_dir.glob("**/config.json")]
    configs = [load_json(fname / "config.json") for fname in ds_dirs]
    ds_names = [config["input"]["dataset_name"] for config in configs]
    sm_paths = [
        ds_dir / "analysis" / sc_analysis_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]
    sm_paths_norm_dirs = [
        ds_dir / "analysis" / sc_analysis_norm_dir for ds_dir in ds_dirs
    ]
    sm_paths_output = [
        ds_dir / "analysis" / sc_analysis_norm_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]

    # Read adatas
    adatas = []
    for sm_path in sm_paths:
        if sm_path.is_file():
            adatas.append(sc.read_h5ad(sm_path))

    for idx, adata in enumerate(adatas):
        print(adata.obs.dataset.unique()[0], "cells: %d, ions: %d" % adata.shape)
        # Save cell id
        adata.obs["cellid"] = adata.obs.index

        adata_corrected = correct_isotopic_abundance(adata)
        print(adata_corrected)

        for fa in adata_corrected.var[adata_corrected.var["M+"] == 0].index:
            isotopologues = adata_corrected.var[
                adata_corrected.var.unlabeled_annotation_id == fa
            ]
            isotopologues = isotopologues[isotopologues["M+"] != 0]
            isotopologues = isotopologues.index
            adata_corrected.obs[fa + "_max"] = (
                np.argmax(adata_corrected[:, isotopologues].layers["corr_norm"], axis=1)
                + 1
            )

        sm_paths_norm_dirs[idx].mkdir(parents=True, exist_ok=True)
        adata_corrected.write(sm_paths_output[idx].with_suffix(".h5ad"))

    print("Single-cell analysis")
    am_analysis_norm_dir = "ablation_mark_analysis_norm"
    sc_analysis_dir = "single_cell_analysis"
    sc_analysis_norm_dir = "single_cell_analysis_norm"

    ds_dirs = [fname.parent for fname in spacem_dir.glob("**/config.json")]
    configs = [load_json(fname / "config.json") for fname in ds_dirs]
    ds_names = [config["input"]["dataset_name"] for config in configs]

    am_sm_paths = [
        ds_dir / "analysis" / am_analysis_norm_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]
    sc_sm_paths = [
        ds_dir / "analysis" / sc_analysis_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]
    sc_paths_output = [
        ds_dir / "analysis" / sc_analysis_norm_dir / "spatiomolecular_adata.h5ad"
        for ds_dir in ds_dirs
    ]
    memb_table_paths = [
        ds_dir / "analysis" / "overlap_analysis2" / "overlap.regions.csv"
        for ds_dir in ds_dirs
    ]

    # Read adatas
    adatas_sc = []
    for sm_path in sc_sm_paths:
        if sm_path.is_file():
            adatas_sc.append(sc.read_h5ad(sm_path))

    adatas_am = []
    for am_sm_path in am_sm_paths:
        if am_sm_path.is_file():
            adatas_am.append(sc.read_h5ad(am_sm_path))

    memb_tables = []
    for df_path in memb_table_paths:
        if df_path.is_file():
            memb_tables.append(pd.read_csv(df_path))

    for save_path, adata_sc, adata_am, overlap_df in zip(
        sc_paths_output, adatas_sc, adatas_am, memb_tables
    ):
        print(adata_sc.obs.dataset.unique()[0], "cells: %d, ions: %d" % adata_sc.shape)
        # Save cell id
        # adata.obs["cellid"] = adata.obs.index

        adata_corrected = median_norm(adata_sc, adata_am, overlap_df)

        save_path.parent.mkdir(parents=True, exist_ok=True)
        adata_corrected.write(save_path.with_suffix(".h5ad"))
        
    print("Concatenate single-cell AnnDatas")
    sc_analysis_dir = "single_cell_analysis_norm"
    adata_concat = concatenate_adatas(spacem_dir, sc_analysis_dir)
    output_file = anndata_path / "adata_concat.h5ad"
    adata_concat.write(output_file.with_suffix(".h5ad"))