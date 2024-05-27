from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

from sciso.imzml_utils import peaks_df_to_images, extract_peaks, extract_tic
from sciso.spacem_ion_images import IonImages, write_ion_metadata

from sciso.extract_isotopologue_ion_images import get_iso_ion_list

from sciso.plot import *

from sciso.sc_analysis import apply_isocorr_image

import argparse
from tqdm import tqdm

import gc



def normalize_fa(imzml_path, metadata_fa, name, output_dir):
    fa_dir = output_dir / name
    fa_dir.mkdir(exist_ok=True)
    fa_dir_svg = fa_dir / "svg"
    fa_dir_svg.mkdir(exist_ok=True)
    
    # Extract and plot TIC
    tic_image = extract_tic(imzml_path).T
    sns.heatmap(tic_image, cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True, vmin=0)
    plt.tight_layout()
    plt.savefig(fa_dir / f"tic_image.png", dpi=300, bbox_inches='tight')
    plt.savefig(fa_dir_svg / f"tic_image.svg", bbox_inches='tight')
    plt.close()
    
    # Extract ion images
    tol = 5
    print(tol)
    coords_df, peaks = extract_peaks(imzml_path, metadata_fa, tol_ppm=tol, tol_mode="orbitrap", base_mz=200)
    print("Convert to pandas")
    images = []
    for peak in peaks:
        image = peaks_df_to_images(coords_df, peak["peaks_df"])
        images.append(image[1].T)
    images = np.array(images)

    print("Shape of extracted images", images.shape)
    
    # Save extracted ion images
    print("Save images as an IonImages object")
    ion_images_obj = IonImages(metadata=metadata_fa, shape=images.shape[1:3], array=images)
    ion_images_path = fa_dir / "ion_images.pickle"
    ion_images_metadata_path = fa_dir / "ion_images_metadata.csv"
    
    with open(ion_images_path, "wb") as f:
        pickle.dump(ion_images_obj, f)
    
    write_ion_metadata(ion_images_obj.metadata, ion_images_metadata_path)
    
    print("Save images as numpy")
    with open(fa_dir / "ion_images.npy", "wb") as f:
        np.save(f, images, allow_pickle=False)
    
    with open(fa_dir / "tic_image.npy", "wb") as f:
        np.save(f, tic_image, allow_pickle=False)
    
    # Plot sum of all peaks 
    sum_image = np.sum(images, axis=0)
    sns.heatmap(sum_image, cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True, vmin=0)
    plt.tight_layout()
    plt.savefig(fa_dir / f"sum_image.png", dpi=300, bbox_inches='tight')
    plt.savefig(fa_dir_svg / f"sum_image.svg", bbox_inches='tight')
    plt.close()
    
    # Plot sum of all peaks normalized by TIC
    sns.heatmap(sum_image / tic_image, cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True, vmin=0)
    plt.tight_layout()
    plt.savefig(fa_dir / f"sum_image_tic_norm.png", dpi=300, bbox_inches='tight')
    plt.savefig(fa_dir_svg / f"sum_image_tic_norm.svg", bbox_inches='tight')
    plt.close()
    
    
    # Normalize images
    formula = metadata_fa["formula"].unique()[0]
    X_corr_abs_ion, X_corr_norm_ion, X_isocor_res = apply_isocorr_image(formula, images)
    ion_images_corr_obj = IonImages(metadata=metadata_fa, shape=images.shape[1:3], array=X_corr_norm_ion)
    ion_images_path = fa_dir / "ion_images_norm.pickle"
    
    with open(ion_images_path, "wb") as f:
        pickle.dump(ion_images_corr_obj, f)
    
    with open(fa_dir / "ion_images_norm.npy", "wb") as f:
        np.save(f, X_corr_norm_ion, allow_pickle=False)
    
    # Plot normalized images
    for idx in range(n):
        plt.figure()
        plt.title(f"{name} M+{idx}")
        if idx == 0:
            vmax = 1
        else:
            vmax = 0.3
        sns.heatmap(X_corr_norm_ion[idx], cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True, vmin=0, vmax=vmax)
        plt.tight_layout()
        plt.savefig(fa_dir / f"norm_{idx}.png", dpi=300, bbox_inches='tight')
        plt.savefig(fa_dir_svg / f"norm_{idx}.svg", bbox_inches='tight')
        plt.close()
    
    
    # Plot labeling degree (1 - M+0)
    sns.heatmap(1 - X_corr_norm_ion[0], cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True, vmin=0, vmax=1)
    plt.tight_layout()
    plt.savefig(fa_dir / f"labeling_degree.png", dpi=300, bbox_inches='tight')
    plt.savefig(fa_dir_svg / f"labeling_degree.svg", bbox_inches='tight')
    plt.close()
    
    
    # Fit binomial for every pixel
    p, uptake, mean_p, mean_M, fit_success = fit_binomial_image(X_corr_norm_ion, formula)
    
    # Plot results of fitting
    for feature, feat_name, img in zip(["p", "uptake", "mean_p", "mean_M"], ["p", "uptake", "average p", "average M"], [p, uptake, mean_p, mean_M]):
        plt.figure()
        plt.title(f"{name} {feat_name}")
        sns.heatmap(img, cmap="cividis", cbar=True, linewidths=0, xticklabels=False, yticklabels=False, square=True)
        plt.tight_layout()
        plt.savefig(fa_dir / f"{feature}.png", dpi=300, bbox_inches='tight')
        plt.savefig(fa_dir_svg / f"{feature}.svg", bbox_inches='tight')
        plt.close()
    
        with open(fa_dir / f"{feature}.npy", "wb") as f:
            np.save(f, img, allow_pickle=False)
        
    # return p, uptake
    



if __name__ == "__main__":
    # Parse parameters
    parser = argparse.ArgumentParser(
        description="""Script to run isotope abundance correction and normalization for imzml.

        The script takes as input path to imzml and list of ions to annotate, 
        performs natural isotope abundance correction for each ion's isotopologue distribution in each ablation mark,
        then saves it as plots and spacem ion image objects.
        """
    )
    parser.add_argument(
        "imzml",
        type=str,
        help="imzml file to process",
    )
    parser.add_argument(
        "metadata",
        type=str,
        help="list of ions to analyze in metaspace output format",
    )
    
    parser.add_argument(
        "output_dir",
        type=str,
        help="path to directory, in which to store resulting plots and extracted ion images",
    )

    args = parser.parse_args()
    
    imzml_path = args.imzml

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Open metadata, generate list of isotopologues
    
    ion_list_path = Path(args.metadata)
    ion_list = pd.read_csv(ion_list_path)
    ion_list.head() 
    metadata_new = get_iso_ion_list(ion_list)
    print(metadata_new)
    
    metadata_new.to_csv(output_dir / "isotopologue_list.csv")
    
    # Iterate over fatty acids
    X_corr_norm_ion = []
    p = [] 
    uptake = []
    
    for name in tqdm(metadata_new.name.unique(), desc="Fatty acids", ascii=' >='):
        metadata_fa = metadata_new[metadata_new["name"] == name]
        n = len(metadata_fa)
        print(name, n)
        
        normalize_fa(imzml_path, metadata_fa, name, output_dir)
        n_unreachable = gc.collect()
        print(f"Number of unreachable objects collected by GC: {n_unreachable}")

        
        
