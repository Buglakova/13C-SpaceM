from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import shutil
import pickle
from collections import Counter

from sciso.imzml_utils import peaks_df_to_images, extract_peaks
from sciso.formula_parser import (
    parse_formula,
    _format_formula,
    format_ion_formula,
    generate_ion_formula,
)
from sciso.spacem_ion_images import IonImages, write_ion_metadata


def generate_isotopologue(elements: Counter, mz: float, n_C13: int):
    """
    Generate formula of isotopologue with given number of 13C.

    Args:
        elements: Counter object containing elements and number of atoms.
        mz: Mass of the molecule.
        n_C13: Number of C atoms to replace.

    Returns:
        Tuple:
            Counter object with 12C replaced with 13C
            Mass of resulting isotopologue
    """

    C12_MASS = 12.0
    C13_MASS = 13.003355
    dm = C13_MASS - C12_MASS
    C13_elements = elements.copy()
    C13_elements["C"] -= n_C13
    C13_elements["[13C]"] += n_C13
    # print(n_C13, mz, mz + dm * n_C13)
    return C13_elements, mz + dm * n_C13


def create_iso_ion_list(ion_list_row: pd.Series):
    """
    For a row of molecule database generate table with all isotopologues.

    Args:
        ion_list_row: 
            Row of molecule database (pd.Series), obligatory column:
            "formula" - formula of the ion

    Returns:
        DataFrame with all isotopologues for a given ion
    """
    assert "formula" in ion_list_row.index, 'ion_list_row must have an "formula" column'

    if "annotation_id" in ion_list_row.index:
        ion_list_row["unlabeled_annotation_id"] = ion_list_row["annotation_id"]

    ion_list_row["M+"] = 0
    ion_list_row["dm"] = 0
    df = [ion_list_row]
    elements = Counter(dict(parse_formula(ion_list_row["formula"])))
    for i in range(1, elements["C"] + 1):
        new_row = ion_list_row.copy()
        iso_elements, iso_mz = generate_isotopologue(elements, ion_list_row["mz"], i)
        new_row["formula"] = _format_formula(iso_elements)
        new_row["mz"] = iso_mz
        new_row["annotation_id"] = format_ion_formula(
            _format_formula(iso_elements), "-H"
        )
        new_row["ion"] = format_ion_formula(
            _format_formula(iso_elements), "-H", charge=-1
        )
        new_row["ionFormula"] = generate_ion_formula(
            _format_formula(iso_elements), "-H"
        )
        new_row["M+"] = i
        new_row["dm"] = iso_mz - ion_list_row["mz"]
        df.append(new_row)
    df = pd.DataFrame(df)
    return df


def get_iso_ion_list(ion_list: pd.DataFrame):
    """
    Iterate over rows of ion list and generate new ion list with
    all possible isotopologues for each ion.

    Args:
        ion_list: DataFrame with obligatory column "formula".

    Returns:
        DataFrame in the same format but with isotopologues.
    """
    iso_dfs = []
    for idx, row in ion_list.iterrows():
        iso_dfs.append(create_iso_ion_list(row))

    iso_df = pd.concat(iso_dfs).reset_index(drop=True)
    return iso_df


if __name__ == "__main__":
    # Parse parameters
    parser = argparse.ArgumentParser(
        description="""Script to create list of isotopologues and extract ion images for them.

        The script takes as input a list of formulas + adducts in the metadata format from SpaceM, imzml file and SpaceM directory where to write the result in SpaceM-compatible format.
        """
    )
    parser.add_argument("ion_list", type=str, help="Path to the ion list")
    parser.add_argument(
        "spacem_dir",
        type=str,
        help="path to SpaceM directory, where to save resulting metadata and ion images",
    )

    args = parser.parse_args()

    ion_list_path = args.ion_list
    spacem_dir = Path(args.spacem_dir)
    imzml_path = spacem_dir / "MassSpectrometry"

    imzml_path = list(imzml_path.glob("*.imzML"))[0]
    spacem_dir = spacem_dir / "analysis" / "metaspace"

    # Create new metadata with isotopologues
    ion_list = pd.read_csv(ion_list_path)
    print(ion_list)

    # Calculate masses for all isotopologues
    metadata_new = get_iso_ion_list(ion_list)

    # Extract ion images
    coords_df, peaks = extract_peaks(
        imzml_path, metadata_new, tol_ppm=3, tol_mode="orbitrap", base_mz=200
    )

    images = []
    for peak in peaks:
        image = peaks_df_to_images(coords_df, peak["peaks_df"])
        images.append(image[1].T)
    images = np.array(images)

    print(images.shape)

    # Write result
    ion_images_obj = IonImages(
        metadata=metadata_new, shape=images.shape[1:3], array=images
    )

    ion_images_path = spacem_dir / "ion_images.pickle"

    # If dataset was already annotated, rename previous annotations
    if ion_images_path.is_file():
        shutil.move(ion_images_path, spacem_dir / "ion_images_old.pickle")

    with open(ion_images_path, "wb") as f:
        pickle.dump(ion_images_obj, f)

    ion_images_metadata_path = spacem_dir / "ion_images_metadata.csv"

    # If dataset was already annotated, rename previous annotations
    if ion_images_metadata_path.is_file():
        shutil.move(
            ion_images_metadata_path, spacem_dir / "ion_images_metadata_old.csv"
        )
    write_ion_metadata(ion_images_obj.metadata, ion_images_metadata_path)
