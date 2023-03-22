import pickle
from collections.abc import Iterable
from pathlib import Path
from typing import Optional, Tuple, Union, List
from dataclasses import dataclass

import numpy as np
import pandas as pd
import json


# Class and helper functions copied from old SpaceM
@dataclass
class IonImages:
    shape: Tuple  # shape of 2D images (height, width)
    metadata: pd.DataFrame  # with columns "annotation_id", "formula", "adduct", "msm"
    array: np.ndarray  # shape=(image_num, y, x)

    def get_for_pixel(
        self, x: int, y: int, ion: Union[int, str] = None
    ) -> Union[np.ndarray, float]:
        pixel_spectrum = self.array[:, y, x]
        if ion is None:
            return pixel_spectrum
        elif isinstance(ion, int):
            return pixel_spectrum[ion]
        else:
            ion_index = np.where((self.metadata.get_for_annotation.values == ion))[0][0]
            return float(pixel_spectrum[ion_index])

    def get_for_pixel_index(
        self, pixel_index: int, start: int = 0, ion: Union[int, str] = None
    ) -> Union[np.ndarray, float]:
        pixel = self.array.reshape((-1, np.prod(self.shape[:2])))[:, pixel_index - start]
        if ion is None:
            # All ion intensities for that pixel
            return pixel
        elif isinstance(ion, int):
            # Ion intensity for that pixel and given ion index
            return float(pixel[ion])
        elif isinstance(ion, Iterable) and all(isinstance(i, int) for i in ion):
            # Ion intensity for that pixel and given ion indices
            return pixel[ion]
        elif isinstance(ion, str):
            # Ion intensity for that pixel and given ion name (tuple of formula and adduct)
            ion_index = np.where((self.metadata.annotation_id.values == ion))[0][0]
            return float(pixel[ion_index])

    def get_for_ion(self, ion: Union[int, str]) -> np.ndarray:
        ion_index = np.where((self.metadata.annotation_id.values == ion))[0][0]
        return self.array[ion_index, :, :]


def slice_region_of_interest(
    ion_images: IonImages, region_of_interest: Tuple[Tuple[int, int], Tuple[int, int]]
) -> IonImages:
    (x_min, y_min), (x_max, y_max) = region_of_interest
    return IonImages(
        shape=(y_max - y_min, x_max - x_min),
        metadata=ion_images.metadata,
        array=ion_images.array[y_min:y_max, x_min:x_max],
    )


def write_ion_metadata(ion_metadata: pd.DataFrame, output_path: Path):
    for col in ["databases", "moleculeNames", "moleculeIds"]:
        if col in ion_metadata.columns:
            ion_metadata = _stringify_column_list(ion_metadata, column=col)
    keep_index = ion_metadata.index.name == "annotation_id"
    ion_metadata.to_csv(output_path, index=keep_index)


def read_ion_metadata(
    ion_metadata_path: Optional[Path] = None,
    mass_spec_annotation_dataset_path: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Reads ion image metadata.

    The metadata was originally stored together with the ion images in a pickle file.

    Args:
        ion_metadata_path: A .csv or a .pickle file. If not known, the parent folder must be given
            as mass_spec_annotation_dataset_path.
        mass_spec_annotation_dataset_path: A folder containing either ion_images_metadata.csv or
            ion_images_metadata.pickle

    Returns:
        A dataframe with index annotation_id, containing metadata for each ion image.
    """
    assert ion_metadata_path is not None or mass_spec_annotation_dataset_path is not None
    if ion_metadata_path is None:
        ion_metadata_path = _find_ion_metadata_file(mass_spec_annotation_dataset_path)
    if ion_metadata_path.suffix.lower() == ".pickle":
        with open(ion_metadata_path, "rb") as f:
            ion_images = pickle.load(f)
            ion_metadata = ion_images.metadata
        if "annotation_id" not in ion_metadata.columns:
            ion_metadata["annotation_id"] = _pd_create_unique_column(
                ion_metadata, columns=["formula", "adduct"], format_string="{formula}{adduct}"
            )
    else:
        ion_metadata = pd.read_csv(ion_metadata_path)
        # The moleculeNames column originally contains a list of strings, which would be stringified
        # when saved as csv. To preserve it, we stringify and later parse it as JSON.
        for col in ["databases", "moleculeNames", "moleculeIds"]:
            if col in ion_metadata.columns:
                ion_metadata = _parse_stringified_columns(ion_metadata, column=col)
    if ion_metadata.index.name != "annotation_id":
        ion_metadata.set_index("annotation_id", drop=True, inplace=True)
    return ion_metadata


def _find_ion_metadata_file(mass_spec_annotation_dataset_path: Path) -> Path:
    ion_metadata_csv_path = mass_spec_annotation_dataset_path / "ion_images_metadata.csv"
    ion_metadata_pickle_path = mass_spec_annotation_dataset_path / "ion_images_metadata.pickle"
    if ion_metadata_csv_path.exists():
        ion_metadata_path = ion_metadata_csv_path
    elif ion_metadata_pickle_path.exists():
        ion_metadata_path = ion_metadata_pickle_path
    else:
        raise RuntimeError(
            "Either ion_images_metadata.csv or ion_images_metadata.pickle must exist"
        )
    return ion_metadata_path


def _pd_create_unique_column(
    df: pd.DataFrame, columns: Union[str, List[str]], format_string: str
) -> pd.Series:
    """
    Joins values from string columns into a single, unique column.

    When there are multiple columns intended as primary keys, a concatenated dataframe may contain
    rows with duplicate values for these keys. This function allows to create a new column without
    duplicates by joining the column values and enumerating duplicates. If there is no duplicate,
    the joined string value is used, otherwise from the second occurrence onwards, an underscore
    and a number strating from one is appended.

    Args:
        df: A dataframe
        columns: The names of columns to use as primary keys.
        format_string: A string used as template containing the column names in curly braces.

    Returns:
        A series of unique strings
    """
    if isinstance(columns, str):
        columns = [columns]
    if df.empty:
        return pd.Series([], dtype=str)
    # Create groups of identical column values.
    gb = df.groupby(columns)
    # For each row, get the size of its group.
    groups_lengths = _pd_groupby_count_element_wise(gb)
    # For each row, get its cardinal number within its group.
    groups_enumerated = gb.cumcount()
    # Join the columns.
    formatted = df.apply(lambda s: format_string.format(**{k: s[k] for k in columns}), axis=1)
    # Append the number if needed, i.e. group larger than one, skipping suffix for first element.
    unique = formatted.where(
        (groups_lengths.values == 1) | (groups_enumerated == 0),
        formatted + "_" + groups_enumerated.astype(str),
    )
    assert unique.is_unique
    return unique


def _stringify_column_list(df: pd.DataFrame, column: str = "moleculeNames") -> pd.DataFrame:
    """
    If the data frame column is a list, stringify it into JSON.

    Args:
        df: A dataframe with a column that contains lists (e.g. of molecule names)
        column: The name of the column containing lists

    Returns:
        A dataframe where the column of molecule names contains a JSON string per row
    """
    if not df.empty and pd.api.types.is_list_like(df[column].iloc[0]):
        df = df.copy()
        df[column] = df[column].apply(lambda l: json.dumps(l))
    return df


def _parse_stringified_columns(df: pd.DataFrame, column: str = "moleculeNames") -> pd.DataFrame:
    """
    If the column of df is stringified, parses it into a list of strings.

    Args:
        df: A dataframe with a column of stringified lists
        column: The name of the column containing stringified lists

    Returns:
        A dataframe where the column of molecule names is a list of strings
    """
    if not df.empty and not pd.api.types.is_list_like(df[column].iloc[0]):
        df = df.copy()
        df[column] = df[column].astype(str).apply(_parse_json_or_py_list)
    return df


def _pd_groupby_count_element_wise(gb: pd.core.groupby.GroupBy) -> pd.Series:
    # Like pandas GroupBy.count() but returns original number of rows, not aggregation.
    return gb.transform(lambda series: len(series)).iloc[:, 0]


def _parse_json_or_py_list(s: str) -> List[str]:
    """
    Tries to turn the stringified list expression s back into
    a string. First tries to deserialize the json string, if
    that fails, tries the legacy method using eval

    Args:
        s: the string that should contain a stringified list

    Returns:
        A list of strings
    """
    try:
        return json.loads(s)
    except ValueError:
        # Fallback: Earlier code did not handle moleculeNames specifically and it became
        # stringified as a Python list, not a JSON string.
        return eval(s)  

