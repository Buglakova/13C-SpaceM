import json
from pathlib import Path
import json
from skimage import io


def load_json(json_path):
    with open(json_path) as json_file:
        json_dict = json.load(json_file)
    return json_dict


def add_spacem_images(dataset_path: Path, cell_segmentation="cell_segmentation_external"):
    ds_name = load_json(dataset_path / "config.json")["input"]["dataset_name"]

    state_path = dataset_path / "analysis"/ "pipeline_state.json"

    state = load_json(state_path)

    # Cell segmentation
    cell_masks = io.imread(dataset_path / state["results"][cell_segmentation]["cell_mask"])

    # Pre-maldi
    premaldi_paths = load_json(dataset_path / state["results"]["transformation"]['cropped_pre_maldi_channels_paths'])
    cropped_premaldi = {}
    for chan in premaldi_paths.keys():
        cropped_premaldi[chan] = io.imread(dataset_path / premaldi_paths[chan])

    # Post-maldi
    cropped_postmaldi = io.imread(dataset_path / "analysis" / "transformation" / "cropped_post_maldi_channels" / "img_t1_z1_c0.tif")

    # AM segmentation
    am_masks = io.imread(dataset_path / state["results"]["ion_image_registration"]['ablation_mark_mask'])

    # Overlap analysis
    overlap_masks = io.imread(dataset_path / state["results"]["overlap_analysis2"]['overlap_labels'])

    return ds_name, cell_masks, cropped_premaldi, cropped_postmaldi, am_masks, overlap_masks