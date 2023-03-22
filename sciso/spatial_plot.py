from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation
from skimage.util import map_array

import numpy as np
import matplotlib.pyplot as plt


def segm_masks_to_contours(masks: np.array):
    contours = find_boundaries(masks, mode="inner")
    contours = contours * masks
    result = np.zeros_like(contours)
    for val in np.unique(contours)[1:]:
        thick_contours = binary_dilation(binary_dilation(contours == val))
        result[thick_contours > 0] =  val
    return result


def map_values(masks, labels, values):
    mapped = map_array(masks, labels, values)
    return mapped


def get_feature_cmap(feature, cmap, rescale=False, lower_percentile=1, upper_percentile=98):
    # Rescale feature between 0 & 1 to make a colormap
    if rescale:
        lower_contrast_limit = np.percentile(feature, lower_percentile)
        upper_contrast_limit = np.percentile(feature, upper_percentile)
        feature_scaled = (feature - lower_contrast_limit) / (upper_contrast_limit - lower_contrast_limit)
    else:
        feature_scaled = feature
    # Cap the measurement between 0 & 1
    feature_scaled[feature_scaled < 0] = 0
    feature_scaled[feature_scaled > 1] = 1

    colors = plt.cm.get_cmap(cmap)(feature_scaled)
    return feature_scaled, colors


def map_obs_feature_to_masks(masks, adata, feature):
    feature_scaled, colors = get_feature_cmap(adata.obs[feature], "Greens")
    mapped_array = map_values(masks, adata.obs.cellid.to_numpy(), feature_scaled.to_numpy())
    return mapped_array


def map_ion_feature_to_masks(masks, adata, feature, layer="corr_norm", int_threshold=None, rescale=False):
    # feature_scaled, colors = get_feature_cmap(adata[:, feature].X.T[0], "Greens")
    feature_scaled, colors = get_feature_cmap(np.array(adata[:, feature].layers["corr_norm"].T[0]), "Greens", rescale=rescale)
    print("count nonzero", np.count_nonzero(feature_scaled))
    if int_threshold is not None:
        feature_scaled[adata[:, feature].X.T[0] < int_threshold] = 0
    print("count nonzero", np.count_nonzero(feature_scaled))
    mapped_array = map_values(masks, adata.obs.cellid.to_numpy(), feature_scaled)
    return mapped_array


def get_ion_image(adata_am, ion, slide, well, image_shape = (100, 100)):
    adata_filtered = adata_am[(adata_am.obs.slide == slide) & (adata_am.obs.well == well)]
    ion_X = adata_filtered[:, ion].layers['corr_norm']
    return ion_X.reshape(image_shape)