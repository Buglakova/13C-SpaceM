import matplotlib.pyplot as plt
from pathlib import Path
import seaborn as sns
import scanpy as sc
import numpy as np
import pandas as pd

from sklearn.metrics import ConfusionMatrixDisplay

from scipy.optimize import minimize
from scipy.stats import binom

from sklearn.linear_model import LogisticRegressionCV
from sklearn.model_selection import train_test_split

from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix, plot_confusion_matrix, roc_curve, roc_auc_score


x_width = 3.3
y_width = 2.5

sns.set_palette("husl")

def set_nature_style():
    mplstyle_path = Path(__file__).parent.resolve() / "nature.mplstyle"
    plt.style.use(mplstyle_path)
    
    
def print_statistics(adata: sc.AnnData, condition_col="condition", ds_col="batch"):
    datasets = adata.obs[ds_col].unique()
    conditions = adata.obs[condition_col].unique()
    
    print("Number of cells per dataset")
    for dataset in datasets:
        print(dataset, len(adata[adata.obs[ds_col] == dataset, :]), " cells")
    
    print("Number of cells per condition")
    for condition in conditions:
        print(condition, len(adata[adata.obs[condition_col] == condition, :]), " cells")
        
    print("Pivot table with condition and labeling state")
    print(adata.obs[["growthConditions", condition_col, ds_col]].pivot_table(columns=["growthConditions", condition_col], aggfunc="count"))
        
    print("Number of ions: ", len(adata.var.ionFormula.unique()))
    print("Number of unique formulas: ", len(adata.var.formula.unique()))
    
    
def get_sum_intensities(adata: sc.AnnData):
    adata_sum = adata.copy()
    for ion in adata.var.unlabeled_annotation_id.unique():
        adata_ion = adata[:, adata.var.unlabeled_annotation_id == ion]
        adata_sum.obs[f"{ion}_sum"] = adata_ion.X.sum(axis=1)
    return adata_sum

def plot_sum_intensities(adata: sc.AnnData, plots_path: Path):
    """
    Plot histograms of total ion intensities per well per ion
    (sum over all isotopologues for each ion)
    for intracellullar AMs.

    Args:
        adata: AnnData for ablation marks
        plots_path:
    """
    adata_am_sum = get_sum_intensities(adata)
    col_names = [f"{ion}_sum" for ion in adata_am_sum.var.unlabeled_annotation_id.unique()]

    for name in col_names:
        g = sns.FacetGrid(adata_am_sum.obs, col="well",  row="slide")
        g.map_dataframe(sns.histplot, x=name)
        g.fig.suptitle(name)
        plt.savefig(plots_path / f"hist_per_slide_{name}.png")
        plt.show()
    
    
def get_iso_dist(adata: sc.AnnData, ion: str, cell_idx: str, layer="X"):
    adata_ion_cell = adata[cell_idx, adata.var.unlabeled_annotation_id == ion]
    if layer == "X":
        spectrum = adata_ion_cell.X[0]
    else:
        spectrum = adata_ion_cell.layers[layer][0]
        
    return spectrum

    
def plot_iso_distribution(adata: sc.AnnData, ion: str, cell_idx: str, layer="X", ax=None, alpha=1, dx=0, color="gray", label=None):
    adata_ion_cell = adata[cell_idx, adata.var.unlabeled_annotation_id == ion]
    spectrum = get_iso_dist(adata, ion, cell_idx, layer=layer)

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(x_width * 2, y_width))
        
    x = np.arange(len(spectrum)) + dx
    if label is None:
        label = cell_idx + " " + ion
    sns.scatterplot(x=x, y=spectrum, s=9, ax=ax, label=label, alpha=alpha, color=color,)
    ax.vlines(x, np.zeros_like(spectrum), spectrum, linewidth=3, alpha=alpha, color=color)
    
    # Set xticks
    xtick_pos = np.arange(0, len(spectrum))
    ax.set_xticks([*xtick_pos, *(xtick_pos + 0.01)])
    xtick_labels_m = ["M+%d"%d for d in range(len(spectrum))]
    xtick_labels_mz = ["\n\n%.3f"%mz for mz in  adata_ion_cell.var.mz]
    xtick_labels = [*xtick_labels_m, *xtick_labels_mz]
    ax.set_xticklabels(xtick_labels, fontsize=6)
    for label in ax.get_xmajorticklabels():
        if 'M+' not in label.get_text(): 
            label.set_fontsize(4)
            label.set_color("silver")
        
    # Set axis labels
    ax.set_xlabel("Isotopologue number")
    ax.set_ylabel("I")
    
    return ax
    
        
def ion_stripplot_bulk(adata, ion_unlabeled_annotation_id, ax=None, color="#18974C"):
    """
    For a given FA plot all isotopologue intensities for all cells
    and add mean lines for each distribution to compare with bulk 

    Args:
        adata: AnnData for cells
        ion_unlabeled_annotation_id: unlabeled ion formula
        ax: _description_. Defaults to None.
    """        
    ion_idx = adata.var.loc[adata.var.unlabeled_annotation_id == ion_unlabeled_annotation_id, "index"]
    # print(ion_idx)
    M = adata.var.loc[adata.var.unlabeled_annotation_id == ion_unlabeled_annotation_id, "M+"]
    mz_mass = adata.var.loc[adata.var.unlabeled_annotation_id == ion_unlabeled_annotation_id, "mz"]
    # print(mz)

    # df = pd.DataFrame(adata.layers["corr_norm"][:, ion_idx], columns=mz)
    df = pd.DataFrame(adata[:, adata.var.unlabeled_annotation_id == ion_unlabeled_annotation_id].layers["corr_norm"], columns=M)
    #     print("unmelted")
    #     print(df.iloc[0: 10])
    df = pd.melt(df, value_vars=df.columns)
    #     print("Melted")
    #     print(df[0: 10])
    # df.mz = [round(x, 2) for x in df.mz.astype(float)]
    # df.mz = mz

    ax = sns.stripplot(data=df, x="M+", y="value", s=0.6, jitter=0.4, zorder=-1, ax=ax, color=color)
    sns.boxplot(data=df, x="M+", y="value", showmeans=True,
            meanline=True,
            meanprops={'color': 'k', 'ls': '-', 'lw': 2},
            medianprops={'visible': False},
            whiskerprops={'visible': False},
            zorder=10,
            showfliers=False,
            showbox=False,
            showcaps=False,
            ax=ax)
    
    ax.set_ylim((0, 1))
    # ax.tick_params(axis="x", labelsize=10);

    ax.set_title(ion_unlabeled_annotation_id)
    
    # Set xticks
    show_mass = False
    if show_mass:
        xtick_pos = M
        ax.set_xticks([*xtick_pos, *(xtick_pos + 0.01)])
        xtick_labels_m = ["M+%d"%d for d in M]
        xtick_labels_mz = ["\n\n%.3f"%mz for mz in mz_mass]
        xtick_labels = [*xtick_labels_m, *xtick_labels_mz]
        ax.set_xticklabels(xtick_labels, fontsize=6)
        for label in ax.get_xmajorticklabels():
            if 'M+' not in label.get_text(): 
                label.set_fontsize(4)
                label.set_color("silver")
    else:
        xtick_pos = M
        ax.set_xticks(xtick_pos)
        xtick_labels_m = ["M+%d"%d for d in M]
        ax.set_xticklabels(xtick_labels_m, fontsize=6)
        ax.tick_params(axis="x", labelrotation=45)
    
    # Set axis labels
    ax.set_xlabel("Isotopologue number")
    ax.set_ylabel("Normalized intensity")
    
    
def model_C18(param, n):
    x = np.arange(0, n)
    uptake = param[0]
    uptake_C16 = param[1]
    p = param[2]
    density = binom.pmf(x, n-1, p)
    dist = np.zeros(n)
    dist[0] += uptake
    dist[1] += (1 - uptake) * uptake_C16 
    dist += (1 - uptake) * (1 - uptake_C16) * density
    return dist


def model_C16(param, n):
    x = np.arange(0, n)
    uptake = param[0]
    p = param[1]
    density = binom.pmf(x, n-1, p)
    dist = np.zeros(n)
    dist[0] += uptake
    dist += (1 - uptake) * density
    return dist


def square_dist(param, iso_even, model):
    n = iso_even.shape[0]
    dist = model(param, n)
    return np.sum(np.square(iso_even - dist))


def fit_binomial(iso_dist, model):
    iso_even = iso_dist[0::2]
    n = iso_even.shape[0]
    x = np.arange(0, n)
    if model == "C16":
        param_min = minimize(square_dist, [iso_even[0], iso_even[1:].argmax() / n], args=(iso_even, model_C16), bounds=((0, 1), (0, 1)))
        # print(param_min)
        return x, n, param_min.x, param_min.success
    else:
        param_min = minimize(square_dist, [iso_even[0], iso_even[1] / iso_even[0], iso_even[2:].argmax() / n], args=(iso_even, model_C18), bounds=((0, 1), (0, 1), (0, 1)))
        return x, n, param_min.x, param_min.success
    
def calculate_mean(iso_dist_even):
    n = len(iso_dist_even) - 1
    x = np.arange(1, len(iso_dist_even) + 1)
    iso_dist_even = (iso_dist_even / np.sum(iso_dist_even))
    mean_M = np.sum(x * iso_dist_even)
    p = mean_M / n
    return mean_M * 2, p
    
    
def fit_binomial_adata(adata):
    layer = "corr_norm"
    for ion in adata.var.unlabeled_annotation_id.unique():
        print("Ion:", ion)
        for cell_idx in adata.obs.index:
            iso_dist = np.array(get_iso_dist(adata, ion, cell_idx, layer))
            if ion.startswith("C18"):
                x, n, param, success = fit_binomial(iso_dist, "C18")
                if not success or np.allclose(iso_dist[3:], np.zeros_like(iso_dist[3:])):
                    param[0] = iso_dist[0]
                    param[1] = iso_dist[2] / iso_dist[0]
                    param[2] = None
                adata.obs.loc[cell_idx, f"{ion}_p"] = param[2]
                adata.obs.loc[cell_idx, f"{ion}_uptake"] = param[0]
                adata.obs.loc[cell_idx, f"{ion}_uptake_palmitate"] = param[1]
                adata.obs.loc[cell_idx, f"{ion}_success"] = int(success)
                # Estimate mean and p without fitting
                iso_dist_even = iso_dist[2::2]
                iso_dist_even[0] = 0
                mean_M, mean_p = calculate_mean(iso_dist_even)
                adata.obs.loc[cell_idx, f"{ion}_mean_p"] = mean_p
                adata.obs.loc[cell_idx, f"{ion}_mean_M"] = mean_M
            else:
                x, n, param, success = fit_binomial(iso_dist, "C16")
                iso_dist_even = iso_dist[2::2]
                if not success or np.allclose(iso_dist[1:], np.zeros_like(iso_dist[1:])):
                    param[0] = iso_dist[0]
                    param[1] = None
                adata.obs.loc[cell_idx, f"{ion}_p"] = param[1]
                adata.obs.loc[cell_idx, f"{ion}_uptake"] = param[0]
                adata.obs.loc[cell_idx, f"{ion}_success"] = int(success)
                # Estimate mean and p without fitting
                iso_dist_even = iso_dist[2::2]
                mean_M, mean_p = calculate_mean(iso_dist_even)
                adata.obs.loc[cell_idx, f"{ion}_mean_p"] = mean_p
                adata.obs.loc[cell_idx, f"{ion}_mean_M"] = mean_M
            
            
def train_classifier(adata, name_prefix, pos_label, plots_path, cond_col="condition"):
    X_train_idx, X_test_idx, y_train_idx, y_test_idx = train_test_split(range(len(adata)), range(len(adata)), random_state=42, test_size=0.3,  shuffle=True)
    X_train = adata.layers["corr_norm"][X_train_idx, :]
    y_train = adata.obs[cond_col].iloc[X_train_idx]
    X_test = adata.layers["corr_norm"][X_test_idx, :]
    y_test = adata.obs[cond_col].iloc[X_test_idx]
    # print(y_test)

    clf = LogisticRegressionCV(random_state=0, max_iter=100000, verbose=1).fit(X_train, y_train)
    
    y_pred = clf.predict(X_test)
    
    print("Accuracy: ", accuracy_score(y_test, y_pred))
    
    
    plot_confusion_matrix(clf, X_test, y_test, cmap="magma")
    plt.savefig(plots_path / f"{name_prefix}_{cond_col}_confusion.png")
    plt.savefig(plots_path / f"{name_prefix}_{cond_col}_confusion.svg")
    plt.show()
    
    fpr, tpr, thresholds = roc_curve(y_test, clf.predict_proba(X_test)[:, 0], pos_label=pos_label)
    print("ROC AUC score: ", roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1]))
    plt.plot(fpr, tpr)
    plt.xlabel("False positive")
    plt.ylabel("True positive")
    plt.plot((0, 1), (0, 1))
    plt.title("ROC AUC score = %.2f"%roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1]))
    plt.savefig(plots_path / f"{name_prefix}_{cond_col}_roc.png")
    plt.savefig(plots_path / f"{name_prefix}_{cond_col}_roc.svg")
    
    return clf