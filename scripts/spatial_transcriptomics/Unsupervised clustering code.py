import os
import scanpy as sc
import imageio.v3 as iio
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import seaborn as sns


BASE_DIR     = # sample directory
IMAGE_DIR    = # image directory
CYTASSIST_RES = # µm per pixel in spatial data
HE_RESOLUTIONS = {} # enter resolution of HE images 


MIN_GENES         = 100
MIN_CELLS         = 3
N_HVG             = 2000
N_NEIGHBORS       = 20
N_PCS             = 30
LEIDEN_RESOLUTION = 0.2


sample_id = ""  # enter sample name 

adata_path = os.path.join(BASE_DIR, f"")
he_path    = os.path.join(IMAGE_DIR, f"")


adata = sc.read_h5ad(adata_path)
he_img = iio.imread(he_path)
he_res = HE_RESOLUTIONS[sample_id]
scale_to_he = CYTASSIST_RES / he_res

print("Filtering cells & genes…")
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
print(f" → {adata.n_obs} spots remain after filtering")

print("Normalizing & log-transform…")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()  # save for later


print("Finding highly variable genes…")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_HVG,
    flavor="seurat"
)
adata = adata[:, adata.var.highly_variable].copy()

print("Scaling, PCA, and neighbor graph…")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS)

print("Running Leiden clustering…")
sc.tl.leiden(
    adata,
    resolution=LEIDEN_RESOLUTION
)
adata.obs["cluster"] = adata.obs["leiden"].astype(str)

n_cl = adata.obs["cluster"].nunique()
print(f"Identified {n_cl} clusters.")

# === spatial overlay ===
cytassist_resolution = # enter resolution of cytassist image
he_res = HE_RESOLUTIONS[sample_id]
scale_to_he = cytassist_resolution / he_res

h, w, _ = he_img.shape
dpi = 300
fig_w, fig_h = w / dpi, h / dpi
bin_size = 2 * scale_to_he   # 2 µm bins in pixels

clusters = adata.obs['cluster'].cat.categories.tolist()
palette = sns.color_palette('tab10', len(clusters))
color_map = {cl: palette[i] for i, cl in enumerate(clusters)}

def plot_overlay(cl_list, outpath):
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.imshow(he_img)
    coords_scaled = adata.obsm['spatial'] * scale_to_he

    for cl in cl_list:
        mask = adata.obs['cluster'] == cl
        pts = coords_scaled[mask.values]
        patches = [
            Rectangle((x - bin_size/2, y - bin_size/2), bin_size, bin_size)
            for x, y in pts
        ]
        coll = PatchCollection(patches,
                               facecolor=color_map[cl],
                               edgecolor='none',
                               alpha=0.6)
        ax.add_collection(coll)

    ax.axis('off')
    plt.tight_layout()
    plt.savefig(outpath, dpi=dpi)
    plt.show()
    plt.close()

# combined
plot_overlay(clusters, os.path.join(sample_dir, f"")) # choose title

# individual
for cl in clusters:
    plot_overlay([cl], os.path.join(sample_dir, f"")) # choose title

print("✔ Spatial overlays saved.")

