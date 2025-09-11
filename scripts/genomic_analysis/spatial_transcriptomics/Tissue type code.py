import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.collections import PatchCollection
import imageio.v3 as iio


base_dir = # enter base directory
output_dir = # enter output directory
image_dir = # enter image directory
os.makedirs(output_dir, exist_ok=True)

he_resolutions = {} # enter resolution of HE images

cytassist_resolution = # enter resolution of cytassist image
REFINED_GENE_PANELS = {
    "Squamous": ["DSG3", "KRT5", "KRT14", "TP63"],
    "Stroma": ["ACTA2", "PDGFRA", "COL1A1", "COL3A1", "FAP", "MMP2", "MMP9"],
    "Barretts": ["MUC2", "TFF3", "REG4", "CDX2"],
    "Tumor": ["MKI67", "SPINK1", "ERBB2", "CLDN4"]
}

COLOR_MAP = {
    'Tumor': 'royalblue',
    'Barretts': 'skyblue',
    'Squamous': 'green',
    'Stroma': 'palegreen',
}
ALPHA_SETTING = # choose alpha setting
sample_ids = [] # choose samples

for sample_id in sample_ids:
    print(f"\n▶ Processing sample: {sample_id}")

    try:
        if sample_id not in he_resolutions:
            continue

        h5ad_file_path = os.path.join(base_dir, ) # enter path
        he_image_path = os.path.join(image_dir, ) # enter path

        adata = sc.read_h5ad(h5ad_file_path)
        he_img = iio.imread(he_image_path)
        img_height, img_width, _ = he_img.shape
        DPI = 100 

        fig_width_in = img_width / DPI
        fig_height_in = img_height / DPI


        base_scale = cytassist_resolution / he_resolutions[sample_id]
        params = special_plot_params.get(sample_id, {})
        
        x_scale_factor = params.get("x_scale_factor", 1.0)
        y_scale_factor = params.get("y_scale_factor", 1.0)
        x_shift = params.get("x_shift", 0.0)
        y_shift = params.get("y_shift", 0.0)
        
        # Final scales for each axis
        scale_to_he_x = base_scale * x_scale_factor
        scale_to_he_y = base_scale * y_scale_factor


        panel_counts = pd.DataFrame(index=adata.obs.index)
        for panel_name, gene_list in REFINED_GENE_PANELS.items():
            genes_present = [g for g in gene_list if g in adata.var_names]
            if genes_present:
                panel_counts[panel_name] = np.array(adata[:, genes_present].X.sum(axis=1)).flatten()
            else:
                panel_counts[panel_name] = 0

        counts_st = panel_counts["Stroma"]
        counts_sq = panel_counts["Squamous"]
        counts_ba = panel_counts["Barretts"]
        counts_tu = panel_counts["Tumor"]
        
        adata.obs.loc[counts_st > 0, 'final_cell_type'] = "Stroma"
        adata.obs.loc[counts_sq > 0, 'final_cell_type'] = "Squamous"
        
        both_sq_st_mask = (counts_sq > 0) & (counts_st > 0)
        adata.obs.loc[both_sq_st_mask, 'final_cell_type'] = np.where(
            counts_sq[both_sq_st_mask] > counts_st[both_sq_st_mask], "Squamous", "Stroma"
        )
        
        adata.obs.loc[counts_ba > 0, 'final_cell_type'] = "Barretts"
        adata.obs.loc[counts_tu > 0, 'final_cell_type'] = "Tumor"

        print("  - Generating plot...")
        fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in))
        
        ax.imshow(he_img)
    
        bin_size_px_x = 2 * scale_to_he_x
        bin_size_px_y = 2 * scale_to_he_y

        for cell_type, color in COLOR_MAP.items():
            mask = adata.obs['final_cell_type'] == cell_type
            if not mask.any():
                continue

            coords = adata[mask].obsm['spatial']
            
            scale_vector = np.array([scale_to_he_x, scale_to_he_y])
            shift_vector = np.array([x_shift, y_shift])
            coords_scaled = coords * scale_vector + shift_vector
            
            patches = [Rectangle((x - bin_size_px_x / 2, y - bin_size_px_y / 2), bin_size_px_x, bin_size_px_y) for x, y in coords_scaled]
            collection = PatchCollection(patches, facecolor=color, alpha=ALPHA_SETTING, edgecolor='none')
            ax.add_collection(collection)

        legend_elements = [Patch(facecolor=color, label=label, alpha=ALPHA_SETTING) for label, color in COLOR_MAP.items()]
        ax.legend(handles=legend_elements, loc='best', frameon=True, fontsize='large')
        title = f"" # add title
        ax.set_title(title, fontsize=16)
        ax.axis("off")
        plt.tight_layout(pad=0)

        output_path = os.path.join(output_dir, f"{") # add title
        
        plt.savefig(output_path, dpi=DPI, bbox_inches="tight", pad_inches=0)
        plt.show()
        plt.close(fig)
