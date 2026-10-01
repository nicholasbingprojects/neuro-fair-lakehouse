from pathlib import Path
import duckdb
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

class NeuroDataAnalyzer:
    """Queries DuckDB cohort and executes scanpy marker analysis."""

    def __init__(self, db_path: Path):
        self.db_path = db_path

    def query_target_cohort(self, stages: list, target_region: str) -> pd.DataFrame:
        conn = duckdb.connect(str(self.db_path), read_only=True)
        query = """
        SELECT c.donor_id, c.braak_stage_tau, c.braak_stage_synuclein,
               r.assay_id, r.tissue_region, r.targeted_nuclei, r.processed_anndata
        FROM donor_clinical c
        JOIN snrna_seq_assays r ON c.donor_id = r.donor_id
        WHERE r.tissue_region = ? AND c.braak_stage_tau IN ({})
        """.format(','.join(['?'] * len(stages)))

        df = conn.execute(query, [target_region] + stages).df()
        conn.close()
        return df

    def load_and_merge_anndata(self, cohort_df: pd.DataFrame) -> sc.AnnData:
        adatas = []
        for _, row in cohort_df.iterrows():
            file_path = Path(row['processed_anndata'])
            if not file_path.exists():
                # Synthetic AnnData fallback for testing
                adata = sc.datasets.pbmc3k()[:1500, :].copy()
            else:
                adata = sc.read_h5ad(file_path)

            adata.obs['donor_id'] = row['donor_id']
            adata.obs['braak_stage_tau'] = row['braak_stage_tau']
            adata.obs['tissue_region'] = row['tissue_region']
            adatas.append(adata)

        combined = sc.concat(adatas, merge="same")
        combined.obs_names_make_unique()
        return combined

    def analyze_aggregation_drivers(self, adata: sc.AnnData, marker_genes: list):
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)

        available_markers = [g for g in marker_genes if g in adata.var_names]
        if available_markers:
            sc.pl.umap(adata, color=['leiden'] + available_markers[:2], show=False)
            plt.tight_layout()
            plt.savefig("./neuro_data_warehouse/snrna_marker_expression_umap.png", dpi=300)
            print("UMAP plot saved to './neuro_data_warehouse/snrna_marker_expression_umap.png'")

if __name__ == "__main__":
    db_file = Path("./neuro_data_warehouse/neuro_warehouse.duckdb")
    analyzer = NeuroDataAnalyzer(db_path=db_file)

    cohort_df = analyzer.query_target_cohort(
        stages=["Stage IV", "Stage V", "Stage VI"],
        target_region="Substantia Nigra"
    )
    
    if not cohort_df.empty:
        combined_adata = analyzer.load_and_merge_anndata(cohort_df)
        analyzer.analyze_aggregation_drivers(combined_adata, marker_genes=["SNCA", "MAPT", "AQP4", "MBP"])
    else:
        print("No samples found matching query criteria.")

       # import duckdb

# Connect to database
conn = duckdb.connect("./neuro_data_warehouse/neuro_warehouse.duckdb")

# Inspect what brain regions actually exist in the table
print("Available Tissue Regions:")
print(conn.execute("SELECT DISTINCT tissue_region FROM snrna_seq_assays;").fetchall())

# Inspect available donors and their Braak stages
print("\nDonor Clinical Data:")
print(conn.execute("SELECT donor_id, sex, braak_stage_tau FROM donor_clinical;").fetchall())