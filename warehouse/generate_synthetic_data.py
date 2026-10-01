import json
import random
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Paths setup
warehouse_dir = Path("./neuro_data_warehouse")
clinical_dir = warehouse_dir / "clinical"
snrna_dir = warehouse_dir / "omics" / "metadata"
anndata_dir = warehouse_dir / "omics" / "processed"
imaging_dir = warehouse_dir / "imaging" / "metadata"

for d in [clinical_dir, snrna_dir, anndata_dir, imaging_dir]:
    d.mkdir(parents=True, exist_ok=True)

# Generate 10 Synthetic Donors
donors = []
regions = ["Substantia Nigra", "Prefrontal Cortex", "Hippocampus"]
braak_stages = ["Stage I", "Stage II", "Stage III", "Stage IV", "Stage V", "Stage VI"]

for i in range(1, 11):
    donor_id = f"DONOR_{i:03d}"
    donors.append({
        "donor_id": donor_id,
        "age_at_death": random.randint(60, 88),
        "sex": random.choice(["Male", "Female"]),
        "post_mortem_interval_hrs": round(random.uniform(3.0, 12.0), 1),
        "braak_stage_tau": random.choice(braak_stages),
        "braak_stage_synuclein": random.choice(braak_stages),
        "unredacted_clinical_note": f"Patient {donor_id} exhibited clinical signs of neurodegeneration.",
        "redacted_clinical_note": "Patient exhibits clinical signs of neurodegeneration.",
        "deidentification_status": "DEIDENTIFIED",
        "hpo_terms": ["HP:0002129", "HP:0001300"]
    })

# Save Clinical JSON
with open(clinical_dir / "clinical_metadata_deidentified.json", "w") as f:
    json.dump(donors, f, indent=2)

# Generate snRNA Assay Metadata + Mock AnnData Files
snrna_assays = []
for idx, d in enumerate(donors):
    for r_idx, region in enumerate(regions[:2]):  # 2 regions per donor
        assay_id = f"ASSAY_RNA_{idx+1:03d}_{r_idx+1}"
        h5ad_filename = f"{assay_id}.h5ad"
        h5ad_path = anndata_dir / h5ad_filename

        # Create lightweight AnnData file for Scanpy downstream analysis
        n_obs, n_vars = 100, 50
        X = np.random.negative_binomial(5, 0.3, size=(n_obs, n_vars))
        obs = pd.DataFrame({
            "donor_id": d["donor_id"],
            "tissue_region": region,
            "cell_type": random.choices(["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"], k=n_obs)
        }, index=[f"cell_{j}" for j in range(n_obs)])
        var = pd.DataFrame(index=[f"GENE_{g}" for g in range(n_vars)])
        adata = sc.AnnData(X=X, obs=obs, var=var)
        adata.write_h5ad(h5ad_path)

        snrna_assays.append({
            "assay_id": assay_id,
            "donor_id": d["donor_id"],
            "tissue_region": region,
            "modality": "snRNA-seq",
            "minseqe": {
                "experiment_type": "10x Chromium Single Cell 3'",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "read_length_bp": 150,
                "alignment_genome": "GRCh38",
                "targeted_nuclei": 10000,
                "mean_reads_per_nucleus": 45000
            },
            "file_paths": {
                "raw_fastq_r1": f"s3://neuro-lakehouse/raw/{assay_id}_R1.fastq.gz",
                "raw_fastq_r2": f"s3://neuro-lakehouse/raw/{assay_id}_R2.fastq.gz",
                "processed_anndata": str(h5ad_path.as_posix())
            }
        })

# Save snRNA JSON
with open(snrna_dir / "snrna_seq_metadata.json", "w") as f:
    json.dump(snrna_assays, f, indent=2)

# Generate Spatial Imaging Metadata
imaging_assays = []
for idx, d in enumerate(donors):
    image_id = f"ASSAY_IMG_{idx+1:03d}"
    imaging_assays.append({
        "image_assay_id": image_id,
        "donor_id": d["donor_id"],
        "tissue_region": random.choice(regions),
        "slice_thickness_um": 10.0,
        "objective": "20x",
        "pixel_resolution_um": 0.5,
        "channels_stained": ["DAPI", "GFAP", "MAP2", "TH"],
        "file_format": "OME-TIFF",
        "file_path": f"s3://neuro-lakehouse/imaging/{image_id}.ome.tif",
        "qc_metrics": {
            "signal_to_noise_ratio": 18.5,
            "autofluorescence_subtracted": True
        }
    })

# Save Imaging JSON
with open(imaging_dir / "spatial_imaging_metadata.json", "w") as f:
    json.dump(imaging_assays, f, indent=2)

print("Synthetic data and mock AnnData objects generated successfully!")