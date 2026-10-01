import os
import json
import random
import hashlib
from pathlib import Path
import pandas as pd
from faker import Faker

fake = Faker('en_AU')
Faker.seed(42)
random.seed(42)

class SyntheticNeuroDataGenerator:
    # Generates synthetic multi-modal neuroscience datasets for testing FAIR ingestion and governance pipelines."""

    def __init__(self, num_donors=15):
        self.num_donors = num_donors
        self.hpo_catalog = {
            "HP:0002129": "Parkinsonism",
            "HP:0001300": "Action tremor",
            "HP:0000726": "Dementia",
            "HP:0002354": "Memory impairment",
            "HP:0002063": "Rigidity",
            "HP:0002172": "Postural instability"
        }
        self.brain_regions = [
            "Substantia Nigra", "Dorsolateral Prefrontal Cortex", 
            "Hippocampus (CA1)", "Caudate Nucleus"
        ]
        self.imaging_targets = ["alpha-synuclein (MJFR1)", "p-tau (AT8)", "NeuN", "Iba1", "DAPI"]

    def _generate_salt_donor_id(self, idx: int) -> str:
        salt = f"hrec_approved_salt_{idx}".encode('utf-8')
        return f"SUBJ_{hashlib.sha256(salt).hexdigest()[:5].upper()}"

    def generate_clinical_metadata(self) -> pd.DataFrame:
        clinical_data = []
        for i in range(self.num_donors):
            donor_id = self._generate_salt_donor_id(i)
            age = random.randint(60, 90)
            
            unredacted_note = (
                f"Patient {fake.name()} (DOB: {fake.date_of_birth(minimum_age=age, maximum_age=age)}, "
                f"Medicare: {fake.bothify(text='#### ##### #')}) was reviewed by Dr. {fake.last_name()} "
                f"at {fake.city()} Brain Clinic. Consent signed on {fake.date_this_decade()}."
            )

            clinical_data.append({
                "donor_id": donor_id,
                "age_at_death": age,
                "sex": random.choice(["Male", "Female"]),
                "post_mortem_interval_hrs": round(random.uniform(2.0, 18.5), 1),
                "hpo_terms": random.sample(list(self.hpo_catalog.keys()), k=random.randint(1, 3)),
                "braak_stage_tau": f"Stage {random.choice(['I', 'II', 'III', 'IV', 'V', 'VI'])}",
                "braak_stage_synuclein": f"Stage {random.randint(0, 6)}",
                "unredacted_clinical_note": unredacted_note
            })
        return pd.DataFrame(clinical_data)

    def generate_snrna_metadata(self, clinical_df: pd.DataFrame) -> list:
        snrna_assays = []
        for _, row in clinical_df.iterrows():
            donor_id = row['donor_id']
            for region in random.sample(self.brain_regions, k=2):
                assay_id = f"ASSAY_snRNA_{donor_id}_{region.replace(' ', '_')[:6]}"
                snrna_assays.append({
                    "assay_id": assay_id,
                    "donor_id": donor_id,
                    "tissue_region": region,
                    "modality": "snRNA-seq",
                    "minseqe": {
                        "experiment_type": "single-nucleus RNA sequencing",
                        "sequencing_platform": "Illumina NovaSeq 6000",
                        "read_length_bp": 150,
                        "alignment_genome": "GRCh38",
                        "targeted_nuclei": random.randint(6000, 12000),
                        "mean_reads_per_nucleus": random.randint(25000, 50000)
                    },
                    "file_paths": {
                        "raw_fastq_r1": f"/hpc/data/raw/snrna/{donor_id}/{assay_id}_R1.fastq.gz",
                        "raw_fastq_r2": f"/hpc/data/raw/snrna/{donor_id}/{assay_id}_R2.fastq.gz",
                        "processed_anndata": f"/hpc/data/processed/snrna/{donor_id}/{assay_id}_filtered.h5ad"
                    }
                })
        return snrna_assays

    def generate_spatial_imaging_metadata(self, clinical_df: pd.DataFrame) -> list:
        imaging_assays = []
        for _, row in clinical_df.iterrows():
            donor_id = row['donor_id']
            region = random.choice(self.brain_regions)
            image_id = f"IMG_OME_{donor_id}_{region.replace(' ', '_')[:6]}"
            imaging_assays.append({
                "image_assay_id": image_id,
                "donor_id": donor_id,
                "tissue_region": region,
                "slice_thickness_um": 10.0,
                "objective": "40x Oil",
                "pixel_resolution_um": 0.325,
                "channels_stained": self.imaging_targets,
                "file_format": "OME-TIFF",
                "file_path": f"/hpc/data/imaging/ome_tiff/{donor_id}/{image_id}.ome.tif",
                "qc_metrics": {
                    "signal_to_noise_ratio": round(random.uniform(12.5, 34.0), 2),
                    "autofluorescence_subtracted": True
                }
            })
        return imaging_assays

if __name__ == "__main__":
    base_dir = Path("./neuro_data_warehouse")
    clinical_dir = base_dir / "clinical"
    imaging_dir = base_dir / "imaging" / "metadata"
    omics_dir = base_dir / "omics" / "metadata"

    for folder in [clinical_dir, imaging_dir, omics_dir]:
        folder.mkdir(parents=True, exist_ok=True)

    generator = SyntheticNeuroDataGenerator(num_donors=10)
    df_clinical = generator.generate_clinical_metadata()
    snrna_metadata = generator.generate_snrna_metadata(df_clinical)
    imaging_metadata = generator.generate_spatial_imaging_metadata(df_clinical)

    df_clinical.to_json(clinical_dir / "clinical_metadata_raw.json", orient="records", indent=2)
    with open(omics_dir / "snrna_seq_metadata.json", "w") as f:
        json.dump(snrna_metadata, f, indent=2)
    with open(imaging_dir / "spatial_imaging_metadata.json", "w") as f:
        json.dump(imaging_metadata, f, indent=2)

    print(f"Synthetic raw dataset exported to '{base_dir.resolve()}'.")