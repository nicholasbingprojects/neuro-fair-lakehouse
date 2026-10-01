import duckdb
from pathlib import Path
import json

class NeuroWarehouseIngestor:
    """Ingests metadata into an embedded DuckDB database."""

    def __init__(self, db_path: Path):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = duckdb.connect(str(self.db_path))
        self._init_schema()

    def _init_schema(self):
        self.conn.execute("""
        CREATE TABLE IF NOT EXISTS donor_clinical (
            donor_id VARCHAR PRIMARY KEY,
            age_at_death INTEGER CHECK (age_at_death <= 89),
            sex VARCHAR(10),
            post_mortem_interval_hrs DOUBLE,
            braak_stage_tau VARCHAR(10),
            braak_stage_synuclein VARCHAR(10),
            redacted_clinical_note TEXT,
            deidentification_status VARCHAR(30)
        );

        CREATE TABLE IF NOT EXISTS donor_hpo_terms (
            donor_id VARCHAR REFERENCES donor_clinical(donor_id),
            hpo_term_id VARCHAR(12),
            PRIMARY KEY (donor_id, hpo_term_id)
        );

        CREATE TABLE IF NOT EXISTS snrna_seq_assays (
            assay_id VARCHAR PRIMARY KEY,
            donor_id VARCHAR REFERENCES donor_clinical(donor_id),
            tissue_region VARCHAR(100),
            modality VARCHAR(30),
            experiment_type VARCHAR(100),
            sequencing_platform VARCHAR(50),
            read_length_bp INTEGER,
            alignment_genome VARCHAR(20),
            targeted_nuclei INTEGER,
            mean_reads_per_nucleus INTEGER,
            raw_fastq_r1 VARCHAR(255),
            raw_fastq_r2 VARCHAR(255),
            processed_anndata VARCHAR(255)
        );

        CREATE TABLE IF NOT EXISTS spatial_imaging_assays (
            image_assay_id VARCHAR PRIMARY KEY,
            donor_id VARCHAR REFERENCES donor_clinical(donor_id),
            tissue_region VARCHAR(100),
            slice_thickness_um DOUBLE,
            objective VARCHAR(20),
            pixel_resolution_um DOUBLE,
            channels_stained VARCHAR[],
            file_format VARCHAR(20),
            file_path VARCHAR(255),
            signal_to_noise_ratio DOUBLE,
            autofluorescence_subtracted BOOLEAN
        );
        """)

    def ingest_clinical_data(self, json_path: Path):
        safe_path = Path(json_path).as_posix()
        if not Path(json_path).exists():
            raise FileNotFoundError(f"Missing clinical metadata file at: {json_path}")

        self.conn.execute(f"""
            INSERT OR REPLACE INTO donor_clinical
            SELECT 
                donor_id, 
                CAST(age_at_death AS INT), 
                sex, 
                CAST(post_mortem_interval_hrs AS DOUBLE),
                braak_stage_tau, 
                braak_stage_synuclein, 
                COALESCE(
                    redacted_clinical_note,
                    unredacted_clinical_note
                ) AS redacted_clinical_note,
                COALESCE(
                    deidentification_status,
                    'RAW_PENDING_REDACTION'
                ) AS deidentification_status
            FROM read_json(
                '{safe_path}', 
                columns={{
                    'donor_id': 'VARCHAR',
                    'age_at_death': 'UBIGINT',
                    'sex': 'VARCHAR',
                    'post_mortem_interval_hrs': 'DOUBLE',
                    'braak_stage_tau': 'VARCHAR',
                    'braak_stage_synuclein': 'VARCHAR',
                    'redacted_clinical_note': 'VARCHAR',
                    'unredacted_clinical_note': 'VARCHAR',
                    'deidentification_status': 'VARCHAR'
                }}
            );
        """)

        self.conn.execute(f"""
            INSERT OR REPLACE INTO donor_hpo_terms
            SELECT 
                donor_id, 
                unnest(hpo_terms) AS hpo_term_id
            FROM read_json(
                '{safe_path}',
                columns={{
                    'donor_id': 'VARCHAR',
                    'hpo_terms': 'VARCHAR[]'
                }}
            );
        """)

    def ingest_snrna_metadata(self, json_path: Path):
        safe_path = Path(json_path).as_posix()
        if not Path(json_path).exists():
            raise FileNotFoundError(f"Missing snRNA metadata file at: {json_path}")

        self.conn.execute(f"""
            INSERT OR REPLACE INTO snrna_seq_assays
            SELECT 
                assay_id, 
                donor_id, 
                tissue_region, 
                modality,
                minseqe.experiment_type, 
                minseqe.sequencing_platform,
                CAST(minseqe.read_length_bp AS INT), 
                minseqe.alignment_genome,
                CAST(minseqe.targeted_nuclei AS INT), 
                CAST(minseqe.mean_reads_per_nucleus AS INT),
                file_paths.raw_fastq_r1, 
                file_paths.raw_fastq_r2, 
                file_paths.processed_anndata
            FROM read_json(
                '{safe_path}',
                columns={{
                    'assay_id': 'VARCHAR',
                    'donor_id': 'VARCHAR',
                    'tissue_region': 'VARCHAR',
                    'modality': 'VARCHAR',
                    'minseqe': 'STRUCT(experiment_type VARCHAR, sequencing_platform VARCHAR, read_length_bp UBIGINT, alignment_genome VARCHAR, targeted_nuclei UBIGINT, mean_reads_per_nucleus UBIGINT)',
                    'file_paths': 'STRUCT(raw_fastq_r1 VARCHAR, raw_fastq_r2 VARCHAR, processed_anndata VARCHAR)'
                }}
            );
        """)

    def ingest_spatial_imaging_metadata(self, json_path: Path):
        safe_path = Path(json_path).as_posix()
        if not Path(json_path).exists():
            raise FileNotFoundError(f"Missing imaging metadata file at: {json_path}")

        self.conn.execute(f"""
            INSERT OR REPLACE INTO spatial_imaging_assays
            SELECT 
                image_assay_id, 
                donor_id, 
                tissue_region, 
                CAST(slice_thickness_um AS DOUBLE), 
                objective,
                CAST(pixel_resolution_um AS DOUBLE), 
                channels_stained, 
                file_format, 
                file_path,
                CAST(qc_metrics.signal_to_noise_ratio AS DOUBLE),
                CAST(qc_metrics.autofluorescence_subtracted AS BOOLEAN)
            FROM read_json(
                '{safe_path}',
                columns={{
                    'image_assay_id': 'VARCHAR',
                    'donor_id': 'VARCHAR',
                    'tissue_region': 'VARCHAR',
                    'slice_thickness_um': 'DOUBLE',
                    'objective': 'VARCHAR',
                    'pixel_resolution_um': 'DOUBLE',
                    'channels_stained': 'VARCHAR[]',
                    'file_format': 'VARCHAR',
                    'file_path': 'VARCHAR',
                    'qc_metrics': 'STRUCT(signal_to_noise_ratio DOUBLE, autofluorescence_subtracted BOOLEAN)'
                }}
            );
        """)

    def verify_counts(self):
        """Prints insertion row counts for validation."""
        tables = ["donor_clinical", "donor_hpo_terms", "snrna_seq_assays", "spatial_imaging_assays"]
        print("\n--- Ingested Record Summary ---")
        for tbl in tables:
            cnt = self.conn.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
            print(f"Table '{tbl}': {cnt} rows")

    def close(self):
        self.conn.close()

if __name__ == "__main__":
    warehouse_dir = Path("./neuro_data_warehouse")
    db_file = warehouse_dir / "neuro_warehouse.duckdb"

    ingestor = NeuroWarehouseIngestor(db_path=db_file)

    # Smart file selection: prefer deidentified file if available, fall back to raw
    deidentified_path = warehouse_dir / "clinical" / "clinical_metadata_deidentified.json"
    raw_path = warehouse_dir / "clinical" / "clinical_metadata_raw.json"
    clinical_file = deidentified_path if deidentified_path.exists() else raw_path

    ingestor.ingest_clinical_data(clinical_file)
    ingestor.ingest_snrna_metadata(warehouse_dir / "omics" / "metadata" / "snrna_seq_metadata.json")
    ingestor.ingest_spatial_imaging_metadata(warehouse_dir / "imaging" / "metadata" / "spatial_imaging_metadata.json")
    
    # Print record count report
    ingestor.verify_counts()
    ingestor.close()

    print(f"\nDatabase build complete: '{db_file.resolve()}'.")