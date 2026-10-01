from pydantic import BaseModel, Field, field_validator
from typing import List, Optional
import re

class MINSEQEMetadata(BaseModel):
    experiment_type: str
    sequencing_platform: str
    read_length_bp: int = Field(..., gt=0)
    alignment_genome: str

class BrainDonorRecord(BaseModel):
    donor_id: str
    hpo_terms: List[str]
    post_mortem_interval_hrs: float
    minseqe_metadata: Optional[MINSEQEMetadata] = None

    @field_validator('donor_id')
    def validate_anonymized_id(cls, v):
        if not re.match(r"^SUBJ_[A-F0-9]{5}$", v):
            raise ValueError("Donor ID must follow 'SUBJ_XXXXX' anonymized format")
        return v

    @field_validator('hpo_terms')
    def validate_hpo_terms(cls, terms):
        for term in terms:
            if not re.match(r"^HP:\d{7}$", term):
                raise ValueError(f"Invalid HPO term format: {term}. Must match 'HP:XXXXXXX'")
        return terms

if __name__ == "__main__":
    sample = {
        "donor_id": "SUBJ_B38C1",
        "hpo_terms": ["HP:0002129", "HP:0001300"],
        "post_mortem_interval_hrs": 4.5
    }
    validated = BrainDonorRecord(**sample)
    print(f"Pydantic schema successfully validated: {validated.donor_id}")