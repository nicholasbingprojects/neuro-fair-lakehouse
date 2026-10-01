import json
import re
from pathlib import Path
import spacy
from spacy.pipeline import EntityRuler
import pandas as pd
class ClinicalDeidentifier:
    """Strips PII (names, dates, contact numbers, Medicare) using spaCy NER and regex patterns."""

    def __init__(self, spacy_model: str = "en_core_web_sm"):
        try:
            self.nlp = spacy.load(spacy_model)
        except OSError:
            raise RuntimeError(f"Run: python -m spacy download {spacy_model}")

        ruler = self.nlp.add_pipe("entity_ruler", before="ner")
        patterns = [
            {"label": "MEDICARE_NUM", "pattern": [{"TEXT": {"REGEX": r"^\d{4}\s?\d{5}\s?\d{1}$"}}]},
            {"label": "PHONE_NUM", "pattern": [{"TEXT": {"REGEX": r"^(\+61|0)[2-9]\d{8}$|^04\d{2}-\d{3}-\d{3}$"}}]},
            {"label": "DATE_STRICT", "pattern": [{"TEXT": {"REGEX": r"^\d{4}-\d{2}-\d{2}$|^\d{2}/\d{2}/\d{4}$"}}]}
        ]
        ruler.add_patterns(patterns)

        self.target_pii_labels = {
            "PERSON", "DATE", "DATE_STRICT", "GPE", "FAC", "MEDICARE_NUM", "PHONE_NUM"
        }

    def scrub_text(self, text: str) -> str:
        doc = self.nlp(text)
        entities = sorted(doc.ents, key=lambda e: e.start_char, reverse=True)
        
        scrubbed_text = text
        for ent in entities:
            if ent.label_ in self.target_pii_labels:
                placeholder = f"[{ent.label_}_REDACTED]"
                scrubbed_text = scrubbed_text[:ent.start_char] + placeholder + scrubbed_text[ent.end_char:]

        scrubbed_text = re.sub(r"Medicare:\s*\d{4}\s?\d{5}\s?\d{1}", "Medicare: [MEDICARE_REDACTED]", scrubbed_text)
        scrubbed_text = re.sub(r"DOB:\s*\d{4}-\d{2}-\d{2}", "DOB: [DOB_REDACTED]", scrubbed_text)
        return scrubbed_text

    def process_file(self, input_path: Path, output_path: Path):
        with open(input_path, "r") as f:
            records = json.load(f)

        processed_records = []
        for rec in records:
            cleaned_rec = rec.copy()
            raw_note = cleaned_rec.pop("unredacted_clinical_note", "")
            cleaned_rec["redacted_clinical_note"] = self.scrub_text(raw_note)
            cleaned_rec["deidentification_status"] = "PASSED_SPACY_PII_V1"
            processed_records.append(cleaned_rec)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(processed_records, f, indent=2)

        print(f"De-identified {len(processed_records)} clinical records -> {output_path.name}")

if __name__ == "__main__":
    base_dir = Path("./neuro_data_warehouse/clinical")
    deidentifier = ClinicalDeidentifier(spacy_model="en_core_web_sm")
    deidentifier.process_file(
        input_path=base_dir / "clinical_metadata_raw.json",
        output_path=base_dir / "clinical_metadata_deidentified.json"
    )