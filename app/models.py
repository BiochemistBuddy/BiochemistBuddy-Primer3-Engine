from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class ProductSizeRange(BaseModel):
    model_config = ConfigDict(extra="forbid")

    min: int = Field(ge=1)
    max: int = Field(ge=1)

    @model_validator(mode="after")
    def validate_order(self) -> "ProductSizeRange":
        if self.min > self.max:
            raise ValueError("product size min must not exceed max")
        return self


class DesignRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sequence_id: str = Field(min_length=1, max_length=200)
    template: str = Field(min_length=1)
    product_size: ProductSizeRange
    max_results: int = Field(default=5, ge=1, le=50)

    @field_validator("template")
    @classmethod
    def normalize_template(cls, value: str) -> str:
        normalized = "".join(value.split()).upper()
        invalid = sorted(set(normalized) - set("ACGT"))
        if invalid:
            raise ValueError(f"template contains unsupported bases: {''.join(invalid)}")
        return normalized


class DesignEnvelope(BaseModel):
    model_config = ConfigDict(extra="forbid")

    request: DesignRequest
    rule_set: dict[str, Any]


class Thermodynamics(BaseModel):
    left_hairpin_delta_g_kcal_per_mol: float
    left_homodimer_delta_g_kcal_per_mol: float
    right_hairpin_delta_g_kcal_per_mol: float
    right_homodimer_delta_g_kcal_per_mol: float
    heterodimer_delta_g_kcal_per_mol: float


class Primer(BaseModel):
    sequence: str
    coordinates: tuple[int, int]
    melting_temperature_celsius: float
    gc_fraction: float


class PrimerPair(BaseModel):
    rank: int
    left: Primer
    right: Primer
    product_size_bases: int
    penalty: float
    thermodynamics: Thermodynamics


class ResultProvenance(BaseModel):
    primer3_py_version: str
    libprimer3_version: str
    sequence_sha256: str
    rule_set_sha256: str
    primer3_global_arguments: dict[str, Any]


class DesignResult(BaseModel):
    sequence_id: str
    pairs: list[PrimerPair]
    rejected_pair_count: int
    primer3_explanations: dict[str, str]
    provenance: ResultProvenance
