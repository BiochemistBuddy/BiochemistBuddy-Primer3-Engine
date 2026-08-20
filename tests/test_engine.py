import random

from app.engine import design_pcr_primers
from app.models import DesignRequest


def _rule_set() -> dict:
    return {
        "workflow": "pcr",
        "common_constraints": {
            "target_binding_length_bases": {"min": 18, "max": 30, "preferred": 22},
            "melting_temperature_celsius": {"min": 60, "max": 64, "preferred": 62},
            "gc_fraction": {"min": 0.35, "max": 0.65, "preferred": 0.5},
            "max_homopolymer_bases": 3,
            "minimum_secondary_structure_delta_g_kcal_per_mol": -9,
        },
        "workflow_constraints": {"max_pair_tm_difference_celsius": 3},
    }


def test_design_is_reproducible_and_records_primer3_provenance() -> None:
    generator = random.Random(7)
    template = "".join(generator.choice("ACGT") for _ in range(800))
    request = DesignRequest(
        sequence_id="synthetic-fixture",
        template=template,
        product_size={"min": 100, "max": 300},
        max_results=2,
    )

    first = design_pcr_primers(request, _rule_set())
    second = design_pcr_primers(request, _rule_set())

    assert first == second
    assert len(first.pairs) == 2
    assert first.provenance.primer3_py_version == "2.3.0"
    assert first.provenance.primer3_global_arguments["PRIMER_PRODUCT_SIZE_RANGE"] == [[100, 300]]
