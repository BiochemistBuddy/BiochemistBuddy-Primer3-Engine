import hashlib
import json
from typing import Any

import primer3
from primer3 import bindings, thermoanalysis

from app.models import (
    DesignRequest,
    DesignResult,
    InternalOligoCandidate,
    InternalOligoDesignResult,
    InternalOligoThermodynamics,
    Primer,
    PrimerPair,
    ResultProvenance,
    Thermodynamics,
)


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _canonical_rule_set(rule_set: dict[str, Any]) -> str:
    return json.dumps(rule_set, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _delta_g_kcal(result: Any) -> float:
    return result.dg / 1000


def _primer(result: dict[str, Any]) -> Primer:
    return Primer(
        sequence=result["SEQUENCE"],
        coordinates=tuple(result["COORDS"]),
        melting_temperature_celsius=result["TM"],
        gc_fraction=result["GC_PERCENT"] / 100,
    )


def _thermodynamics(left_sequence: str, right_sequence: str) -> Thermodynamics:
    return Thermodynamics(
        left_hairpin_delta_g_kcal_per_mol=_delta_g_kcal(bindings.calc_hairpin(left_sequence)),
        left_homodimer_delta_g_kcal_per_mol=_delta_g_kcal(bindings.calc_homodimer(left_sequence)),
        right_hairpin_delta_g_kcal_per_mol=_delta_g_kcal(bindings.calc_hairpin(right_sequence)),
        right_homodimer_delta_g_kcal_per_mol=_delta_g_kcal(
            bindings.calc_homodimer(right_sequence)
        ),
        heterodimer_delta_g_kcal_per_mol=_delta_g_kcal(
            bindings.calc_heterodimer(left_sequence, right_sequence)
        ),
    )


def _passes_delta_g_threshold(values: Thermodynamics, threshold: float) -> bool:
    return all(value > threshold for value in values.model_dump().values())


def _primer3_arguments(request: DesignRequest, rule_set: dict[str, Any]) -> dict[str, Any]:
    if rule_set.get("workflow") != "pcr":
        raise ValueError("PCR design requires a PCR rule set")
    try:
        common = rule_set["common_constraints"]
        workflow = rule_set["workflow_constraints"]
        length = common["target_binding_length_bases"]
        temperature = common["melting_temperature_celsius"]
        gc_fraction = common["gc_fraction"]
        return {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_NUM_RETURN": request.max_results,
            "PRIMER_PRODUCT_SIZE_RANGE": [[request.product_size.min, request.product_size.max]],
            "PRIMER_MIN_SIZE": length["min"],
            "PRIMER_OPT_SIZE": length.get("preferred", length["min"]),
            "PRIMER_MAX_SIZE": length["max"],
            "PRIMER_MIN_TM": temperature["min"],
            "PRIMER_OPT_TM": temperature.get("preferred", temperature["min"]),
            "PRIMER_MAX_TM": temperature["max"],
            "PRIMER_MIN_GC": gc_fraction["min"] * 100,
            "PRIMER_OPT_GC_PERCENT": gc_fraction.get("preferred", gc_fraction["min"]) * 100,
            "PRIMER_MAX_GC": gc_fraction["max"] * 100,
            "PRIMER_MAX_POLY_X": common["max_homopolymer_bases"],
            "PRIMER_PAIR_MAX_DIFF_TM": workflow["max_pair_tm_difference_celsius"],
        }
    except (KeyError, TypeError) as error:
        raise ValueError("PCR rule set is missing required constraints") from error


def design_pcr_primers(request: DesignRequest, rule_set: dict[str, Any]) -> DesignResult:
    global_arguments = _primer3_arguments(request, rule_set)
    raw = bindings.design_primers(
        seq_args={"SEQUENCE_ID": request.sequence_id, "SEQUENCE_TEMPLATE": request.template},
        global_args=global_arguments,
    )

    threshold = rule_set["common_constraints"][
        "minimum_secondary_structure_delta_g_kcal_per_mol"
    ]
    pairs: list[PrimerPair] = []
    rejected_pair_count = 0
    for index in range(raw["PRIMER_PAIR_NUM_RETURNED"]):
        left = raw["PRIMER_LEFT"][index]
        right = raw["PRIMER_RIGHT"][index]
        pair = raw["PRIMER_PAIR"][index]
        thermodynamics = _thermodynamics(left["SEQUENCE"], right["SEQUENCE"])
        if not _passes_delta_g_threshold(thermodynamics, threshold):
            rejected_pair_count += 1
            continue
        pairs.append(
            PrimerPair(
                rank=len(pairs) + 1,
                left=Primer(
                    sequence=left["SEQUENCE"],
                    coordinates=tuple(left["COORDS"]),
                    melting_temperature_celsius=left["TM"],
                    gc_fraction=left["GC_PERCENT"] / 100,
                ),
                right=Primer(
                    sequence=right["SEQUENCE"],
                    coordinates=tuple(right["COORDS"]),
                    melting_temperature_celsius=right["TM"],
                    gc_fraction=right["GC_PERCENT"] / 100,
                ),
                product_size_bases=pair["PRODUCT_SIZE"],
                penalty=pair["PENALTY"],
                thermodynamics=thermodynamics,
            )
        )

    explanations = {
        key: raw[key]
        for key in ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN")
        if key in raw
    }
    return DesignResult(
        sequence_id=request.sequence_id,
        pairs=pairs,
        rejected_pair_count=rejected_pair_count,
        primer3_explanations=explanations,
        provenance=ResultProvenance(
            primer3_py_version=primer3.__version__,
            libprimer3_version=thermoanalysis.get_libprimer3_version(),
            sequence_sha256=_sha256_text(request.template),
            rule_set_sha256=_sha256_text(_canonical_rule_set(rule_set)),
            primer3_global_arguments=global_arguments,
        ),
    )


def _internal_oligo_arguments(
    request: DesignRequest, constraints: dict[str, Any]
) -> dict[str, Any]:
    try:
        primer = constraints["primer"]
        internal = constraints["internal_oligo"]
        return {
            "PRIMER_TASK": "generic",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 1,
            "PRIMER_NUM_RETURN": request.max_results,
            "PRIMER_PRODUCT_SIZE_RANGE": [[request.product_size.min, request.product_size.max]],
            "PRIMER_MIN_SIZE": primer["length_bases"]["min"],
            "PRIMER_OPT_SIZE": primer["length_bases"]["preferred"],
            "PRIMER_MAX_SIZE": primer["length_bases"]["max"],
            "PRIMER_MIN_TM": primer["melting_temperature_celsius"]["min"],
            "PRIMER_OPT_TM": primer["melting_temperature_celsius"]["preferred"],
            "PRIMER_MAX_TM": primer["melting_temperature_celsius"]["max"],
            "PRIMER_MIN_GC": primer["gc_fraction"]["min"] * 100,
            "PRIMER_OPT_GC_PERCENT": primer["gc_fraction"]["preferred"] * 100,
            "PRIMER_MAX_GC": primer["gc_fraction"]["max"] * 100,
            "PRIMER_MAX_POLY_X": primer["max_homopolymer_bases"],
            "PRIMER_PAIR_MAX_DIFF_TM": primer["max_pair_tm_difference_celsius"],
            "PRIMER_INTERNAL_MIN_SIZE": internal["length_bases"]["min"],
            "PRIMER_INTERNAL_OPT_SIZE": internal["length_bases"]["preferred"],
            "PRIMER_INTERNAL_MAX_SIZE": internal["length_bases"]["max"],
            "PRIMER_INTERNAL_MIN_TM": internal["melting_temperature_celsius"]["min"],
            "PRIMER_INTERNAL_OPT_TM": internal["melting_temperature_celsius"]["preferred"],
            "PRIMER_INTERNAL_MAX_TM": internal["melting_temperature_celsius"]["max"],
            "PRIMER_INTERNAL_MIN_GC": internal["gc_fraction"]["min"] * 100,
            "PRIMER_INTERNAL_OPT_GC_PERCENT": internal["gc_fraction"]["preferred"] * 100,
            "PRIMER_INTERNAL_MAX_GC": internal["gc_fraction"]["max"] * 100,
            "PRIMER_INTERNAL_MAX_POLY_X": internal["max_homopolymer_bases"],
        }
    except (KeyError, TypeError) as error:
        raise ValueError("Internal-oligo constraints are incomplete") from error


def _internal_thermodynamics(
    left_sequence: str, right_sequence: str, internal_sequence: str
) -> InternalOligoThermodynamics:
    return InternalOligoThermodynamics(
        hairpin_delta_g_kcal_per_mol=_delta_g_kcal(bindings.calc_hairpin(internal_sequence)),
        homodimer_delta_g_kcal_per_mol=_delta_g_kcal(bindings.calc_homodimer(internal_sequence)),
        left_heterodimer_delta_g_kcal_per_mol=_delta_g_kcal(
            bindings.calc_heterodimer(left_sequence, internal_sequence)
        ),
        right_heterodimer_delta_g_kcal_per_mol=_delta_g_kcal(
            bindings.calc_heterodimer(right_sequence, internal_sequence)
        ),
    )


def design_with_internal_oligo(
    request: DesignRequest, constraints: dict[str, Any]
) -> InternalOligoDesignResult:
    global_arguments = _internal_oligo_arguments(request, constraints)
    raw = bindings.design_primers(
        seq_args={"SEQUENCE_ID": request.sequence_id, "SEQUENCE_TEMPLATE": request.template},
        global_args=global_arguments,
    )
    threshold = constraints["minimum_secondary_structure_delta_g_kcal_per_mol"]
    candidates: list[InternalOligoCandidate] = []
    rejected_candidate_count = 0
    for index in range(raw["PRIMER_PAIR_NUM_RETURNED"]):
        left = raw["PRIMER_LEFT"][index]
        right = raw["PRIMER_RIGHT"][index]
        internal = raw["PRIMER_INTERNAL"][index]
        pair = raw["PRIMER_PAIR"][index]
        primer_thermodynamics = _thermodynamics(left["SEQUENCE"], right["SEQUENCE"])
        internal_thermodynamics = _internal_thermodynamics(
            left["SEQUENCE"], right["SEQUENCE"], internal["SEQUENCE"]
        )
        if not _passes_delta_g_threshold(primer_thermodynamics, threshold) or not all(
            value > threshold for value in internal_thermodynamics.model_dump().values()
        ):
            rejected_candidate_count += 1
            continue
        candidates.append(
            InternalOligoCandidate(
                rank=len(candidates) + 1,
                left=_primer(left),
                right=_primer(right),
                internal=_primer(internal),
                product_size_bases=pair["PRODUCT_SIZE"],
                penalty=pair["PENALTY"],
                primer_thermodynamics=primer_thermodynamics,
                internal_thermodynamics=internal_thermodynamics,
            )
        )
    explanations = {
        key: raw[key]
        for key in (
            "PRIMER_LEFT_EXPLAIN",
            "PRIMER_RIGHT_EXPLAIN",
            "PRIMER_INTERNAL_EXPLAIN",
            "PRIMER_PAIR_EXPLAIN",
        )
        if key in raw
    }
    return InternalOligoDesignResult(
        sequence_id=request.sequence_id,
        candidates=candidates,
        rejected_candidate_count=rejected_candidate_count,
        primer3_explanations=explanations,
        provenance=ResultProvenance(
            primer3_py_version=primer3.__version__,
            libprimer3_version=thermoanalysis.get_libprimer3_version(),
            sequence_sha256=_sha256_text(request.template),
            rule_set_sha256=_sha256_text(_canonical_rule_set(constraints)),
            primer3_global_arguments=global_arguments,
        ),
    )
