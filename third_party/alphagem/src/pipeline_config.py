"""Versioned AlphaGEM parameter profiles.

The paper profile is based on the published Methods and Supplementary Table 1.
The GitHub profile preserves the public main-branch defaults at commit a600b7e
so that repaired paper runs can always be compared with the released code.
"""

from dataclasses import dataclass, replace
from typing import Dict


@dataclass(frozen=True)
class PLMSearchThresholds:
    ss_predictor: float
    coverage: float
    identity: int
    bbh_coverage: float
    bbh_identity: int
    direction: str = "bidirect"
    bbh_backend: str = "diamond"


@dataclass(frozen=True)
class StructureThresholds:
    usalign: float
    tm_score: float
    coverage: float
    plddt: float
    transporter_tm_score: float
    transporter_coverage: float
    high_confidence_tm_score: float
    high_confidence_coverage: float
    dbscan_eps: float


@dataclass(frozen=True)
class PipelineProfile:
    name: str
    method_contract: str
    plmsearch: PLMSearchThresholds
    structure: StructureThresholds
    prottrans_distance: float
    rhea_ss_predictor: float
    ensemble_score: float
    gapfill_strategy: str
    gapfill_search_backend: str
    gapfill_bitscore_threshold: float
    gapfill_bitscore_operator: str
    embedding_semantics: str
    eggnog_database_mode: str

    @property
    def is_paper_method(self) -> bool:
        return self.method_contract == "paper"

PAPER_PROFILE = PipelineProfile(
    name="paper",
    method_contract="paper",
    plmsearch=PLMSearchThresholds(
        ss_predictor=0.8,
        coverage=0.6,
        identity=30,
        bbh_coverage=0.5,
        bbh_identity=60,
        bbh_backend="diamond",
    ),
    structure=StructureThresholds(
        usalign=0.5,
        tm_score=0.7,
        coverage=0.8,
        plddt=70,
        transporter_tm_score=0.7,
        transporter_coverage=0.8,
        high_confidence_tm_score=0.9,
        high_confidence_coverage=0.9,
        dbscan_eps=1.0,
    ),
    prottrans_distance=1.1,
    rhea_ss_predictor=0.9,
    ensemble_score=0.4,
    gapfill_strategy="selected_medium",
    gapfill_search_backend="tblastn_genome",
    gapfill_bitscore_threshold=200.0,
    gapfill_bitscore_operator="gt",
    embedding_semantics="corrected",
    eggnog_database_mode="paper_taxon",
)


GITHUB_CURRENT_PROFILE = PipelineProfile(
    name="github_current",
    method_contract="released_code",
    plmsearch=PLMSearchThresholds(
        ss_predictor=0.95,
        coverage=0.9,
        identity=70,
        bbh_coverage=0.8,
        bbh_identity=70,
        direction="direct",
        bbh_backend="blastp",
    ),
    structure=StructureThresholds(
        usalign=0.5,
        tm_score=0.7,
        coverage=0.8,
        plddt=70,
        transporter_tm_score=0.7,
        transporter_coverage=0.8,
        high_confidence_tm_score=0.9,
        high_confidence_coverage=0.95,
        dbscan_eps=1.0,
    ),
    prottrans_distance=1.2,
    rhea_ss_predictor=0.9,
    ensemble_score=0.4,
    gapfill_strategy="legacy_dual_medium",
    gapfill_search_backend="legacy_auto",
    gapfill_bitscore_threshold=200.0,
    gapfill_bitscore_operator="ge",
    embedding_semantics="legacy_repo",
    eggnog_database_mode="packaged_fallback",
)


LEGACY_MAIN_PROFILE = PipelineProfile(
    name="legacy_main",
    method_contract="released_code",
    plmsearch=PLMSearchThresholds(
        ss_predictor=0.8,
        coverage=0.5,
        identity=20,
        bbh_coverage=0.25,
        bbh_identity=60,
        direction="direct",
        bbh_backend="blastp",
    ),
    structure=StructureThresholds(
        usalign=0.5,
        tm_score=0.7,
        coverage=0.8,
        plddt=70,
        transporter_tm_score=0.7,
        transporter_coverage=0.8,
        high_confidence_tm_score=0.9,
        high_confidence_coverage=0.9,
        dbscan_eps=0.8,
    ),
    prottrans_distance=1.2,
    rhea_ss_predictor=0.9,
    ensemble_score=0.4,
    gapfill_strategy="legacy_dual_medium",
    gapfill_search_backend="legacy_auto",
    gapfill_bitscore_threshold=200.0,
    gapfill_bitscore_operator="ge",
    embedding_semantics="legacy_repo",
    eggnog_database_mode="packaged_fallback",
)

AUTHOR_ARTIFACT_PROFILE = PipelineProfile(
    name="author_artifact_legacy",
    method_contract="released_code",
    plmsearch=PLMSearchThresholds(
        ss_predictor=0.8,
        coverage=0.6,
        identity=30,
        bbh_coverage=0.5,
        bbh_identity=60,
        direction="direct",
        bbh_backend="blastp",
    ),
    structure=StructureThresholds(
        usalign=0.5,
        tm_score=0.7,
        coverage=0.8,
        plddt=70,
        transporter_tm_score=0.7,
        transporter_coverage=0.8,
        high_confidence_tm_score=0.9,
        high_confidence_coverage=0.9,
        dbscan_eps=1.0,
    ),
    prottrans_distance=1.2,
    rhea_ss_predictor=0.9,
    ensemble_score=0.4,
    gapfill_strategy="legacy_dual_medium",
    gapfill_search_backend="legacy_auto",
    gapfill_bitscore_threshold=200.0,
    gapfill_bitscore_operator="ge",
    embedding_semantics="legacy_repo",
    eggnog_database_mode="packaged_fallback",
)


PROFILES: Dict[str, PipelineProfile] = {
    profile.name: profile
    for profile in (
        PAPER_PROFILE,
        AUTHOR_ARTIFACT_PROFILE,
        GITHUB_CURRENT_PROFILE,
        LEGACY_MAIN_PROFILE,
    )
}

PAPER_STRUCTURE_BY_DOMAIN = {
    "mammal": (0.7, 0.85, 0.9, 70),
    "fungus": (0.5, 0.7, 0.8, 70),
    "bacterium": (0.5, 0.65, 0.75, 70),
}

REFERENCE_DOMAINS = {
    "human": "mammal",
    "yeast": "fungus",
    "ecoli": "bacterium",
    "strco": "bacterium",
    "synechocystis": "bacterium",
}


def get_profile(name: str) -> PipelineProfile:
    try:
        return PROFILES[name]
    except KeyError as exc:
        choices = ", ".join(sorted(PROFILES))
        raise ValueError(f"Unknown profile {name!r}; choose one of: {choices}") from exc


def resolve_profile_for_reference(
    profile: PipelineProfile, refname: str
) -> PipelineProfile:
    """Resolve Supplementary Table 1 thresholds before CLI overrides."""
    if not profile.is_paper_method:
        return profile
    try:
        domain = REFERENCE_DOMAINS[refname]
    except KeyError as exc:
        raise ValueError(f"Unknown reference domain for {refname!r}") from exc
    usalign, tm_score, coverage, plddt = PAPER_STRUCTURE_BY_DOMAIN[domain]
    structure = replace(
        profile.structure,
        usalign=usalign,
        tm_score=tm_score,
        coverage=coverage,
        plddt=plddt,
        transporter_tm_score=tm_score,
        transporter_coverage=coverage,
    )
    return replace(profile, structure=structure)


def override_profile(
    profile: PipelineProfile,
    *,
    ss_predictor=None,
    plm_coverage=None,
    plm_identity=None,
    bbh_coverage=None,
    bbh_identity=None,
    usalign=None,
    tm_score=None,
    structure_coverage=None,
    plddt=None,
    transporter_tm_score=None,
    transporter_coverage=None,
    high_confidence_tm_score=None,
    high_confidence_coverage=None,
    dbscan_eps=None,
    prottrans_distance=None,
    rhea_ss_predictor=None,
    ensemble_score=None,
) -> PipelineProfile:
    """Return a derived profile without mutating a named reproducibility profile."""

    override_values = (
        ss_predictor,
        plm_coverage,
        plm_identity,
        bbh_coverage,
        bbh_identity,
        usalign,
        tm_score,
        structure_coverage,
        plddt,
        transporter_tm_score,
        transporter_coverage,
        high_confidence_tm_score,
        high_confidence_coverage,
        dbscan_eps,
        prottrans_distance,
        rhea_ss_predictor,
        ensemble_score,
    )
    if all(value is None for value in override_values):
        return profile

    plmsearch = replace(
        profile.plmsearch,
        **{
            key: value
            for key, value in {
                "ss_predictor": ss_predictor,
                "coverage": plm_coverage,
                "identity": plm_identity,
                "bbh_coverage": bbh_coverage,
                "bbh_identity": bbh_identity,
            }.items()
            if value is not None
        },
    )
    structure = replace(
        profile.structure,
        **{
            key: value
            for key, value in {
                "usalign": usalign,
                "tm_score": tm_score,
                "coverage": structure_coverage,
                "plddt": plddt,
                "transporter_tm_score": transporter_tm_score,
                "transporter_coverage": transporter_coverage,
                "high_confidence_tm_score": high_confidence_tm_score,
                "high_confidence_coverage": high_confidence_coverage,
                "dbscan_eps": dbscan_eps,
            }.items()
            if value is not None
        },
    )
    return replace(
        profile,
        name=f"{profile.name}+overrides",
        plmsearch=plmsearch,
        structure=structure,
        **{
            key: value
            for key, value in {
                "prottrans_distance": prottrans_distance,
                "rhea_ss_predictor": rhea_ss_predictor,
                "ensemble_score": ensemble_score,
            }.items()
            if value is not None
        },
    )
