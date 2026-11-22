"""Canonical tool and derived-resource identities from the author archive.

The archive identity is the SHA-256 shipped with the AlphaGEM repository.
Bundle digests were computed from regular files in the verified archive.  Each
digest hashes sorted lines of the form
``<file_sha256>  <repo-relative-posix-path>\n``.

``paper_declared_version`` records what the paper says and is ``None`` when
the paper does not declare that resource-layer version.  It does not imply
that every bundled project exposes machine-readable release metadata.
"""

PAPER_RESOURCE_ARCHIVE_SHA256 = (
    "1ef308718a647ce339fe1423aac573489567f3c0d3e2f155e42fa1a954bd5875"
)

PAPER_METHOD_EVIDENCE_CONTRACT = {
    "claim": "paper_declared_method_contract",
    "paper_declared": {
        "cobrapy_version": "0.29.1",
        "solver_backend": None,
        "pipeline_python_version": None,
        "statistical_analysis_python_version": "3.10.15",
        "statistical_analysis_python_scope_only": True,
    },
    "author_environment": {
        "source": "AlphaGEM.yml",
        "python_version": "3.10.15",
        "cobrapy_version": "0.29.1",
        "solver_dependency": "gurobi",
        "solver_dependency_version": "11.0.3",
        "selected_solver_backend": None,
    },
    "evidence_boundary": (
        "The paper names COBRApy 0.29.1 but does not name a solver backend "
        "or a Python version for the AlphaGEM reconstruction pipeline. "
        "Python 3.10.15 is stated only for statistical analysis. The Python "
        "and Gurobi dependency pins are author-environment evidence, not "
        "paper-declared method requirements or proof that Gurobi was selected."
    ),
    "author_released_resource_boundaries": {
        "rhea_embeddings": (
            "Bitwise reuse of the six precomputed author-archive artifacts "
            "is verified. Their upstream UniProt snapshot, embedding-model "
            "identity, pooling semantics, and sharding rule were not "
            "disclosed, so exact regeneration from released inputs is not "
            "claimed."
        ),
    },
}

PAPER_EGGNOG_DATABASE_CONTRACT = {
    "path": "tools/eggnog_mapper/data/eggnog.db",
    "expected_database_version": "5.0.2",
    "bytes": 41370988544,
    "sha256": (
        "4d6857c129423889edcd607995b28d97944787a588224f31dd7f745fe9ef1e4f"
    ),
    "official_download_url": (
        "http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz"
    ),
    "identity_source": "locally_verified_materialized_database",
    "download_provenance": "original_download_chain_not_retained",
}

PAPER_EGGNOG_TAXON_DATABASE_CONTRACTS = {
    "yeast": {
        "filename": "fungi.dmnd",
        "database_role": "eggNOG_5.0_Fungi_taxonomic_proteome_database",
        "diamond_version": "2.1.10",
        "diamond_format_version": 3,
        "bytes": 879076733,
        "sha256": (
            "b9c737213a3713af370dbc03e9d798ea8ef9a576754d5707d7a7c75ec9546c16"
        ),
        "sequences": 1774481,
        "letters": 814424912,
        "source_urls": {
            "proteomes": (
                "http://eggnog5.embl.de/download/eggnog_5.0/"
                "e5.proteomes.faa"
            ),
            "taxid_info": (
                "http://eggnog5.embl.de/download/eggnog_5.0/"
                "e5.taxid_info.tsv"
            ),
            "upstream_create_dbs": (
                "https://raw.githubusercontent.com/eggnogdb/"
                "eggnog-mapper/2.1.12/create_dbs.py"
            ),
        },
        "source_sha256": {
            "proteomes": (
                "efc3685ebcbfb1d73359231dbb823cd4cf8bfc80c1ca1b98c7117d1c4674b8f2"
            ),
            "taxid_info": (
                "c3839fada0046a1c83d2eda6625a2c806931e3e982b84c43b7147dbe8474a125"
            ),
            "selected_uncompressed_fasta": (
                "0f64c4114fbcd4bbd78bb7493854226c0baeb5f5dfc9fe6768da8e41a0d7da21"
            ),
            "selected_gzip_fasta": (
                "5d20ecbbe32efe644f0b564d485097c8adddf4c6c04dfe115128b48294f829d3"
            ),
        },
        "selection_contract": (
            "eggnog-mapper 2.1.12 create_dbs.py: retain taxids whose "
            "Named Lineage contains the exact token Fungi, then retain "
            "e5.proteomes.faa records whose header taxid precedes the first dot"
        ),
        "identity_source": (
            "linux_x86_64_diamond_2_1_10_build_from_verified_selected_fasta"
        ),
    },
}

PAPER_RHEA_ARTIFACT_CONTRACT = {
    "identity_source": "verified_author_data_source_archive",
    "source_archive_sha256": PAPER_RESOURCE_ARCHIVE_SHA256,
    "paper_declared_version": None,
    "upstream_uniprot_snapshot": "unknown_not_disclosed",
    "embedding_model_identity": "unknown_not_disclosed",
    "pooling_semantics": "unknown_not_disclosed",
    "sharding_rule": "unknown_not_disclosed",
    "exact_regeneration_from_released_inputs": False,
    "evidence_boundary": (
        "Bitwise reuse of the author-provided derived artifacts is verified; "
        "their upstream generation was not reproducible from the released "
        "inputs."
    ),
    "embedding_set_audit": {
        "total_key_count": 231709,
        "duplicate_key_count": 0,
        "mapping_unique_accession_count": 231709,
        "embedding_only_accession_count": 0,
        "mapping_only_accession_count": 0,
        "embedding_accession_extraction": "pipe_delimited_field_2",
        "mapping_accession_column": "ID",
        "tensor_shape_dtype_audit": "first_tensor_per_shard_only",
    },
    "resources": {
        "rhea_uniprot_mapping": {
            "path": "rhea/rhea2uniprot_sprot.tsv",
            "archive_member_path": "rhea/rhea2uniprot_sprot.tsv",
            "artifact_role": "author_provided_rhea_uniprot_mapping",
            "bytes": 8104031,
            "sha256": (
                "b5e3927284ae0f6fed10a68c7709d75f8445043bf12d0dfa11cc2b71739c78ef"
            ),
            "audited_data_row_count": 367662,
            "audited_unique_accession_count": 231709,
        },
        "rhea_embedding_1": {
            "path": "rhea/rhea_protembedding1.pkl",
            "archive_member_path": "rhea/rhea_protembedding1.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 247866034,
            "sha256": (
                "cefab43d24d2b0540e83edae34d7f399c53531ff296d8f058de48ad6b8ff480d"
            ),
            "audited_key_count": 44784,
        },
        "rhea_embedding_2": {
            "path": "rhea/rhea_protembedding2.pkl",
            "archive_member_path": "rhea/rhea_protembedding2.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 244121600,
            "sha256": (
                "0d4cd5c8e27a323f650f3bc464623c006337f610ac2545f20cac9b7cf4e9dcd4"
            ),
            "audited_key_count": 44059,
        },
        "rhea_embedding_3": {
            "path": "rhea/rhea_protembedding3.pkl",
            "archive_member_path": "rhea/rhea_protembedding3.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 290133188,
            "sha256": (
                "68caa6f3daeaa4d76ba68962f949da89c5b1227c7e6a11bf8f72409d6bdcac4b"
            ),
            "audited_key_count": 52359,
        },
        "rhea_embedding_4": {
            "path": "rhea/rhea_protembedding4.pkl",
            "archive_member_path": "rhea/rhea_protembedding4.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 206795600,
            "sha256": (
                "0e559f8ec29be6d87c2f989a960683c15dd498948262375b4a4568d1800b3ff4"
            ),
            "audited_key_count": 37325,
        },
        "rhea_embedding_5": {
            "path": "rhea/rhea_protembedding5.pkl",
            "archive_member_path": "rhea/rhea_protembedding5.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 225585252,
            "sha256": (
                "e91071ad0e09ccbe21d85a620c7985f5e31b456cd83769c45d52db829974cea3"
            ),
            "audited_key_count": 40745,
        },
        "rhea_embedding_6": {
            "path": "rhea/rhea_protembedding6.pkl",
            "archive_member_path": "rhea/rhea_protembedding6.pkl",
            "artifact_role": "author_provided_precomputed_derived_embedding",
            "bytes": 68877921,
            "sha256": (
                "2febbcbaa297dba9b737f8edfab9c93062476292e94cdbc34ab7cf8d2c9c60fa"
            ),
            "audited_key_count": 12437,
        },
    },
}

PAPER_TOOL_BUNDLE_CONTRACTS = {
    "orthofinder": {
        "path": "tools/OrthoFinder",
        "paper_declared_version": "3.0.1b1",
        "excluded_relative_roots": (),
        "file_count": 4250,
        "bytes": 147374050,
        "sha256": (
            "797a04e48b83357812106dd306ac90519c3cb7db5c13933d407b368aa0b30a62"
        ),
    },
    "usalign": {
        "path": "tools/USalign/USalign",
        "paper_declared_version": "20241108",
        "excluded_relative_roots": (),
        "file_count": 1,
        "bytes": 3705360,
        "sha256": (
            "3c9ce64ad76a11cc599ce7768fb9b8e66df8ff302998cb439a11d7df84576840"
        ),
    },
    "eggnog_mapper": {
        "path": "tools/eggnog_mapper",
        "paper_declared_version": "2.1.12",
        "excluded_relative_roots": ("data",),
        "file_count": 227,
        "bytes": 190298115,
        "sha256": (
            "cec52ef7b1b2c317abf7b6c1e8fed8ae10beedcd1af24ead2bb54b2b1c22335f"
        ),
    },
    "eggnog_static_data": {
        "path": "tools/eggnog_mapper/data",
        "paper_declared_version": None,
        "tool_expected_database_version": "5.0.2",
        "excluded_relative_roots": ("eggnog.db", "eggnog.db.gz", "fungi.dmnd"),
        "file_count": 9,
        "bytes": 16224797454,
        "sha256": (
            "2128f9ae1e5cab076409822c04fc5b83f9c75d977513c87f92a098b6a1edf5f7"
        ),
    },
    "clean": {
        "path": "tools/CLEAN/app",
        "paper_declared_version": "1.0.1",
        "excluded_relative_roots": (
            "data/inputs",
            "results",
            "data/esm_data",
            "build",
            "dist",
            "src/CLEAN.egg-info",
        ),
        "file_count": 563,
        "bytes": 351497433,
        "sha256": (
            "34767e3a698d453775d4043872533a2d84a6b70162468ead1b26d42bcf07c52c"
        ),
    },
    "deepec_transformer": {
        "path": "tools/DeepProZyme",
        "paper_declared_version": "1.0",
        "excluded_relative_roots": (),
        "file_count": 32,
        "bytes": 1808103255,
        "sha256": (
            "9fca01c3907766b469b9e163c3301883f63e305cbbe45bda659ea6f26b72a4be"
        ),
    },
    "plmsearch": {
        "path": "plmsearchtools",
        "paper_declared_version": "2024.05.08",
        "excluded_relative_roots": (),
        "file_count": 23,
        "bytes": 118757175,
        "sha256": (
            "7734b4faa238f1a50aab2016a65a071ab34527c19c268805365a1e92fbba6a15"
        ),
    },
}

# setup.sh installs CLEAN from an isolated copy and removes packaging/build
# artifacts. The runtime tree is therefore not byte-identical to the verified
# author archive. Keep both identities instead of rewriting archive evidence.
PAPER_RUNTIME_TOOL_BUNDLE_CONTRACTS = {
    name: dict(contract)
    for name, contract in PAPER_TOOL_BUNDLE_CONTRACTS.items()
}
PAPER_RUNTIME_TOOL_BUNDLE_CONTRACTS["clean"] = {
    **PAPER_TOOL_BUNDLE_CONTRACTS["clean"],
    "identity_scope": "post_install_runtime_tree",
    "source_archive_identity": {
        "file_count": PAPER_TOOL_BUNDLE_CONTRACTS["clean"]["file_count"],
        "bytes": PAPER_TOOL_BUNDLE_CONTRACTS["clean"]["bytes"],
        "sha256": PAPER_TOOL_BUNDLE_CONTRACTS["clean"]["sha256"],
    },
    "file_count": 506,
    "bytes": 351120440,
    "sha256": (
        "a4ff2fbe213f86671559a7a3be2bbf53b8edbf89f6c65680e73b5c07375fc1d2"
    ),
}
