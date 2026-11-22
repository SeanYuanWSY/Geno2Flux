AlphaGEM
======
_AlphaGEM Enables Precise Genome-Scale Metabolic Modelling by Integrating Protein Structure Alignment with deep-learning-based Dark Metabolism Mining_
_doi: https://doi.org/10.1101/2025.07.21.665674_

Paper-reproduction branch
======

This fork repairs the public entry point and keeps paper-method and released-code
semantics separate. The default `--profile paper` encodes the paper-declared
method contract, with any necessary operationalizations labeled explicitly:

- PLMSearch `0.8`, pairwise blastp coverage `>=0.6` and identity `>=30%`;
- DIAMOND 2.1.10 reciprocal best hits, coverage `>=0.5` and identity `>=60%`;
- domain-specific structure thresholds from Supplementary Table 1;
- 2-D PCA ProtTrans filtering with one `mu=1.1`;
- mass **and charge** balance for reaction-pool additions;
- exactly the selected minimal or full medium during gap filling;
- explicit raw target-genome `tblastn` validation with bit score `>200`.

The paper names COBRApy 0.29.1 for gap filling but does not name a solver
backend. `AlphaGEM.yml` pins Gurobi 11.0.3 as an available dependency and
Python 3.10.15 as the author environment; neither proves that Gurobi was the
selected backend, and neither is a paper-declared pipeline requirement. The
paper states Python 3.10.15 only for its statistical-analysis section. Every formal
paper-profile invocation therefore requires an explicit `--solver` so the
observed backend is auditable, but a successful Gurobi or another MILP-capable
backend is recorded as an observed execution choice rather than proof of an
exact paper environment.

Every formal paper-profile invocation also requires a versioned
`--target-input-contract` JSON. The contract pins the exact FASTA
bytes, record count, normalized protein-ID digest, canonical sequence digest,
case name, and provenance. Formal runs reject `--solver auto`, changed target
inputs, duplicate normalized protein IDs, parameter overrides, dirty/uncommitted
source trees, and solvers that cannot enforce every declared gap-fill
tolerance. Override profiles remain available only as explicitly unverified
sensitivity/control runs.

The gap-fill stage additionally requires `--gapfill-genome`,
`--gapfill-coordinate-map`, `--gapfill-coordinate-segments`, and
`--gapfill-coordinate-contract`. The last JSON pins all three coordinate
artifacts plus the target FASTA to the official source URLs, audit command,
audit-script hash, resource-manifest hash, and 1-based segment-union semantics.
The first file is the raw target-organism genome FASTA stated in the Methods.
The envelope TSV
maps every target protein to exactly one genomic locus, with columns
`target_protein_id`, `genome_seqid`, `start`, `end`, `strand`,
`refseq_protein_id`, `refseq_cds_id`, and `mapping_evidence`. The segment TSV
contains one row per CDS segment, using `segment_index_genomic_order`,
`segment_start`, `segment_end`, `phase`, and `locus_tag` in addition to the
protein, RefSeq, seqid, and strand fields. TBLASTN hits are matched against the
CDS-segment union, never the intron-spanning envelope. These explicit contracts
restore the hidden coordinate-to-protein step missing from the released
implementation. An incomplete, inconsistent, or ambiguous mapping is rejected.
`--gapfill-cds` remains only as a labeled legacy proxy and is not accepted by
the paper profile.

`--profile author_artifact_legacy` is a sensitivity-control profile for the
released implementation: direct NCBI BLASTP additions, legacy long-sequence ESM
pooling, and the historical minimal-plus-full gap-fill union. It must not be
reported as the paper Methods result. `github_current` and `legacy_main` preserve
the two public-code default sets for further diagnosis.

Every formal run writes `working/<name>/alphagem_run_manifest.json`,
`gapfill_audit.json`, and `model_validation.json`. Resume is rejected when the
FASTA, profile, parameters, source commit, runtime-package versions, or
upstream artifact hashes differ. The manifest records the source commit,
code-file digest, Git cleanliness, Python runtime, and key dependency versions.
The manifest labels the result as a `paper_declared_method_contract`, never as
an exact paper-runtime reproduction, and keeps paper declarations, author
environment pins, and observed runtime/solver details separate.
Existing ESM embeddings are reused only when their provenance sidecar matches
the staged FASTA and model weights. BBH selection preserves the released
implementation: after coverage filtering, the highest percent identity is
selected, with equal-identity hits resolved by their stable row order in the
DIAMOND output rather than by bit score or lexical subject order.

The six precomputed Rhea-search embedding shards are a different evidence
class. Their bytes are pinned to the verified author data archive and are
checked against an allowlist before the released PLMSearch code unpickles
them. The paper and released resources do not disclose the upstream UniProt
snapshot, embedding-model identity, pooling semantics, or sharding rule, and
they do not include the shard FASTAs. The formal manifest therefore records
bitwise reuse of author-provided derived artifacts and explicitly does not
claim that their upstream generation was reproduced. Synthetic provenance
sidecars are neither required nor accepted as evidence for these shards.

For the C. albicans PLMSearch reproduction on the RTX PRO 6000 server:

    python AlphaGEM.py \
      --profile paper \
      --mode plmsearch \
      --name candida-paper-plm \
      --refname yeast \
      --fasta /absolute/path/to/candida.fasta \
      --target-input-contract /absolute/path/to/candida.input-contract.json \
      --gapfill-genome /absolute/path/to/candida.genomic.fna \
      --gapfill-coordinate-map /absolute/path/to/candida.coordinate-map.tsv \
      --gapfill-coordinate-segments /absolute/path/to/candida.coordinate-segments.tsv \
      --gapfill-coordinate-contract /absolute/path/to/candida.coordinate-contract.json \
      --solver gurobi \
      --deepec-python /absolute/path/to/deepectransformer/bin/python \
      --blast-threads 104 \
      --esm-model-path /absolute/path/to/esm1b_t33_650M_UR50S.pt \
      --esm-token-budget 98304

For `--mode structure_alignment`, the paper profile additionally requires a
fully downloaded local `Rostlab/prot_t5_xl_uniref50` snapshot:

    --prottrans-model-path /absolute/path/to/prot_t5_xl_uniref50

Every file in that snapshot is hashed into the homology-stage manifest, and
the loader runs with `local_files_only=True`.

Run the legacy control under a different job name so its artifacts cannot mix
with the paper run:

    python AlphaGEM.py \
      --profile author_artifact_legacy \
      --mode plmsearch \
      --name candida-author-legacy \
      --refname yeast \
      --fasta /absolute/path/to/candida.fasta \
      --solver gurobi \
      --blast-threads 104 \
      --esm-model-path /absolute/path/to/esm1b_t33_650M_UR50S.pt

Run the repair contract tests before launching:

    python -m unittest discover -s tests -v

 >## Platform: Ubuntu only (tested on Ubuntu 20.04/22.04). #

 
Installation:
======
## Install codes and datas: #
You can install AlphaGEM through git
   
    git clone https://github.com/hongzhonglu/AlphaGEMs.git AlphaGEM

You can download supplyment resources from figshare:

[Data Source Figshare](https://dx.doi.org/10.6084/m9.figshare.30630524)

## Prepare envrionment: #
It need to create two enviroment as DeepECtransformer uses the old version of pytorch.

* __1.create environment for AlphaGEM:__

      conda env create -f AlphaGEM.yml

* __2.create enviroment for deepectransformer:__
  
      conda env create -f deepecenvironment.yml
  
* __3.activate AlphaGEM environment:__
  
      conda activate AlphaGEM
  
* __4.install PyTorch/CUDA/cuDNN:__
 
   Please install a compatible combination of `PyTorch`, `CUDA`, and `cuDNN` for your system yourself. (Recommended `PyTorch` version: >= 2.4.0)
  
   Use the official selector to obtain the correct command for your OS and CUDA version: [PyTorch Local Installation](https://pytorch.org/get-started/locally/)
  
   >Tested platform for this project: `PyTorch 2.9.0` + `CUDA 13.0` on `Nvidia GeForce RTX 5080` and `RTX 5070ti`.

* __5.install ESM (provides `esm` package used by embeddings):__

      pip install fair-esm==2.0.0

* __6.prepare the source data:__
  
      python setup.py install

* __7.Gurobi license (if you don’t have one yet)__

  This is a released-code/author-environment recommendation, not a solver
  requirement stated by the paper. To reproduce the author environment with
  Gurobi, obtain a valid license; otherwise select another explicitly recorded
  MILP-capable optimizer. To retrieve a Gurobi license:

      grbgetkey <YOUR-LICENSE-KEY>
  
  For details, see the Gurobi licensing guide: https://www.gurobi.com/documentation/

Test:
======
   
   ## Create GEM for _C. albicans_ #
   
   Released-code smoke only (not the formal paper profile), structure mode:
   
      python AlphaGEM.py --profile github_current --name candida --mode structure_alignment --fasta ./working/candida/candida.fasta --refname yeast
      
   Released-code smoke only, PLMSearch mode:
   
      python AlphaGEM.py --profile github_current --name candida --mode plmsearch --fasta ./working/candida/candida.fasta --refname yeast
 
Usage:
======

 ## Parameters: #
  - **--mode**: Workflow mode. Default: `structure_alignment`. Choices: `structure_alignment`, `plmsearch`
  - **--refname**: Reference species. Choices: `ecoli`, `yeast`, `strco`, `human`, `synechocystis`. Required
  - **--name**: Job/species name used to create `working/<name>`. Required
  - **--fasta**: Target species FASTA file path. Required
  - **--target-input-contract**: JSON pinning the exact target FASTA identity and provenance. Required by every formal paper run
  - **--solver**: Explicit COBRA/optlang backend. `auto` is rejected by formal paper runs
  - **--gapfill-genome**: Raw target-organism genome FASTA used by the paper `tblastn` gap-fill validation. Required when the paper run includes `gapfill`
  - **--gapfill-coordinate-map**: Complete audited target-protein-to-genome envelope TSV. Required with `--gapfill-genome`
  - **--gapfill-coordinate-segments**: Complete one-row-per-CDS-segment TSV used for exon-union overlap. Required with `--gapfill-genome`
  - **--gapfill-coordinate-contract**: JSON pinning coordinate artifacts, exact semantics, official sources, and audit provenance. Required with `--gapfill-genome`
  - **--gapfill-cds**: Legacy CDS proxy input; never treated as the paper Methods path
  - **--prottrans-model-path**: Local `Rostlab/prot_t5_xl_uniref50` snapshot required by the paper structure-homology stage
  - **--list**: Path to a list file of structures and gene names. Default: empty
  - **--structure**: Directory containing structure files. Default: empty
  - **--cleanuse**: Whether CLEAN has been used. To enable, pass `--cleanuse True`; to disable, omit the flag
  - **--USalignfilter**: USalign filter threshold. Default: 0.5
  - **--TMscore**: Structure alignment filter threshold(foldseek). Default: 0.7
  - **--upTMscore**: Safe TMscore(foldseek). Default: 0.9
  - **--TMscoretrans**: Transporter TMscore filter (foldseek). Default: 0.7
  - **--coverage**: Coverage filter threshold(foldseek). Default: 0.8
  - **--upcoverage**: Safe coverage(foldseek). Default: 0.9
  - **--coveragetrans**: Transporter coverage filter(foldseek). Default: 0.8
  - **--pLDDT**: pLDDT filter threshold. Default: 70
  - **--esp**: Clustering parameter. Default: 1
  - **--grothmedium**: Growth medium. Default: `min`; choices: `min`, `full`

 ## Tutorials: #

 ## Examples: #
  - **Structure alignment mode (default)**
  
    > **Pre-requisites for Structure Alignment:**
    > If you perform structure alignment, please provide:
    > 1. A directory containing structure data via the `--structure` parameter.
    > 2. A genome-wide table file via the `--list` parameter. 
    >
    > **List File Format Requirements:**
    > *   The table must contain the following columns: **Entry**, **Entry Name**, **Structure**.
    > *   The content of the **Entry** column must strictly match the gene IDs provided in your `--fasta` file.

        python AlphaGEM.py \
          --mode "structure alignment" \
          --name my_species \
          --refname yeast \
          --fasta ./my_species.fasta \
          --structure ./structure_dir \
          --list ./my_species.xlsx \
          --cleanuse True \
          --TMscore 0.7 --coverage 0.8 --pLDDT 70

  - **PLMSearch mode**

        python AlphaGEM.py \
          --mode plmsearch \
          --name my_species \
          --refname ecoli \
          --fasta ./my_species.fasta
