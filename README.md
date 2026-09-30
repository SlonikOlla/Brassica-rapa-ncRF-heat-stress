# Brassica rapa ncRNA-fragment heat-stress reanalysis

Reproducibility package for the reanalysis of authentic 17–26-nt structural noncoding-RNA fragments in the *Brassica rapa* R-o-18 heat-stress developmental series.

## Scope

This release contains analysis scripts, derived tables, quality-control outputs, endpoint-family differential-analysis results, normalization sensitivity analyses, exact-fragment query sets, ncRF–mRNA complementarity results, snoR21 randomization controls, corrected GO-enrichment outputs, and supplementary tables.

**The manuscript and all rendered manuscript/supplementary figures are intentionally excluded.** Source-data tables used to generate quantitative panels are retained under `reproducibility/figure_source_data/` because they are analysis outputs rather than figure artwork.

## Public input data

- Small RNA: GEO GSE58896 / SRA SRP043646 / BioProject PRJNA253869
- Matched RNA-seq: GEO GSE58895
- SuperSeries: GSE58897

The primary small-RNA reanalysis uses the original GEO processed per-sample sequence-count files preserving adapter-resolved 17–26-nt inserts. Fixed 29-nt SRA FASTQ sequences were not trimmed to manufacture shorter ncRFs.

## Main directories

- `results/ncRF_17_26/qc/` — library complexity and replicate concordance
- `results/ncRF_17_26/trna/` — tRNA exact-read and endpoint assignments
- `results/ncRF_17_26/rrna/` — rRNA assignments and LSU 3′ nested-fragment analyses
- `results/ncRF_17_26/sno_snRNA/` — snoRNA/snRNA endpoint analyses
- `results/ncRF_17_26/unified/` — unified endpoint-family matrices and stage-specific differential results
- `results/ncRF_17_26/normalization/` — normalization sensitivity analyses
- `results/ncRF_17_26/targets/` — complementarity, RNA-seq integration, snoR21 controls and GO analysis
- `reproducibility/figure_source_data/` — frozen quantitative source tables used for manuscript panels; no rendered figures
- `Brassica_rapa_Supplementary_Tables_S1-S8.xlsx` — supplementary tables
- `FILE_MANIFEST.txt` — complete release inventory

## Key software

The analysis used tRNAscan-SE 2.0.12, Rfam 14, Infernal 1.1.5, edgeR/TMM, NCBI BLAST+ 2.17.0, DIAMOND 2.1.22, R and custom Python scripts. Exact package/runtime versions are preserved where available in scripts and analysis outputs.

## Statistical scope

Small-RNA groups contain two libraries per condition. Differential analysis is therefore interpreted conservatively. ncRF abundance changes are relative/compositional because external spike-ins were not available. The snoR21 ncRF–mRNA analysis is computational and associative; it does not establish direct binding, AGO loading, cleavage, or causality.

## Citation

Please cite the associated article once published. Until then, cite this GitHub release by version/tag and release URL.

## License

No license is asserted in this release package. Add the project license before public reuse terms are declared.

## Files omitted from the GitHub-sized package

Three bulky/regenerable classes are intentionally omitted from this GitHub package: the full rRNA observed-mapping table, the global all-observed-sequences table, and the local BLAST nucleotide database index files. These are intermediate products rather than manuscript figures or primary summary results. They can be regenerated from the public input data and the included analysis workflow. This keeps every repository file below GitHub's 100 MB per-file limit.
