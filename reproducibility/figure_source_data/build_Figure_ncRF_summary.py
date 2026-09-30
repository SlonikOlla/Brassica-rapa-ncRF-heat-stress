from pathlib import Path
import csv

OUT = Path(
    "results/ncRF_17_26/figures/source_data/"
    "Figure_ncRF_central_summary.tsv"
)

rows = [
    # Panel A: stage-wide endpoint-family remodeling
    ["A","Leaves","tested",107],
    ["A","Leaves","strong_UP",0],
    ["A","Leaves","strong_DOWN",3],

    ["A","Inflorescence_meristem","tested",58],
    ["A","Inflorescence_meristem","strong_UP",0],
    ["A","Inflorescence_meristem","strong_DOWN",0],

    ["A","Pollen","tested",42],
    ["A","Pollen","strong_UP",0],
    ["A","Pollen","strong_DOWN",1],

    ["A","Unfertilized_ovules","tested",38],
    ["A","Unfertilized_ovules","strong_UP",0],
    ["A","Unfertilized_ovules","strong_DOWN",0],

    ["A","Fertilized_ovules","tested",52],
    ["A","Fertilized_ovules","strong_UP",0],
    ["A","Fertilized_ovules","strong_DOWN",0],

    ["A","Embryo","tested",69],
    ["A","Embryo","strong_UP",1],
    ["A","Embryo","strong_DOWN",1],

    ["A","Endosperm","tested",70],
    ["A","Endosperm","strong_UP",31],
    ["A","Endosperm","strong_DOWN",2],

    ["A","Progeny","tested",73],
    ["A","Progeny","strong_UP",0],
    ["A","Progeny","strong_DOWN",0],

    # Panel B: Endosperm strong ncRFs by class
    ["B","rRNA","strong_UP",3],
    ["B","rRNA","strong_DOWN",0],
    ["B","snoRNA","strong_UP",11],
    ["B","snoRNA","strong_DOWN",0],
    ["B","snRNA","strong_UP",1],
    ["B","snRNA","strong_DOWN",1],
    ["B","tRNA","strong_UP",16],
    ["B","tRNA","strong_DOWN",1],

    # Panel C: LSU 3'-terminal family
    ["C","LSU_3prime","CT_CPM",151],
    ["C","LSU_3prime","TR_CPM",5043],
    ["C","LSU_3prime","fold_change",33.3],

    # Panel D: major Endosperm endpoint-family effects
    ["D","LSU_3prime","log2FC",4.862],
    ["D","SSU_5prime","log2FC",3.231],
    ["D","U4_3prime","log2FC",5.485],
    ["D","U2_3prime","log2FC",-3.358],
    ["D","snoU61_5prime","log2FC",3.621],
    ["D","snoR98_3prime","log2FC",4.705],

    # Panel E: snoR21 high-complementarity RNA-seq targets
    ["E","snoR21","targets",25],
    ["E","snoR21","negative",19],
    ["E","snoR21","strong_DOWN",9],
    ["E","snoR21","strong_UP",2],
    ["E","snoR21","fraction_negative",0.76],
    ["E","snoR21","median_log2FC",-0.6245615],

    # Panel F: null-model tests
    ["F","expression_matched","fraction_negative_P",0.0078992],
    ["F","expression_matched","median_logFC_P",0.0159984],
    ["F","dinucleotide_sequence","fraction_negative_P",0.0059940],
    ["F","dinucleotide_sequence","median_logFC_P",0.0079920],
    ["F","combined","fraction_negative_P",0.0379620],
    ["F","combined","median_logFC_P",0.0529471],

    # Panel G: collapsed GO enrichment
    ["G","ABC_xenobiotic_transporter","FDR",0.003316],
    ["G","xenobiotic_transport","FDR",0.003316],
    ["G","ATP_hydrolysis_activity","FDR",0.006770],
    ["G","ABC_transporter_activity","FDR",0.009722],
    ["G","ATP_binding","FDR",0.009722],
]

with open(OUT,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["panel","feature","metric","value"])
    w.writerows(rows)

print("ROWS",len(rows))
print("WROTE",OUT)

print()
print("===== CENTRAL FIGURE CONTENT =====")
for r in rows:
    print(*r,sep="\t")
