from pathlib import Path
import csv
import statistics

RAND = Path("results/ncRF_17_26/targets/randomization")
OUTDIR = Path("results/ncRF_17_26/figures/source_data")

EXPR = RAND / "snoR21_authentic_expression_matched_10000.tsv"
SEQ  = RAND / "snoR21_dinuc_null_1000_target_summary.tsv"
COMB = RAND / "snoR21_sequence_null_genelevel_combined.tsv"

OUT = OUTDIR / "Figure_ncRF_panelF_fraction_negative_nulls.tsv"
OUTC = OUTDIR / "Figure_ncRF_panelF_combined_effect_null.tsv"

AUTH_FRAC = 19/25

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f,delimiter="\t"))

expr=read_tsv(EXPR)
seq=read_tsv(SEQ)
comb=read_tsv(COMB)

if len(expr) != 10000:
    raise RuntimeError(f"Expected 10000 expression nulls; found {len(expr)}")

if len(seq) != 1000:
    raise RuntimeError(f"Expected 1000 sequence nulls; found {len(seq)}")

if len(comb) != 1000:
    raise RuntimeError(f"Expected 1000 combined nulls; found {len(comb)}")

# ------------------------------------------------------------
# Expression-matched null
# ------------------------------------------------------------

expr_frac=[
    float(r["fraction_negative"])
    for r in expr
]

expr_mean=statistics.mean(expr_frac)

expr_extreme=sum(x >= AUTH_FRAC for x in expr_frac)
expr_p=(expr_extreme+1)/(len(expr_frac)+1)

# ------------------------------------------------------------
# Dinucleotide-preserving sequence null
# ------------------------------------------------------------

seq_frac=[
    float(r["fraction_negative"])
    for r in seq
]

seq_mean=statistics.mean(seq_frac)

seq_extreme=sum(x >= AUTH_FRAC for x in seq_frac)
seq_p=(seq_extreme+1)/(len(seq_frac)+1)

# ------------------------------------------------------------
# Combined sequence + expression adjustment
# ------------------------------------------------------------

combined_effect=[]

for r in comb:
    seqf=float(r["seq_fraction_negative"])
    matf=float(r["matched_fraction_negative"])

    combined_effect.append(
        (r["null_id"],seqf-matf)
    )

auth_effect=AUTH_FRAC-expr_mean

combined_extreme=sum(
    effect >= auth_effect
    for _,effect in combined_effect
)

combined_p=(
    (combined_extreme+1) /
    (len(combined_effect)+1)
)

# ------------------------------------------------------------
# Long-form source data for distributions
# ------------------------------------------------------------

with open(OUT,"w",newline="") as f:
    fields=["null_type","replicate","fraction_negative"]
    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
    w.writeheader()

    for i,x in enumerate(expr_frac,1):
        w.writerow({
            "null_type":"expression_matched",
            "replicate":i,
            "fraction_negative":x
        })

    for i,x in enumerate(seq_frac,1):
        w.writerow({
            "null_type":"dinucleotide_sequence",
            "replicate":i,
            "fraction_negative":x
        })

with open(OUTC,"w",newline="") as f:
    fields=[
        "null_id",
        "sequence_fraction_negative",
        "matched_fraction_negative",
        "sequence_minus_matched"
    ]

    w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
    w.writeheader()

    for r,(null_id,effect) in zip(comb,combined_effect):
        w.writerow({
            "null_id":null_id,
            "sequence_fraction_negative":
                r["seq_fraction_negative"],
            "matched_fraction_negative":
                r["matched_fraction_negative"],
            "sequence_minus_matched":effect
        })

# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

print("===== AUTHENTIC snoR21 =====")
print("N_TARGETS",25)
print("N_NEGATIVE",19)
print("FRACTION_NEGATIVE",f"{AUTH_FRAC:.6f}")

print()
print("===== EXPRESSION-MATCHED NULL =====")
print("N",len(expr_frac))
print("MEAN",f"{expr_mean:.9f}")
print("MEDIAN",f"{statistics.median(expr_frac):.9f}")
print("EXTREME",expr_extreme)
print("EMPIRICAL_P",f"{expr_p:.9f}")

print()
print("===== DINUCLEOTIDE SEQUENCE NULL =====")
print("N",len(seq_frac))
print("MEAN",f"{seq_mean:.9f}")
print("MEDIAN",f"{statistics.median(seq_frac):.9f}")
print("EXTREME",seq_extreme)
print("EMPIRICAL_P",f"{seq_p:.9f}")

effects=[x for _,x in combined_effect]

print()
print("===== COMBINED SEQUENCE + EXPRESSION =====")
print("AUTHENTIC_EFFECT",f"{auth_effect:.9f}")
print("NULL_EFFECT_MEAN",f"{statistics.mean(effects):.9f}")
print("NULL_EFFECT_MEDIAN",f"{statistics.median(effects):.9f}")
print("EXTREME",combined_extreme)
print("EMPIRICAL_P",f"{combined_p:.9f}")

# Verify frozen results.
expected=[
    (expr_p,0.00789921,"expression"),
    (seq_p,0.00599401,"sequence"),
    (combined_p,0.03796204,"combined")
]

for observed,expected_p,name in expected:
    if abs(observed-expected_p) > 1e-6:
        raise RuntimeError(
            f"{name} P mismatch: {observed} vs {expected_p}"
        )

print()
print("VERIFICATION_OK")
print("WROTE",OUT)
print("WROTE",OUTC)
