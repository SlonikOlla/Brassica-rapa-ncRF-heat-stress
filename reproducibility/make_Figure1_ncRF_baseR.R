#!/usr/bin/env Rscript

SRC <- "results/ncRF_17_26/figures/source_data"
OUT <- "results/ncRF_17_26/figures/final"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

rd <- function(x) {
    read.delim(
        file.path(SRC,x),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )
}

A <- rd("Figure_ncRF_panelA_stage_remodeling.tsv")
B <- rd("Figure_ncRF_panelB_Endosperm_classes.tsv")
C <- rd("Figure_ncRF_panelC_LSU3prime_CANONICAL.tsv")
D <- rd("Figure_ncRF_panelD_Endosperm_effects.tsv")
E <- rd("Figure_ncRF_panelE_snoR21_targets.tsv")
F <- rd("Figure_ncRF_panelF_fraction_negative_nulls.tsv")
FC <- rd("Figure_ncRF_panelF_combined_effect_null.tsv")

# ------------------------------------------------------------
# Frozen-result checks
# ------------------------------------------------------------

stopifnot(
    nrow(A)==8,
    sum(A$strong_total)==39,
    A$strong_total[A$stage=="Endosperm"]==33,
    nrow(C)==10,
    nrow(D)==33,
    sum(D$direction=="UP")==31,
    sum(D$direction=="DOWN")==2,
    sum(D$normalization_robust=="YES")==26,
    nrow(E)==25,
    sum(E$log2FC<0)==19,
    sum(E$response=="strong_DOWN")==9,
    sum(E$response=="strong_UP")==2,
    sum(F$null_type=="expression_matched")==10000,
    sum(F$null_type=="dinucleotide_sequence")==1000,
    nrow(FC)==1000
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

panel <- function(letter) {
    mtext(
        letter,
        side=3,
        line=0.35,
        adj=-0.10,
        font=2,
        cex=1.35
    )
}

endpoint_label <- function(x) {
    x <- gsub("3prime_CCA","3'-CCA",x,fixed=TRUE)
    x <- gsub("5prime","5'",x,fixed=TRUE)
    x <- gsub("3prime","3'",x,fixed=TRUE)
    x
}

# Colors chosen consistently by biological direction/class.
col_up <- "#B2182B"
col_down <- "#2166AC"
col_neutral <- "#666666"
col_robust <- "#222222"
col_sensitive <- "#AAAAAA"

# ============================================================
# Drawing function
# ============================================================

draw_figure <- function() {

    layout(
        matrix(1:6,nrow=3,byrow=TRUE),
        widths=c(1,1),
        heights=c(0.82,1.42,1.12)
    )

    par(
        oma=c(1.2,1.0,2.5,0.5),
        mgp=c(2.0,0.55,0),
        tcl=-0.25,
        las=1
    )

    # ========================================================
    # A — stage specificity
    # ========================================================

    par(mar=c(6.0,4.2,2.5,1.0))

    stage_labels <- c(
        "Leaves",
        "Inflorescence\nmeristem",
        "Pollen",
        "Unfertilized\novules",
        "Fertilized\novules",
        "Embryo",
        "Endosperm",
        "Progeny"
    )

    ymax <- max(A$strong_UP)+5
    ymin <- -max(A$strong_DOWN)-3

    bp <- barplot(
        rbind(A$strong_UP,-A$strong_DOWN),
        beside=FALSE,
        col=c(col_up,col_down),
        border=NA,
        names.arg=rep("",nrow(A)),
        ylim=c(ymin,ymax),
        ylab="Strong endpoint families",
        main="Stage-specific ncRF remodeling",
        cex.main=0.95
    )

    abline(h=0,col="black",lwd=0.8)

    axis(
        1,
        at=bp,
        labels=stage_labels,
        las=2,
        tick=FALSE,
        cex.axis=0.66
    )

    for (i in seq_len(nrow(A))) {
        if (A$strong_UP[i] > 0) {
            text(
                bp[i],
                A$strong_UP[i]+0.7,
                A$strong_UP[i],
                cex=0.72
            )
        }

        if (A$strong_DOWN[i] > 0) {
            text(
                bp[i],
                -A$strong_DOWN[i]-0.6,
                A$strong_DOWN[i],
                cex=0.72
            )
        }
    }

    legend(
        "topleft",
        legend=c("UP","DOWN"),
        fill=c(col_up,col_down),
        border=NA,
        bty="n",
        cex=0.72,
        horiz=TRUE
    )

    panel("A")

    # ========================================================
    # B — Endosperm class composition
    # ========================================================

    par(mar=c(4.2,5.0,2.5,1.0))

    orderB <- match(
        c("tRNA","snoRNA","rRNA","snRNA"),
        B$ncRNA_class
    )

    BB <- B[orderB,]

    lim <- max(BB$strong_UP)+3

    yy <- barplot(
        rbind(BB$strong_UP,-BB$strong_DOWN),
        horiz=TRUE,
        beside=FALSE,
        col=c(col_up,col_down),
        border=NA,
        names.arg=BB$ncRNA_class,
        xlim=c(-3,lim),
        xlab="Strong endpoint families",
        main="Endosperm: 31 UP / 2 DOWN",
        cex.main=0.95
    )

    abline(v=0,lwd=0.8)

    for (i in seq_len(nrow(BB))) {
        if (BB$strong_UP[i]>0) {
            text(
                BB$strong_UP[i]+0.45,
                yy[i],
                BB$strong_UP[i],
                cex=0.75
            )
        }

        if (BB$strong_DOWN[i]>0) {
            text(
                -BB$strong_DOWN[i]-0.35,
                yy[i],
                BB$strong_DOWN[i],
                cex=0.75
            )
        }
    }

    text(
        0.98,0.08,
        "26/33 robust to alternative normalization",
        adj=c(1,0),
        xpd=NA,
        cex=0.68
    )

    panel("B")

    # ========================================================
    # C — LSU 3' nested terminal series
    # ========================================================

    par(mar=c(4.3,11.2,2.7,1.0))

    C <- C[order(C$length),]

    mat <- as.matrix(
        C[,c("CT1","CT2","TR1","TR2")]
    )

    z <- log10(mat+1)

    # image() expects x by y; transpose and reverse rows.
    pal <- hcl.colors(100,"YlOrRd",rev=FALSE)

    image(
        x=1:4,
        y=1:nrow(C),
        z=t(z[nrow(C):1,]),
        col=pal,
        axes=FALSE,
        xlab="Endosperm library",
        ylab="",
        main="Nested LSU 3'-terminal ncRFs",
        cex.main=0.95
    )

    axis(
        1,
        at=1:4,
        labels=c("CT1","CT2","TR1","TR2")
    )

    seq_labels <- paste0(
        C$length," nt  ",C$sequence
    )

    axis(
        2,
        at=1:nrow(C),
        labels=rev(seq_labels),
        las=1,
        tick=FALSE,
        cex.axis=0.55
    )

    box()

    mtext(
        "Canonical series: 99.7% CT; 98.8% TR LSU 3' reads",
        side=1,
        line=2.4,
        cex=0.67
    )

    mtext(
        "Family-level TMM: log2FC = +4.862; FDR = 2.76e-5",
        side=1,
        line=3.2,
        cex=0.67
    )

    panel("C")

    # ========================================================
    # D — all 33 Endosperm effects
    # ========================================================

    par(mar=c(4.2,8.2,2.7,1.0))

    DD <- D[order(D$log2FC),]

    labsD <- paste(
        DD$family,
        endpoint_label(DD$endpoint)
    )

    yy <- seq_len(nrow(DD))

    plot(
        DD$log2FC,
        yy,
        type="n",
        yaxt="n",
        ylab="",
        xlab="log2 fold change (TMM)",
        main="33 strong Endosperm endpoint families",
        cex.main=0.95,
        xlim=range(c(DD$log2FC,0))+c(-0.4,0.4)
    )

    abline(v=0,lty=2,col="grey50")

    segments(
        0,yy,
        DD$log2FC,yy,
        col="grey75",
        lwd=0.8
    )

    point_col <- ifelse(
        DD$direction=="UP",
        col_up,
        col_down
    )

    pchD <- ifelse(
        DD$normalization_robust=="YES",
        19,
        21
    )

    points(
        DD$log2FC,
        yy,
        pch=pchD,
        col=point_col,
        bg="white",
        cex=0.82
    )

    axis(
        2,
        at=yy,
        labels=labsD,
        las=1,
        tick=FALSE,
        cex.axis=0.50
    )

    legend(
        "bottomright",
        legend=c(
            "Robust to alternative normalization",
            "Primary TMM only"
        ),
        pch=c(19,21),
        pt.bg=c(NA,"white"),
        bty="n",
        cex=0.60
    )

    panel("D")

    # ========================================================
    # E — snoR21 candidate targets
    # ========================================================

    par(mar=c(4.2,7.7,2.7,1.0))

    EE <- E[order(E$log2FC),]
    yy <- seq_len(nrow(EE))

    plot(
        EE$log2FC,
        yy,
        type="n",
        yaxt="n",
        ylab="",
        xlab="mRNA log2 fold change",
        main="snoR21 3' ncRF candidate-target responses",
        cex.main=0.92,
        xlim=range(c(EE$log2FC,0))+c(-0.25,0.25)
    )

    abline(v=0,lty=2,col="grey50")

    segments(
        0,yy,
        EE$log2FC,yy,
        col="grey75"
    )

    ecol <- ifelse(
        EE$log2FC<0,
        col_down,
        col_up
    )

    points(
        EE$log2FC,
        yy,
        pch=19,
        col=ecol,
        cex=0.78
    )

    gene_labels <- sub(
        "\\.1_BraROA$","",
        EE$gene_id
    )

    axis(
        2,
        at=yy,
        labels=gene_labels,
        las=1,
        tick=FALSE,
        cex.axis=0.53
    )

    legend(
        "bottomright",
        legend=c(
            "19/25 negative (76%)",
            "9 strong DOWN",
            "2 strong UP"
        ),
        bty="n",
        cex=0.65
    )

    panel("E")

    # ========================================================
    # F — null testing
    # ========================================================

    par(mar=c(5.1,4.4,2.7,1.0))

    expr <- F$fraction_negative[
        F$null_type=="expression_matched"
    ]

    seqnull <- F$fraction_negative[
        F$null_type=="dinucleotide_sequence"
    ]

    comb <- FC$sequence_minus_matched

    # Three standardized vertical distributions.
    # First two are fraction negative.
    # Third is excess fraction negative relative to matched genes.

    plot(
        NA,
        xlim=c(0.45,3.55),
        ylim=c(
            min(c(expr,seqnull,comb,0.250772))-0.05,
            max(c(expr,seqnull,0.76))+0.08
        ),
        xaxt="n",
        xlab="",
        ylab="Fraction negative / excess fraction negative",
        main="snoR21 sequence-specific inverse association",
        cex.main=0.90
    )

    boxplot(
        expr,
        at=1,
        add=TRUE,
        axes=FALSE,
        outline=FALSE,
        width=0.48
    )

    boxplot(
        seqnull,
        at=2,
        add=TRUE,
        axes=FALSE,
        outline=FALSE,
        width=0.48
    )

    boxplot(
        comb,
        at=3,
        add=TRUE,
        axes=FALSE,
        outline=FALSE,
        width=0.48
    )

    points(
        c(1,2),
        c(0.76,0.76),
        pch=18,
        cex=1.35,
        col=col_up
    )

    points(
        3,
        0.250772,
        pch=18,
        cex=1.35,
        col=col_up
    )

    abline(h=0,lty=3,col="grey70")

    axis(
        1,
        at=1:3,
        labels=c(
            "Expression-\nmatched",
            "Dinucleotide\nsequence",
            "Sequence minus\nmatched"
        ),
        tick=FALSE,
        cex.axis=0.72
    )

    text(
        1,0.815,
        "P = 0.0079",
        cex=0.68
    )

    text(
        2,0.815,
        "P = 0.0060",
        cex=0.68
    )

    text(
        3,0.32,
        "P = 0.038",
        cex=0.68
    )

    legend(
        "topleft",
        legend="Authentic snoR21",
        pch=18,
        col=col_up,
        bty="n",
        cex=0.68
    )

    mtext(
        "Combined statistic = sequence-selected fraction minus expression-matched fraction",
        side=1,
        line=3.7,
        cex=0.60
    )

    panel("F")

    mtext(
        "Endosperm-specific remodeling of ncRNA-derived terminal fragments",
        outer=TRUE,
        side=3,
        line=0.7,
        font=2,
        cex=1.25
    )
}

# ------------------------------------------------------------
# Vector PDF
# ------------------------------------------------------------

pdf(
    file.path(
        OUT,
        "Figure1_ncRF_Endosperm_remodeling.pdf"
    ),
    width=12,
    height=14,
    useDingbats=FALSE
)

draw_figure()
dev.off()

# ------------------------------------------------------------
# 600-dpi PNG
# ------------------------------------------------------------

png(
    file.path(
        OUT,
        "Figure1_ncRF_Endosperm_remodeling.png"
    ),
    width=12,
    height=14,
    units="in",
    res=600,
    type="cairo"
)

draw_figure()
dev.off()

cat(
    "WROTE",
    file.path(
        OUT,
        "Figure1_ncRF_Endosperm_remodeling.pdf"
    ),
    "\n"
)

cat(
    "WROTE",
    file.path(
        OUT,
        "Figure1_ncRF_Endosperm_remodeling.png"
    ),
    "\n"
)
