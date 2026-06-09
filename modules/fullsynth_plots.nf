process FULLSYNTH_PLOTS {

    publishDir "$params.outdir"

    input:
    path(deseq_synths,  stageAs: 'deseq/outframe_?.csv')
    path(edger_synths,  stageAs: 'edger/outframe_?.csv')
    path(wilcox_synths, stageAs: 'wilcox/outframe_?.csv')
    path(svc_synths,    stageAs: 'svc/outframe_?.csv')

    output:
    path "synth_plot.pdf"

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("reshape2", quietly = TRUE)
    library("viridis")

    read_fullsynths = function(dir, label){
        lapply(list.files(dir, full.names = TRUE), read.csv) %>%
            bind_rows() %>%
            mutate(method = label) %>%
            melt(id.vars = c("method", "samplenum", "depth"))
    }

    deseq.synth.input  = read_fullsynths("deseq",  "DESeq2")
    edger.synth.input  = read_fullsynths("edger",  "edgeR")
    wilcox.synth.input = read_fullsynths("wilcox", "Wilcoxon")

    if("$params.mlstep" == "true"){
        svc.synth.input = read_fullsynths("svc", "SVC")
        gg.synth.input = rbind(deseq.synth.input,
                               edger.synth.input,
                               wilcox.synth.input,
                               svc.synth.input)
    } else {
        gg.synth.input = rbind(deseq.synth.input,
                               edger.synth.input,
                               wilcox.synth.input)
    }

    gg.synth.input = gg.synth.input %>%
        mutate(samplenum = as.factor(samplenum),
               depth     = as.factor(depth)) %>%
        group_by(method, samplenum, depth, variable) %>%
        summarize(mean = mean(value)) %>%
        mutate(mean = signif(mean, 3))

    # Take FDP as inverse so higher values are better (matches power and normAUC)
    gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")] =
        1 - gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")]
    gg.synth.input\$variable = fct_recode(gg.synth.input\$variable, "1-FDP" = "FDP")

    gg.synth = ggplot(gg.synth.input, aes(x = depth, y = samplenum)) +
        geom_tile(aes(fill = mean), color = "gray20", size = .75, width = 1, height = 1) +
        geom_label(aes(label = mean), color = "white",
                   lineheight = 0.75, size = 3.5, alpha = 0.7, fill = "gray10") +
        scale_x_discrete("Sequencing depth",  expand = c(0, 0)) +
        scale_y_discrete("Replicates per treatment", expand = c(0, 0)) +
        scale_fill_viridis(option = "plasma", limits = c(0, 1), expand = c(0, 0)) +
        facet_grid(variable ~ method) +
        theme_bw() +
        theme(panel.grid      = element_blank(),
              legend.title    = element_text(face = "bold"),
              axis.text       = element_text(face = "bold", size = 11, color = "grey40"),
              axis.title      = element_text(face = "bold", size = 12))

    ggsave(gg.synth,
           filename = "synth_plot.pdf",
           device = "pdf", bg = "transparent",
           width  = 10 * length(unique(gg.synth.input\$method)),
           height = 30, units = "cm")
    """
}
