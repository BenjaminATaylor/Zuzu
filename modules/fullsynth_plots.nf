process FULLSYNTH_PLOTS {

    publishDir "$params.outdir"

    input:
    path(deseq_synths,  stageAs: 'deseq/outframe_?.csv')
    path(edger_synths,  stageAs: 'edger/outframe_?.csv')
    path(wilcox_synths, stageAs: 'wilcox/outframe_?.csv')
    path(svc_synths,    stageAs: 'svc/outframe_?.csv')

    output:
    path "synth_plot.pdf"
    path "synth_results.csv"
    path "synth_summary.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("reshape2", quietly = TRUE)
    library("viridis")

    fullsynth.columns = c(
        "method", "power", "FDP", "youden_index",
        "true_positives", "false_positives", "true_negatives", "false_negatives",
        "n_tested", "n_called", "n_truth", "n_background", "metric_status",
        "replicate", "refnum", "altnum", "totalnum", "depth"
    )

    # Read each method only when it follows the common synthetic result schema
    read_fullsynths = function(dir, label){
        files = list.files(dir, full.names = TRUE)
        if(length(files) == 0){
            stop(sprintf("No %s synthetic result files were supplied", label))
        }

        frames = lapply(files, function(file){
            frame = read.csv(file, check.names = FALSE)
            missing.columns = setdiff(fullsynth.columns, names(frame))
            unexpected.columns = setdiff(names(frame), fullsynth.columns)
            if(length(missing.columns) > 0 || length(unexpected.columns) > 0){
                stop(sprintf(
                    paste0(
                        "%s synthetic result %s does not match the common schema; ",
                        "missing: %s; unexpected: %s"
                    ),
                    label,
                    basename(file),
                    if(length(missing.columns) == 0) "none" else paste(missing.columns, collapse = ", "),
                    if(length(unexpected.columns) == 0) "none" else paste(unexpected.columns, collapse = ", ")
                ))
            }
            frame = frame[,fullsynth.columns,drop=FALSE]
            if(any(frame\$refnum < 1 | frame\$altnum < 1 |
                   frame\$totalnum != frame\$refnum + frame\$altnum)){
                stop(sprintf("%s synthetic result %s has invalid group sizes", label, basename(file)))
            }
            frame
        })

        bind_rows(frames) %>%
            mutate(method = label)
    }

    # Combine per-method results into one replicate-level table
    deseq.synth.input  = read_fullsynths("deseq",  "DESeq2")
    edger.synth.input  = read_fullsynths("edger",  "edgeR")
    wilcox.synth.input = read_fullsynths("wilcox", "Wilcoxon")

    if("$params.mlstep" == "true"){
        svc.synth.input = read_fullsynths("svc", "SVC-RFECV candidates")
        fullsynth.results = bind_rows(deseq.synth.input,
                                      edger.synth.input,
                                      wilcox.synth.input,
                                      svc.synth.input)
    } else {
        fullsynth.results = bind_rows(deseq.synth.input,
                                      edger.synth.input,
                                      wilcox.synth.input)
    }

    fullsynth.results = fullsynth.results %>%
        arrange(method, replicate, totalnum, refnum, altnum, depth)
    write.csv(fullsynth.results, "synth_results.csv", row.names = FALSE)

    # Reshape and summarise benchmark metrics across replicates
    fullsynth.long = melt(
        fullsynth.results,
        id.vars = c("method", "replicate", "refnum", "altnum", "totalnum", "depth", "metric_status"),
        measure.vars = c("power", "FDP", "youden_index")
    )

    fullsynth.summary = fullsynth.long %>%
        group_by(method, refnum, altnum, totalnum, depth, variable) %>%
        summarize(
            mean = if(any(is.finite(value))){
                mean(value[is.finite(value)])
            } else {
                NA_real_
            },
            n_valid = sum(is.finite(value)),
            n_total = n(),
            n_no_calls = sum(metric_status == "no_calls"),
            .groups = "drop"
        )
    write.csv(fullsynth.summary, "synth_summary.csv", row.names = FALSE)

    # Label and order each requested reference/alternative sample-size pair
    sample.size.levels = fullsynth.summary %>%
        distinct(refnum, altnum, totalnum) %>%
        arrange(totalnum, refnum, altnum) %>%
        transmute(label = paste0(refnum, " ref + ", altnum, " alt")) %>%
        pull(label)

    gg.synth.input = fullsynth.summary %>%
        mutate(mean      = signif(mean, 3),
               sample_size = factor(
                   paste0(refnum, " ref + ", altnum, " alt"),
                   levels = sample.size.levels
               ),
               depth     = as.factor(depth))

    # Take FDP as inverse so higher values are better (matching power and
    # Youden's index).
    gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")] =
        1 - gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")]
    gg.synth.input\$variable = fct_recode(
        gg.synth.input\$variable,
        "1-FDP" = "FDP",
        "Youden's index" = "youden_index"
    )
    gg.synth.input = gg.synth.input %>%
        mutate(plot_label = ifelse(
            is.na(mean),
            paste0("NA\n(n=0/", n_total, ")"),
            paste0(mean, "\n(n=", n_valid, "/", n_total, ")")
        ))

    gg.synth = ggplot(gg.synth.input, aes(x = depth, y = sample_size)) +
        geom_tile(aes(fill = mean), color = "gray20", size = .75, width = 1, height = 1) +
        geom_label(aes(label = plot_label), color = "white",
                   lineheight = 0.75, size = 3.5, alpha = 0.7, fill = "gray10") +
        scale_x_discrete("Sequencing depth",  expand = c(0, 0)) +
        scale_y_discrete("Sample size (reference + alternative)", expand = c(0, 0)) +
        scale_fill_viridis(option = "plasma", limits = c(-1, 1), expand = c(0, 0)) +
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
