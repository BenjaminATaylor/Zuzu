process FULLSYNTH_PLOTS {

    publishDir "$params.outdir"

    input:
    val deseq_synths
    val edger_synths
    val wilcox_synths
    val svc_synths


    output:
    path "synth_plot.pdf"

    """
    #!/usr/bin/env Rscript
    library("tidyverse", quietly = TRUE)
    library("reshape2", quietly = TRUE)
    library("viridis")

    # reading function for lists of fullsynth output files
    read_fullsynths = function(inlist, label){
        
        this.inlist <- str_remove_all(inlist,"[\\\\[\\\\] ]") %>% 
            strsplit(split = ",") %>% unlist()
        
        out <- sapply(this.inlist, read.csv) %>% 
            t() %>% 
            `row.names<-`(NULL) %>% 
            data.frame() %>%
            mutate_at(c("power","FDP","normAUC","samplenum","depth"),as.numeric) %>%
            mutate(method = label) %>% 
            melt(id.vars = c("method","samplenum", "depth"))

        return(out) 
        
    }

    deseq.synth.input = read_fullsynths("$deseq_synths", "DESeq2")
    edger.synth.input = read_fullsynths("$edger_synths", "edgeR")
    wilcox.synth.input = read_fullsynths("$wilcox_synths", "Wilcoxon")
    # Only include ML outputs if specified by user
    if("$params.mlstep" == "true"){
        svc.synth.input = read_fullsynths("$svc_synths", "SVC")

        gg.synth.input = rbind(deseq.synth.input, 
                            edger.synth.input,
                            wilcox.synth.input,
                            svc.synth.input) %>% 
            mutate(samplenum = as.factor(samplenum)) %>%
            mutate(depth = as.factor(depth)) %>%
            group_by(method,samplenum,depth,variable) %>% 
            summarize(mean=mean(value)) %>%
            mutate(mean = signif(mean,3))
    } else {
        gg.synth.input = rbind(deseq.synth.input, 
                        edger.synth.input,
                        wilcox.synth.input) %>% 
        mutate(samplenum = as.factor(samplenum)) %>%
        mutate(depth = as.factor(depth)) %>%
        group_by(method,samplenum,depth,variable) %>% 
        summarize(mean=mean(value)) %>%
        mutate(mean = signif(mean,3))
        
    }


    # Take FDP as inverse (1-FDP), so that higher values are better (matching to power and normAUC)
    gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")] = 1-(gg.synth.input\$mean[which(gg.synth.input\$variable == "FDP")])
    gg.synth.input\$variable = fct_recode(gg.synth.input\$variable, "1-FDP"="FDP")
    
    # Plot as heatmap    
    gg.synth = ggplot(gg.synth.input, aes(x = depth, y = samplenum)) +
        geom_tile(aes(fill = mean),color = "gray20", size=.75, width=1, height = 1)+
        geom_label(aes(label= mean), color = "white",
                    lineheight = 0.75, size = 3.5, alpha = 0.7, fill = "gray10") +
        scale_x_discrete("Sequencing depth",
                        expand = c(0,0)) +
        scale_y_discrete("Replicates per treatment",
                        expand = c(0,0)) +
        scale_fill_viridis(option = "plasma", 
                            limits = c(0,1),
                            expand = c(0,0)) +
        facet_grid(variable~method) +
        theme_bw() +
        theme(panel.grid = element_blank(),
                legend.title = element_text(face = "bold"),
                axis.text = element_text(face = "bold",size =11, color = "grey40"),
                axis.title = element_text(face = "bold", size =12))

    #set save width based on number of methods tried
    widthvar = 10*(length(unique(gg.synth.input\$method)))

    ggsave(gg.synth, 
        filename = "synth_plot.pdf",
        device = "pdf", bg = "transparent",
        width = widthvar, height = 30, units = "cm")
    """
}
