options(scipen = 100)
library(viridis)
library(factoextra)
library(tidyverse)
library(ade4)
library(ggplot2)
library(plotly)
# library(rayshader)
library(Biostrings)
library(SynMut)
library(seqinr)
# library(future.apply)
library(ggrepel)
library(ggpubr)
library(ggtree)
library(ggalt)
# library(systemPipeR)
# plan("multiprocess")


# input data --------------------------------------------------------------
meta_data <- read_csv("../data/meta_omsn_tidy_2020_02_20_pan.csv")
names(meta_data)[c(6, 9, 14, 20, 21)] <-
    c("Gene", "virus_species_ori", "Host_ori", "Host", "Virus Species")
cds_omsn <- readDNAStringSet("../data/cleaned_cds_osmn_2020_02_20_pan.fasta")
meta_data$width <- width(cds_omsn)
###delete guangdong pangolin
cds_omsn <- cds_omsn[!grepl("pangolin/Guangdong/P2S/2019", meta_data$cds_name)]
meta_data <- meta_data[!grepl("pangolin/Guangdong/P2S/2019", meta_data$cds_name),]

# Filter only betacov
idx <- meta_data$`Virus Genus` == "Betacoronavirus"
meta_data <- meta_data[idx,]
cds_omsn <- cds_omsn[idx]

count_tmp <- table(meta_data$`Virus Species`)
major_names <- names(count_tmp[count_tmp >= 40])
meta_data$species_major <- meta_data$`Virus Species`
meta_data$species_major[!meta_data$species_major %in% major_names] <-
    "Other species"
meta_data$species_major[grep("unclassified", meta_data$species_major)] <-
    "Other species"
## reoreder
meta_data$Host <- factor(
    meta_data$Host,
    levels = c(
        "avian",
        "bat",
        "camel",
        "feline",
        "human",
        "pangolin",
        "rodent",
        "swine",
        "others"
    )
)

meta_data$strain_name[meta_data$strain_name == "2019-nCoV"] <- "SARS-CoV-2"
meta_data$`Strain Name`[meta_data$`Strain Name` == "2019-nCoV"] <- "SARS-CoV-2"
meta_data$`Virus Species`[meta_data$`Virus Species` == "2019-nCoV"] <- "SARS-CoV-2"
meta_data$Type[meta_data$Type == "2019-nCoV"] <- "SARS-CoV-2"
meta_data$species_major[meta_data$species_major == "2019-nCoV"] <- "SARS-CoV-2"

# seq_mut <- seq_mut_orf_orf1ab
# gene <- "orf1ab"
CA_analysis <- function(seq_mut, gene) {
    stopifnot(gene %in% unique(meta_data$Gene))
    # bind to ori data
    name_tmp = names(seq_mut)
    meta_data_mut <- meta_data[seq_along(seq_mut),]
    meta_data_mut <- meta_data_mut %>% mutate_all(function(x){NA})
    meta_data_mut$cds_name <- name_tmp
    meta_data_mut$id <- names(seq_mut)
    meta_data_mut$organism <- name_tmp
    meta_data_mut$strain_name <- name_tmp
    meta_data_mut$Gene <- gene
    meta_data_mut$gene_symbol <- toupper(strsplit(gene, "")[[1]][1])
    meta_data_mut$`Strain Name` <- name_tmp
    meta_data_mut$virus_species_ori <- name_tmp
    meta_data_mut$`Virus Genus` <- meta_data$`Virus Genus`[1]
    meta_data_mut$`Virus Family` <- meta_data$`Virus Family`[1]
    meta_data_mut$Host_ori <- "mutant"
    meta_data_mut$`GenBank Host` <- "mutant"
    meta_data_mut$Host <- "mutant"
    meta_data_mut$`Virus Species` <- name_tmp
    meta_data_mut$important <- 2
    meta_data_mut$Type <- "Mutant"
    meta_data_mut$species_major <- "Other species" 

    meta_data_mut <- bind_rows(meta_data, meta_data_mut)
    cds_omsn_mut <- c(cds_omsn, seq_mut)
    idx_meta_o <- meta_data_mut$Gene==gene

    df_cu_ori <- as_tibble(get_cu(cds_omsn_mut))
    df_rscu <- as_tibble(get_rscu(cds_omsn_mut))

    ## remove stp codon
    df_rscu <- df_rscu %>% select(-TAA, - TGA, - TAG)
    df_cu <- df_cu_ori %>% select(-TAA, - TGA, - TAG)

    # within and between CA ---------------------------------------------------
    fac_codon <- factor(Biostrings::GENETIC_CODE[colnames(df_cu)])
    ## orf1ab
    ttuco.coa_o <- dudi.coa(t(df_cu[idx_meta_o,]), scannf = FALSE, nf = 5)
    ttuco.wca_o <- wca(ttuco.coa_o, fac_codon, scan = FALSE, nf = 5)
    ttuco.bca_o <- bca(ttuco.coa_o, fac_codon, scan = FALSE, nf = 5)
    ### variability at the synonymous level & at the amino acid level
    # 100 * sum(ttuco.wca_o$eig) / sum(ttuco.coa_o$eig)
    # 100 * sum(ttuco.bca_o$eig) / sum(ttuco.coa_o$eig)
    tuco <- df_cu[idx_meta_o,]
    # (nrow(tuco) - 1) * (ncol(tuco) - 1) / sum(tuco) -> exptoti
    # (nrow(tuco) - 1) * (ncol(tuco) - length(levels(fac_codon))) / sum(tuco) ->
    #     exptotiw
    # (nrow(tuco) - 1) * (length(levels(fac_codon)) - 1) / sum(tuco) -> exptotib
    # 100 * (sum(ttuco.wca_o$eig) - exptotiw) / (sum(ttuco.coa_o$eig) - exptoti)
    # 100 * (sum(ttuco.bca_o$eig) - exptotib) / (sum(ttuco.coa_o$eig) - exptoti)

    
    # Synonymous codon usage (WCA) --------------------------------------------
    ## F1
    ### Coding sequences point of view

    F1 <- rep(NA, nrow(meta_data_mut))
    F1[idx_meta_o] <- ttuco.wca_o$co[, 1]
    F2 <- rep(NA, nrow(meta_data_mut))
    F2[idx_meta_o] <- ttuco.wca_o$co[, 2]
    F3 <- rep(NA, nrow(meta_data_mut))
    F3[idx_meta_o] <- ttuco.wca_o$co[, 3]
    tmp_df <- bind_cols(meta_data_mut, F1 = F1, F2 = F2, F3 = F3)
    tmp_df$Type <- factor(tmp_df$Type, levels = c("Others", "SARS-CoV-2", "Reference", "Mutant"))
    idx_o <- which(tmp_df$id == "MN908947" & tmp_df$Gene == gene)
    tmp_df <- tmp_df %>% mutate(dist_o = sqrt((F1 - F1[idx_o]) ^ 2 + (F2 - F2[idx_o]) ^ 2))
    tmp_df$Gene <- factor(tmp_df$Gene, levels = c("orf1ab", "spike", "membrane", "nucleocapsid"))
    tmp_df <- tmp_df %>% filter(Gene==gene)

    kmodel <- tmp_df %>% mutate(l_num = 1:n()) %>%
        group_by(Gene) %>% select(l_num, F1, F2) %>% nest() %>%
        mutate(kmodel = lapply(data, function(x) {
            tmp <- kmeans(x[, 2:3], centers = 7, nstart = 2000, iter.max = 1000)
            factor(tmp$cluster)
        })) %>% unnest(c(data, kmodel)) %>% arrange(l_num) %>%
        select(kmodel)

    tmp_df$species_major <- factor(tmp_df$species_major,
        levels = c(unique(tmp_df$species_major)[2:7],unique(tmp_df$species_major)[1]))

    plot_wca_host <- bind_cols(tmp_df, kmodel) %>%
    ggplot() +
        geom_encircle(aes(x = F1, y = F2, group = kmodel), linetype = 2, alpha = 0.6,
            expand = 0.01, spread = 0.001) +
        geom_point(aes(x = F1, y = F2, color = Host, shape = Host),
                alpha = 0.5) +
        scale_shape_manual(name = "Host", values = seq_along(unique(tmp_df$Host))) +
        facet_wrap(vars(Gene)) +
        geom_text_repel(aes(x = F1, y = F2, label = "SARS-CoV-2"),
                        data = filter(tmp_df, id == "MN908947"),
                        nudge_x = 0.5, nudge_y = 0.5) +
        geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
                        data = tmp_df %>% filter(Host=="mutant"),
                        nudge_x = -1, direction = "y", size = 3) +
        scale_color_viridis_d() +
        ggtitle("WCA (synonymous codon usage)")
    ggsave(paste0("../results/wca_mutant_", gene, "_Host.pdf"), plot=plot_wca_host)

    plot_wca_species <- bind_cols(tmp_df, kmodel) %>%
        # filter(Gene == "membrane", species_major == "SARS-CoV-2")  %>% 
    ggplot() +
        geom_point(aes(x = F1, y = F2, color = species_major),
                alpha = 0.8) +
        scale_shape_manual(name = "Virus species", values = 1:7) +
        facet_wrap(vars(Gene)) +
        geom_encircle(aes(x = F1, y = F2, group = species_major, fill = species_major), 
            linetype = 2, alpha = 0.5, expand = 0.01, spread = 0.001) +
        geom_text_repel(aes(x = F1, y = F2, label = "SARS-CoV-2"),
                        data = filter(tmp_df, id == "MN908947"),
                        nudge_x = 0.5, nudge_y = 0.5) +
        geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
                        data = tmp_df %>% filter(Host=="mutant"),
                        nudge_x = -1, direction = "y", size = 3) +
        scale_color_viridis_d(name = "Virus species") +
        scale_fill_viridis_d(name = "Virus species") +    
        ggtitle("WCA (synonymous codon usage)")
    ggsave(paste0("../results/wca_mutant_", gene, "_Species.pdf"), plot=plot_wca_species)

    tmp_df <- tmp_df %>% mutate(Type = ifelse(important == 2, "Mutant", "WT"))

    plot_wca_mut <-
        plot_ly(
            filter(tmp_df, Gene == gene),
            x = ~F1,
            y = ~F2,
            z = ~F3,
            color = filter(tmp_df, Gene == gene) %>% .$Host,
            text = filter(tmp_df, Gene == gene) %>% .$strain_name,
            colors = viridis_pal(option = "D")(6),
            marker = list(size = 5, opacity = 0.5)
        ) %>%
        add_markers()
    p <- plotly_build(plot_wca_mut)
    sapply(seq_along(p$x$data), function(i){
        x = p$x$data[[i]]
        tmp = grepl("Mutant", x$text)
        p$x$data[[i]]$marker$symbol <<- ifelse(tmp, "square", "circle")
    })
    htmlwidgets::saveWidget(p, paste0("../results/3d_mutant_", gene, ".html"))

}















# # Amino acid usage (BCA) --------------------------------------------------

# ### Coding sequences point of view
# #### F1 F2
# F1 <- rep(NA, nrow(meta_data))
# F1[idx_meta_o] <- ttuco.bca_o$co[, 1]
# F1[idx_meta_s] <- ttuco.bca_s$co[, 1]
# F1[idx_meta_m] <- ttuco.bca_m$co[, 1]
# F1[idx_meta_n] <- ttuco.bca_n$co[, 1]
# F2 <- rep(NA, nrow(meta_data))
# F2[idx_meta_o] <- ttuco.bca_o$co[, 2]
# F2[idx_meta_s] <- ttuco.bca_s$co[, 2]
# F2[idx_meta_m] <- ttuco.bca_m$co[, 2]
# F2[idx_meta_n] <- ttuco.bca_n$co[, 2]
# F3 <- rep(NA, nrow(meta_data))
# F3[idx_meta_o] <- ttuco.bca_o$co[, 3]
# F3[idx_meta_s] <- ttuco.bca_s$co[, 3]
# F3[idx_meta_m] <- ttuco.bca_m$co[, 3]
# F3[idx_meta_n] <- ttuco.bca_n$co[, 3]
# tmp_df <- bind_cols(meta_data, F1 = F1, F2 = F2, F3 = F3)

# tmp_df$Type <- factor(tmp_df$Type, levels = c("Others", "SARS-CoV-2", "Reference"))
# idx_o <- which(tmp_df$id == "MN908947" & tmp_df$Gene == "orf1ab")
# idx_s <- which(tmp_df$id == "MN908947" & tmp_df$Gene == "spike")
# idx_m <- which(tmp_df$id == "MN908947" & tmp_df$Gene == "membrane")
# idx_n <- which(tmp_df$id == "MN908947" & tmp_df$Gene == "nucleocapsid")
# tmp_df <- tmp_df %>% mutate(dist_o = sqrt((F1 - F1[idx_o]) ^ 2 + (F2 - F2[idx_o]) ^ 2))
# tmp_df <- tmp_df %>% mutate(dist_s = sqrt((F1 - F1[idx_s]) ^ 2 + (F2 - F2[idx_s]) ^ 2))
# tmp_df <- tmp_df %>% mutate(dist_m = sqrt((F1 - F1[idx_m]) ^ 2 + (F2 - F2[idx_m]) ^ 2))
# tmp_df <- tmp_df %>% mutate(dist_n = sqrt((F1 - F1[idx_n]) ^ 2 + (F2 - F2[idx_n]) ^ 2))
# tmp_df$Gene <- factor(tmp_df$Gene, levels = c("orf1ab", "spike",
#                                               "membrane", "nucleocapsid"))

# kmodel <- tmp_df %>% mutate(l_num = 1:n()) %>%
#     group_by(Gene) %>% select(l_num, F1, F2) %>% nest() %>%
#     mutate(kmodel = lapply(data, function(x) {
#         tmp <- kmeans(x[, 2:3], centers = 7, nstart = 2000, iter.max = 1000)
#         factor(tmp$cluster)
#     })) %>% unnest(c(data, kmodel)) %>% arrange(l_num) %>%
#     select(kmodel)


# tmp_df$species_major <- factor(tmp_df$species_major,
#     levels = c(unique(tmp_df$species_major)[2:7],unique(tmp_df$species_major)[1]))

# (plot_bca_host <- bind_cols(tmp_df, kmodel) %>%
# ggplot() +
#     geom_point(aes(x = F1, y = F2, color = Host, shape = kmodel),
#                alpha = 0.5) +
#     geom_encircle(aes(x = F1, y = F2, group = kmodel), linetype = 2, alpha = 0.6,
#         expand = 0.01, spread = 0.001) +
#     facet_wrap(vars(Gene)) +
#     scale_shape_manual(name = "Cluster", values = 1:7) +
#     geom_text_repel(aes(x = F1, y = F2, label = "SARS-CoV-2"),
#                     data = filter(tmp_df, id == "MN908947"),
#                     nudge_x = 0.5, nudge_y = 0.5) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_o],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_m],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_n],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     scale_color_viridis_d() +
#     ggtitle("B. BCA (amino acid usage)"))

# (plot_bca_species <- bind_cols(tmp_df, kmodel) %>%
# ggplot() +
#     geom_point(aes(x = F1, y = F2, color = species_major),
#                alpha = 0.8) +
#     # scale_shape_manual(name = "Virus species", values = 1:7) +
#     facet_wrap(vars(Gene)) +
#     geom_encircle(aes(x = F1, y = F2, group = species_major, fill = species_major), 
#         linetype = 2, alpha = 0.5, expand = 0.01, spread = 0.001) +
#     geom_text_repel(aes(x = F1, y = F2, label = "SARS-CoV-2"),
#                     data = filter(tmp_df, id == "MN908947"),
#                     nudge_x = 0.5, nudge_y = 0.5) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_o],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_s],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_m],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     geom_text_repel(aes(x = F1, y = F2, label = strain_name), alpha = 0.6,
#                     data = tmp_df %>% filter(Gene == tmp_df$Gene[idx_n],
#                                              important == 1),
#                     nudge_x = -1, direction = "y", size = 3) +
#     scale_color_viridis_d(name = "Virus Species") +
#     scale_fill_viridis_d(name = "Virus Species") +
#     ggtitle("B. BCA (amino acid usage)"))


# # plotting together -------------------------------------------------------

# ggarrange(plot_wca_host, plot_bca_host, common.legend = T, ncol = 1,
#           legend = "right")
# ggsave("../results/Figure_4.tiff", width = 10, height = 15,
#        device = "tiff", dpi = 300, scale = 0.9,
#        compress = "lzw")
# # ggsave("../results/Figure_3.svg", width = 10, height = 15, dpi = 300, scale = 0.9)

# ggarrange(plot_wca_species, plot_bca_species, common.legend = T, nrow = 2,
#           legend = "right")
# ggsave("../results/Figure_5.tiff",
#        width = 14, height = 15,
#        device = "tiff", dpi = 300, scale = 0.9,
#        compress = "lzw")
# # ggsave("../results/Figure_4.svg", width = 12, height = 15, dpi = 300, scale = 0.9)


# ## 3D
# plot_bca_spike <-
#     plot_ly(
#         filter(tmp_df, Gene == "spike"),
#         x = ~F1,
#         y = ~F2,
#         z = ~F3,
#         color = filter(tmp_df, Gene == "spike") %>% .$Host,
#         text = filter(tmp_df, Gene == "spike") %>% .$strain_name,
#         colors = viridis_pal(option = "D")(6),
#         marker = list(size = 5, opacity = 0.5)
#     ) %>%
#     add_markers()

# htmlwidgets::saveWidget(plot_wca_spike, "../results/plot_bca_spike.html")

# plot_bca_membrane <-
#     plot_ly(
#         filter(tmp_df, Gene == "membrane"),
#         x = ~F1,
#         y = ~F2,
#         z = ~F3,
#         color = filter(tmp_df, Gene == "membrane") %>% .$Host,
#         text = filter(tmp_df, Gene == "membrane") %>% .$strain_name,
#         colors = viridis_pal(option = "D")(6),
#         marker = list(size = 5, opacity = 0.5)
#     ) %>%
#     add_markers()

# htmlwidgets::saveWidget(plot_bca_membrane, "../results/plot_bca_membrane.html")

# plot_bca_orf1ab <-
#     plot_ly(
#         filter(tmp_df, Gene == "orf1ab"),
#         x = ~F1,
#         y = ~F2,
#         z = ~F3,
#         color = filter(tmp_df, Gene == "orf1ab") %>% .$Host,
#         text = filter(tmp_df, Gene == "orf1ab") %>% .$strain_name,
#         colors = viridis_pal(option = "D")(6),
#         marker = list(size = 5, opacity = 0.5)
#     ) %>%
#     add_markers()

# htmlwidgets::saveWidget(plot_bca_orf1ab, "../results/plot_bca_orf1ab.html")

# plot_bca_nucleocapsid <-
#     plot_ly(
#         filter(tmp_df, Gene == "nucleocapsid"),
#         x = ~F1,
#         y = ~F2,
#         z = ~F3,
#         color = filter(tmp_df, Gene == "nucleocapsid") %>% .$Host,
#         text = filter(tmp_df, Gene == "nucleocapsid") %>% .$strain_name,
#         colors = viridis_pal(option = "D")(6),
#         marker = list(size = 5, opacity = 0.5)
#     ) %>%
#     add_markers()

# htmlwidgets::saveWidget(plot_bca_nucleocapsid, "../results/plot_bca_nucleocapsid.html")

# # scree plot --------------------------------------------------------------

# df_eig_all <- tibble(`Global codon usage (orf1ab)` = tuco.coa_o$eig[1:15],
#                      `Global codon usage (spike)` = tuco.coa_s$eig[1:15],
#                      `Global codon usage (membrane)` = tuco.coa_m$eig[1:15],
#                      `Global codon usage (nucleocapsid)` = tuco.coa_n$eig[1:15],
#                      `Synonymous codon usage (orf1ab)` = ttuco.wca_o$eig[1:15],
#                      `Synonymous codon usage (spike)` = ttuco.wca_s$eig[1:15],
#                      `Synonymous codon usage (membrane)` = ttuco.wca_m$eig[1:15],
#                      `Synonymous codon usage (nucleocapsid)` = ttuco.wca_n$eig[1:15],
#                      `Amino acid usage (orf1ab)` = ttuco.bca_o$eig[1:15],
#                      `Amino acid usage (spike)` = ttuco.bca_s$eig[1:15],
#                      `Amino acid usage (membrane)` = ttuco.bca_m$eig[1:15],
#                      `Amino acid usage (nucleocapsid)` = ttuco.bca_n$eig[1:15])
# df_eig_all <- df_eig_all %>% mutate_all(function(x) {
#     100 * x / sum(x)
# })
# df_eig_all <- bind_cols(df_eig_all, id = 1:15)

# df_eig_all <- df_eig_all %>% gather(key = group, value = value, - id, factor_key = T)

# ggplot(df_eig_all) +
#     geom_col(aes(x = id, y = value)) +
#     facet_wrap(vars(group), ncol = 4)
# ggsave("../results/Figure_S8.tiff",
#        width = 12, height = 12,
#        device = "tiff", dpi = 300, scale = 0.8,
#        compress = "lzw")


# # acession ID -------------------------------------------------------------
# write_csv(meta_data, "../results/meta_data_beta_cov_id.csv")


# # codon distance ----------------------------------------------------------
# which(meta_data$Type == "SARS-CoV-2" & meta_data$Gene == "spike")
# df_cu[3001,]
# which(grepl("402131", meta_data$id) & meta_data$Gene == "spike")
# df_cu[3074,]
# which(grepl("173883", meta_data$id) & meta_data$Gene == "spike")
# df_cu[2682,]

# sum((df_rscu[3001, 4] - df_rscu[2682, 4]) ^ 2)
# sum((df_rscu[3001, 4] - df_rscu[3074, 4]) ^ 2)

# fviz_cos2(tuco.coa_s, choice = "col", axes = 1:2)
# fviz_contrib(tuco.coa_s, choice = "col", axes = 1:2)
