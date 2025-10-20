# Adapated from: https://github.com/Leo-Poon-Lab/SARS-CoV-2-sVNT-diversity/blob/main/Antigenic_diversity.R
# Quantify antigenic diversity on different RBDs for different groups,
# and comparison between groups.

library(readxl)
library(tidyverse)
library(parallel) # mclapply might still be used for iterating, or can be replaced by lapply
library(ggrepel)
library(patchwork)
library(ggplot2)

# read data
df_inhibition <- read_excel("data/Antigen diversity vs magnitude.xlsx") # Ensure this path is correct
df_inhibition <- df_inhibition %>% pivot_longer(-c(group, sample), names_to = "RBD")
if(all(df_inhibition$value < 1, na.rm = TRUE)) {
  df_inhibition$value <- df_inhibition$value * 100 # Convert to percentage if all values are < 1
}
df_inhibition$Response <- ifelse(df_inhibition$value<20, "negative", "positive")
df_inhibition$group_sample <- paste0(df_inhibition$group, "_", df_inhibition$sample)

RBDs_all <- unique(df_inhibition$RBD)
RBD_SARSCoV2 <- c(
  "WT SARS-CoV-2",
  "Alpha",
  "Beta",
  "Delta",
  "Delta Plus",
  "Lambda",
  "Omicron BA.1",
  "Omicron BA.2",
  "Omicron BA.5",
  "XBB.1",
  "EG.5",
  "JN.1"
)
RBD_Sarbecovirus <- c(
  "SARS-CoV-1", "GD-1", "WIV-1",
  "Rs2018B", "RsSHC014", "LyRall",
  "Pangolin MERS MjHKU4r-CoV", "Human MERS-CoV",
  "Camelid MERS-CoV", "Bat MERS BtCoV-422"
)
stopifnot(sort(c(RBD_SARSCoV2, RBD_Sarbecovirus)) == sort(RBDs_all)) # Ensure all RBDs are accounted for
df_inhibition$RBD_group <- ifelse(df_inhibition$RBD %in% RBD_SARSCoV2, "SARS-CoV-2", "Sarbecovirus")

groups_all <- unique(df_inhibition$group)
group_sample_all <- unique(df_inhibition$group_sample)

dir.create("results/antigenic_diversity", showWarnings = FALSE) # Added showWarnings = FALSE

# When analyzing antigenic diversity, we have two questions to be answered:
# 1. Which groups have higher magnitude responses?
# 2. Which groups have boarder antibody responses?
# (Definitions remain the same)

# Modified function to calculate metrics directly for a group
calculate_group_metrics <- function(data_group){
  # data_group is df_inhibition_i (a subset of df_inhibition for a specific group and panel)
  
  ### 1. Magnitude of responses
  mean_res <- mean(data_group$value, na.rm = TRUE)
  
  positive_responses <- data_group$value[data_group$Response=="positive"]
  mean_pos_res <- if(length(positive_responses) > 0) mean(positive_responses, na.rm = TRUE) else NA
  
  ### 2. Breadth of antibody response (Diversity)
  hs <- NA # Default NA
  hpi <- NA # Default NA

  if(nrow(data_group) > 0) {
    df_tmp <- data_group
    df_tmp$RBD2 <- df_tmp$RBD
    df_tmp$RBD2[df_tmp$Response=="negative"] <- "Negatives"
    
    # Consider all unique RBDs present in the original full panel plus "Negatives"
    # This ensures diversity is comparable even if some RBDs are missing in a subset
    # all_possible_responses <- c(RBDs_all, "Negatives") # If RBDs_all is globally available and relevant
    # df_tmp$RBD2 <- factor(df_tmp$RBD2, levels = all_possible_responses)
    
    df_tmp_freq <- df_tmp %>% group_by(RBD2) %>% summarise(N=n(), .groups = 'drop')
    N_total_group <- nrow(df_tmp) # Total samples in this specific group for this panel
    
    if (N_total_group > 0 && nrow(df_tmp_freq) > 0) {
      df_tmp_freq$p <- df_tmp_freq$N / N_total_group
      
      # Shannon entropy
      p_for_hs <- df_tmp_freq$p[df_tmp_freq$p > 0]
      hs <- if(length(p_for_hs) > 0) sum(-p_for_hs * log(p_for_hs)) else 0
      
      # HPI (Simpson's diversity index like), "antigenic diveristy (pi)" which is simialr to nucleotide diversity pi, to estimate the diversity
      if (N_total_group < 2) {
        hpi <- 0 # Or NA, if preferred for N < 2
      } else {
        sum_n_n_minus_1 <- sum(df_tmp_freq$N * (df_tmp_freq$N - 1), na.rm = TRUE)
        denominator_hpi <- N_total_group * (N_total_group - 1)
        if(denominator_hpi > 0){
            hpi <- 1 - (sum_n_n_minus_1 / denominator_hpi)
        } else { # Should be caught by N_total_group < 2, but as a safeguard
            hpi <- 0
        }
      }
    } else { # If N_total_group is 0 or df_tmp_freq is empty (should not happen if N_total_group > 0)
        hs <- 0 # Or NA
        hpi <- 0 # Or NA
    }
  }
  
  return(list(mean_res=mean_res, mean_pos_res=mean_pos_res, hs=hs, hpi=hpi))
}

# --- Simplified Metric Calculations ---
# Full spectrum/panel
list_metrics_full <- lapply(group_sample_all, function(group_sample_i){
  df_inhibition_i <- df_inhibition %>% filter(group_sample == group_sample_i)
  metrics <- if(nrow(df_inhibition_i) == 0) {
    list(mean_res=NA, mean_pos_res=NA, hs=NA, hpi=NA)
  } else {
    calculate_group_metrics(df_inhibition_i)
  }
  return(c(group_sample=group_sample_i, e_mean_res=metrics$mean_res, e_mean_pos_res=metrics$mean_pos_res, e_hs=metrics$hs, e_hpi=metrics$hpi))
})
df_summary_metrics_full <- as_tibble(do.call(rbind, list_metrics_full))
numeric_cols <- c("e_mean_res", "e_mean_pos_res", "e_hs", "e_hpi")
df_summary_metrics_full[numeric_cols] <- lapply(df_summary_metrics_full[numeric_cols], as.numeric)
write_tsv(df_summary_metrics_full, "results/antigenic_diversity/df_summary_metrics_full.tsv")

# Partial spectrum/panel metrics by RBD_group
# SARS-CoV-2 panel
# Partial spectrum/panel metrics by group_sample instead of just group
# SARS‐CoV‐2 panel
df_inhibition_partial_scov2 <- df_inhibition %>% filter(RBD_group == "SARS-CoV-2")
list_metrics_scov2 <- lapply(group_sample_all, function(gs){
  df_i <- df_inhibition_partial_scov2 %>% filter(group_sample == gs)
  metr <- if (nrow(df_i) == 0) {
    list(mean_res=NA, mean_pos_res=NA, hs=NA, hpi=NA)
  } else {
    calculate_group_metrics(df_i)
  }
  c(group_sample = gs,
    e_mean_res    = metr$mean_res,
    e_mean_pos_res= metr$mean_pos_res,
    e_hs          = metr$hs,
    e_hpi         = metr$hpi)
})
df_summary_metrics_scov2 <- as_tibble(do.call(rbind, list_metrics_scov2))
df_summary_metrics_scov2[numeric_cols] <- lapply(df_summary_metrics_scov2[numeric_cols], as.numeric)
write_tsv(df_summary_metrics_scov2, "results/antigenic_diversity/df_summary_metrics_scov2.tsv")

# Sarbecovirus panel
df_inhibition_partial_sarbecov <- df_inhibition %>% filter(RBD_group == "Sarbecovirus")
list_metrics_sarbecov <- lapply(group_sample_all, function(gs){
  df_i <- df_inhibition_partial_sarbecov %>% filter(group_sample == gs)
  metr <- if (nrow(df_i) == 0) {
    list(mean_res=NA, mean_pos_res=NA, hs=NA, hpi=NA)
  } else {
    calculate_group_metrics(df_i)
  }
  c(group_sample = gs,
    e_mean_res    = metr$mean_res,
    e_mean_pos_res= metr$mean_pos_res,
    e_hs          = metr$hs,
    e_hpi         = metr$hpi)
})
df_summary_metrics_sarbecov <- as_tibble(do.call(rbind, list_metrics_sarbecov))
df_summary_metrics_sarbecov[numeric_cols] <- lapply(df_summary_metrics_sarbecov[numeric_cols], as.numeric)
write_tsv(df_summary_metrics_sarbecov, "results/antigenic_diversity/df_summary_metrics_sarbecov.tsv")

# --- Plotting ---
df_meta <- df_inhibition %>% select(group, sample, group_sample) %>% distinct()
# parent_groups_all <- unique(df_meta$parent_group3[!is.na(df_meta$parent_group3)]) # Or your predefined order

# Full spectrum/panel plots
df_plot_summary_full <- read_tsv("results/antigenic_diversity/df_summary_metrics_full.tsv")
df_plot_summary_full <- left_join(df_plot_summary_full, df_meta, by = "group_sample") # Assumes df_meta is loaded

# Define colors for each group
colors_t_full <- c(
  "PBS" = "#000000",               # Black
  "PBS+PBS" = "#808080",           # Grey
  "cb1_e3 PFU" = "#4682B4",        # Steel Blue (lighter)
  "cb1_e4 PFU" = "#104E8B",        # Medium Blue
  "cb1_e5 PFU" = "#00008B",        # Dark Blue
  "cb1 + RG_WT SARS-CoV-2" = "#800020", # Burgundy/Wine
  "cb1 + VoC_EG.5" = "#228B22",    # Forest Green
  "cb1 + VoC_JN.1" = "#483D8B",    # Dark Purple/Slate Blue
  "PBS + VoC_JN.1" = "#8B4513"     # Dark Brown/Saddle Brown
)

df_plot_summary_full$group <- factor(df_plot_summary_full$group, levels=unique(df_plot_summary_full$group))

p1 <- ggplot(df_plot_summary_full, aes(x=group, y=e_mean_res, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Average %Inhibition")+
  scale_fill_manual(name="Group", values=colors_t_full)+
  scale_color_manual(name="Group", values=colors_t_full, guide="none")+
  theme_bw()+ xlab("")+ggtitle("All antigens") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/Mean_inhibition_responses_full_direct_calc.pdf", width = 8, height=5, plot = p1)

p2 <- ggplot(df_plot_summary_full, aes(x=group, y=e_hs, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Shannon entropy")+
  scale_fill_manual(name="Group", values=colors_t_full)+
  scale_color_manual(name="Group", values=colors_t_full, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("All antigens") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hs_full_direct_calc.pdf", width = 8, height=5, plot = p2)

p3 <- ggplot(df_plot_summary_full, aes(x=group, y=e_hpi, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab(expression("RBD cross-reactivity ("~pi~")"))+
  scale_fill_manual(name="Group", values=colors_t_full)+
  scale_color_manual(name="Group", values=colors_t_full, guide="none")+
  theme_bw()+ xlab("")+ggtitle("All antigens") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hpi_full_direct_calc.pdf", width = 8, height=5, plot = p3)

# Combined plots for full spectrum
p_out <- (p1+ggtitle("A"))/(p3+ggtitle("B"))/(p2+ggtitle("C")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("results/antigenic_diversity/Antigenic_response_combined_add_hs_full_spectrum.pdf", width = 8, height=12)

# Partial spectrum/panel plots (SARS-CoV-2 and Sarbecoviruses)
df_summary_scov2_plot <- read_tsv("results/antigenic_diversity/df_summary_metrics_scov2.tsv")
df_summary_scov2_plot$RBDs <- "SARS-CoV-2"
df_summary_sarbecov_plot <- read_tsv("results/antigenic_diversity/df_summary_metrics_sarbecov.tsv")
df_summary_sarbecov_plot$RBDs <- "Sarbecoviruses"

# Combine partial datasets
df_plot_summary_partial <- bind_rows(df_summary_scov2_plot, df_summary_sarbecov_plot)
df_plot_summary_partial <- left_join(df_plot_summary_partial, df_meta, by = "group_sample")
df_plot_summary_partial$group <- factor(df_plot_summary_partial$group, levels=unique(df_plot_summary_partial$group))

# Define colors for partial spectrum plots similar to full spectrum
colors_t_partial <- colors_t_full

# Plot for mean inhibition - SARS-CoV-2
p1_scov2 <- ggplot(df_plot_summary_partial %>% filter(RBDs == "SARS-CoV-2"), aes(x=group, y=e_mean_res, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Average %Inhibition")+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("SARS-CoV-2") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/Mean_inhibition_responses_scov2_direct_calc.pdf", width = 8, height=5, plot = p1_scov2)

# Plot for mean inhibition - Sarbecoviruses
p1_sarbeco <- ggplot(df_plot_summary_partial %>% filter(RBDs == "Sarbecoviruses"), aes(x=group, y=e_mean_res, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Average %Inhibition")+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("Sarbecoviruses") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/Mean_inhibition_responses_sarbeco_direct_calc.pdf", width = 8, height=5, plot = p1_sarbeco)

# Plot for Shannon entropy - SARS-CoV-2
p2_scov2 <- ggplot(df_plot_summary_partial %>% filter(RBDs == "SARS-CoV-2"), aes(x=group, y=e_hs, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Shannon entropy")+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("SARS-CoV-2") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hs_scov2_direct_calc.pdf", width = 8, height=5, plot = p2_scov2)

# Plot for Shannon entropy - Sarbecoviruses
p2_sarbeco <- ggplot(df_plot_summary_partial %>% filter(RBDs == "Sarbecoviruses"), aes(x=group, y=e_hs, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab("Shannon entropy")+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("Sarbecoviruses") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hs_sarbeco_direct_calc.pdf", width = 8, height=5, plot = p2_sarbeco)

# Plot for HPI - SARS-CoV-2
p3_scov2 <- ggplot(df_plot_summary_partial %>% filter(RBDs == "SARS-CoV-2"), aes(x=group, y=e_hpi, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab(expression("RBD cross-reactivity ("~pi~")"))+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("SARS-CoV-2") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hpi_scov2_direct_calc.pdf", width = 8, height=5, plot = p3_scov2)

# Plot for HPI - Sarbecoviruses
p3_sarbeco <- ggplot(df_plot_summary_partial %>% filter(RBDs == "Sarbecoviruses"), aes(x=group, y=e_hpi, fill=group)) +
  geom_point(aes(color=group), size=3, shape=21, stroke=0.6) +
  scale_x_discrete(expand=c(0,1.5), drop=F) + ylab(expression("RBD cross-reactivity ("~pi~")"))+
  scale_fill_manual(name="Group", values=colors_t_partial)+
  scale_color_manual(name="Group", values=colors_t_partial, guide="none")+
  theme_bw()+ xlab("")+ ggtitle("Sarbecoviruses") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 1), legend.position = "top")+
  guides(fill=guide_legend(ncol=3,byrow=F))+
  NULL
ggsave("results/antigenic_diversity/hpi_sarbeco_direct_calc.pdf", width = 8, height=5, plot = p3_sarbeco)

# Combined plots - side by side arrangement
p_out1 <- (p1)/(p1_scov2)/(p1_sarbeco) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("results/antigenic_diversity/Mean_inhibition_combined.pdf", width = 8, height=12, plot = p_out1)

p_out3 <- (p3)/(p3_scov2)/(p3_sarbeco) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("results/antigenic_diversity/hpi_combined.pdf", width = 8, height=12, plot = p_out3)

p_out4 <- (p2)/(p2_scov2)/(p2_sarbeco) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("results/antigenic_diversity/hs_combined.pdf", width = 8, height=12, plot = p_out4)

# --- 2D Illustration (Response magnitude vs diversity) ---
colors_for_2d <- colors_t_partial 

plot_2d_no_ci <- function(df_summary, x_var_metric, y_var_metric, color_by_var="group", label_points_by="group") {
  
  x_col_name <- paste0("e_", x_var_metric)
  y_col_name <- paste0("e_", y_var_metric)
  
  x_axis_label <- if(grepl("pos_res", x_var_metric)) "Average %inhibition (positive responses)" else "Average %inhibition (all responses)"
  y_axis_label <- if(y_var_metric == "hpi") expression("RBD cross-reactivity ("~pi~")") else if(y_var_metric == "hs") "Shannon Entropy" else y_var_metric

  # Ensure colors_for_2d and parent_groups_all are correctly scoped or passed.
  # Using colors_for_2d from the parent environment for simplicity here.
  
  # Create labels for legend dynamically based on available unique parent_group3 and colors
  unique_parent_groups <- unique(df_summary[[color_by_var]][!is.na(df_summary[[color_by_var]])])
  
  # Make sure parent_groups_all is defined and includes all unique_parent_groups for correct labeling.
  # If parent_groups_all is not defined or doesn't match, legend labels might be generic.
  # For simplicity, we'll assume parent_groups_all is correctly defined to match unique_parent_groups.
  
  legend_labels <- waiver() # Default

  p <- ggplot(df_summary, aes_string(x=x_col_name, y=y_col_name)) +
    geom_point(aes_string(color=color_by_var, fill=color_by_var), alpha=0.8, size=2.5) +
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16), guide = "none") + # open_circle is logical
    geom_text_repel(aes_string(color=color_by_var, label=label_points_by), 
                    bg.color = "white", bg.r = 0.05, size=1.5, 
                    segment.size=0.3, segment.alpha=0.9,
                    arrow = arrow(length = unit(0.01, "npc")), 
                    point.padding=0.3, box.padding=0.5, # Adjusted padding
                    show.legend=FALSE, max.overlaps = Inf, force=10, min.segment.length = 0.1)+
    scale_fill_manual(name="Group", values=colors_for_2d, labels = legend_labels, drop=!any(is.na(legend_labels)))+
    scale_color_manual(name="Group", values=colors_for_2d, labels = legend_labels, drop=!any(is.na(legend_labels)))+
    ylab(y_axis_label)+
    xlab(x_axis_label)+
    theme_bw()+
    theme(legend.position = "top")+
    guides(fill=guide_legend(ncol=3,byrow=F),color=guide_legend(ncol=3,byrow=F))+
    NULL
  return(p)
}

p_scov2_all_res <- plot_2d_no_ci(df_plot_summary_partial %>% filter(RBDs=="SARS-CoV-2"), "mean_res", "hpi")
ggsave("results/antigenic_diversity/2D_res_hpi_scov2_no_ci.pdf", width=8, height=6)
p_sarbeco_all_res <- plot_2d_no_ci(df_plot_summary_partial %>% filter(RBDs!="SARS-CoV-2"), "mean_res", "hpi")
ggsave("results/antigenic_diversity/2D_res_hpi_sarbeco_no_ci.pdf", width=8, height=6)

p_2d_out_all_res <- (p_scov2_all_res+ggtitle("A (SARS-CoV-2)"))/(p_sarbeco_all_res+ggtitle("B (Sarbecoviruses)")) + plot_layout(guides="collect")&theme(legend.position = "bottom")
ggsave("results/antigenic_diversity/2D_all_res_hpi_combined_no_ci.pdf", width = 8, height=10)

# # Point to Root Distance Calculation (without CIs)
# # Use df_plot_summary_partial (which is derived from the summary metric files)
# df_dist_no_ci <- df_plot_summary_partial %>%
#   filter(!is.na(e_mean_res) & !is.na(e_hpi)) %>% # Ensure values are not NA
#   # Group by relevant columns to ensure one row per group/RBD combination
#   # 'label' and 'hex_code' should be present from the join with df_meta
#   group_by(RBDs, parent_group3, group, label, hex_code) %>% 
#   summarise(
#     e_mean_res = first(e_mean_res), 
#     e_hpi = first(e_hpi),
#     .groups = 'drop'
#   ) %>%
#   mutate(
#     e_mean_res_scaled = e_mean_res / 100,
#     e_hpi_adj = ifelse(e_hpi < 0, 0, e_hpi), # Adjust negative HPI to 0 if necessary
#     dist_point_root = sqrt(e_mean_res_scaled^2 + e_hpi_adj^2)
#   ) %>%
#   select(RBDs, parent_group3, group, label, hex_code, e_mean_res, e_hpi, dist_point_root)

# # Ensure writexl is loaded if not already: library(writexl)
# writexl::write_xlsx(df_dist_no_ci, "results/antigenic_diversity/df_dist_root_point_no_ci.xlsx")