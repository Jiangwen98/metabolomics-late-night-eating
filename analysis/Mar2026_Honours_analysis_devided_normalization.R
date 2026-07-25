# ---------------------------------------------------------------------------
# NOTE ON THIS PUBLIC COPY
# Trial participant identifiers have been REDACTED and replaced with
# <..._REDACTED...> placeholders. Absolute local paths have been made relative.
# No participant data is included in this repository. See README.md.
# ---------------------------------------------------------------------------

#### the original code already did the baseline-devided AUC nAUC contrasts,
#### now add-in the contrasts of baseline contrasts for thesis - March 2026 ####

####### !!!For this set of data, all postprandial timepoint was normalized by deviding the level by the level at 0min baseline.

##### Basic set-up
if(!file.exists("otp")){dir.create("otp")}
library(tidyverse)

###### 1. Read in and data-cleaning
### Read in (0,0), (1,1) and (0,1) data
# ##* (0,0) data: analyte area only
# qt <- read.delim("mrmkit_output/jiangwen2505_quant_table_4.txt", as.is = T, check.names = F)
###* (1,1) data: analyte height / ISTD height
# qt <- read.delim("mrmkit_output/jiangwen2505_quant_table_2.txt", as.is = T, check.names = F)
###* (0,1) data: analyte area / ISTD area
qt <- read.delim("mrmkit_output/jiangwen2505_quant_table_3.txt", as.is = T, check.names = F)

### Read in (0,0), (1,1) and (0,1) RAW data
###* (0,0) data: analyte area only
# qt_raw <- read.delim("mrmkit_output/jiangwen2505_quant_table_5.txt", as.is = T, check.names = F)
###* (1,1) data: analyte height / ISTD height
# qt_raw <- read.delim("mrmkit_output/jiangwen2505_quant_table_6.txt", as.is = T, check.names = F)
###* (0,1) data: analyte area / ISTD area
qt_raw <- read.delim("mrmkit_output/jiangwen2505_quant_table_7.txt", as.is = T, check.names = F)

### duplicates removal

### Histogram of before and after batch-correction - supplementary

# MWUT between different columns

### Quality control

## S/N filter: S/N ≥ 5 pass
# extract the S/N row
sn_row <- qt_raw[10, ][-(1:3)]
sn_pass <- names(sn_row[which(sn_row >= 5)])

qtc <- qt[, c(1:3, which(colnames(qt) %in% sn_pass))]

## CoV and D-ratio filter: CoV ≤ 30% or D-ratio ≤ 50% pass
# extract CoV and D-ratio rows
cov_row <- qtc[7, ][-(1:3)]
dr_row <- qtc[9, ][-(1:3)]
cov_pass <- names(cov_row[which(cov_row <= 30)])
dr_pass <- names(dr_row[which(dr_row <= 50)])
cvd_pass <- union(cov_pass, dr_pass)

qtc <- qtc[, c(1:3, which(colnames(qtc) %in% cvd_pass))]
### manual review
###column cleaning

### Remove BQC and other useless rows
qtc <- qtc[-(1:10), ]
qtc <- qtc[-which(qtc$type == "BQC"), ]
qtc <- qtc[-3]

### Add-in sample group information
grp <- read.delim("Honours_group_mapping.txt", 
                  as.is = T, check.names = F)
qtc_info <- NULL

for(i in 1:nrow(qtc)){
  row <- qtc[i, ]
  fname <- row$filename
  name <- substr(fname, start = 1, stop = unlist(gregexpr("\\.", fname))-1)
  for(j in 1:nrow(grp)){
    if(grp[j, 1] == name){
      row <- add_column(row, grp[j, ], .after = "batch")
    }
  }
  rm(j)
  
  qtc_info <- rbind(qtc_info, row)
}
rm(i, row, fname, name)

## Remove useless column
qtc_info <- qtc_info[-(1:3)]

### Name cleaning
colnames(qtc_info) <- colnames(qtc_info) %>% str_replace_all(" ", ".") %>% str_replace_all("-", "_") %>% str_replace_all("/", ".or.") %>% str_replace_all(":", ",")

###### create the quant table that only includes 0-min baseline level ######
qtc_bl <- qtc_info[which(qtc_info$Time == 0),]
# 30A is missing, cause the sample is missing
write.csv(qtc_bl, file = paste0("otp/Baseline_raw_value", ".csv"), row.names = T)
###### create the quant table that only includes 0-min baseline level ######



### Now normalize all postprandial timepoint with being devided by 0min
### add-in the placeholder for missing samples - except for 30A, this is entirely missing, because 30A0 is already NA...
a <- c("4C45", 4, "MN", 45, rep(NA, (ncol(qtc_info)-4)))
b <- c("25A15", 25, "AM", 15, rep(NA, (ncol(qtc_info)-4)))

qtc_info <- filter(qtc_info, !(Subject == 30 & Timing == "AM"))

qtc_info <- rbind(qtc_info, a,b)
qtc_raw <- qtc_info ## keep a copy of unnormalized (by deviding baseline) qtc

norm_berry <- qtc_info %>% group_by(Timing,Subject) %>% group_split()

qtc_norm <- data.frame(matrix(ncol = ncol(qtc_raw), nrow = 0))
for (d in 1:length(norm_berry)){
  df.now <- norm_berry[[d]]
  bl.row <- df.now[which(df.now$Time == 0), ]
  
  df.norm <- df.now[1:4]
  for (m in 5:ncol(df.now)){
    mtbl.now <- colnames(df.now)[[m]]
    col.norm <- as.numeric(df.now[[m]])/as.numeric(bl.row[[m]])
    
    df.norm <- cbind(df.norm, col.norm)
    colnames(df.norm)[[ncol(df.norm)]] <- mtbl.now
    
    rm(m, mtbl.now, col.norm)
  }
  
  qtc_norm <- rbind(qtc_norm, df.norm)
  
  rm(d, df.now, bl.row, df.norm)
}
 
qtc_info <- qtc_norm ### pass that normalized qtc to qtc_info

# write.csv(qtc_info, file = paste0("otp/", "Baseline_devided_quant_auto01", ".csv"), row.names = F)
##### 2. Level-time line plots -- 
qtc_plot <- qtc_info

# #### add in lactate_b to pyruvate_1 ratio column
# qtc_plot$lactate_b_to_pyruvate_1_ratio <- qtc_plot$lactate_b/qtc_plot$pyruvate_1
# #### add in lactate_b to pyruvate_1 ratio column

# ### generate the ones that with all PM and MN devided by AM
# qtc_plot <- filter(qtc_plot, !(Subject == 30)) #the entire 30 can't be used, cause 30A are all NAs...
# qtc_plot$Time <- factor(qtc_plot$Time, levels = c("0", "15", "30", "45", "60", "90", "120", "150", "180"))
# qtc_plot <- qtc_plot[order(qtc_plot$Time), ]
# plot_berry <- qtc_plot %>% group_by(Subject, Timing) %>% group_split()
# 
# plot_now <- data.frame(matrix(ncol = ncol(qtc_plot), nrow = 0))
# for(q in c(1,4,7,10,13,16,19,22)){
#   am.df <- plot_berry[[q]]
#   mn.df <- plot_berry[[(q+1)]]
#   pm.df <- plot_berry[[(q+2)]]
#   n.sub <- unique(am.df$Subject)
#   
#   am.now <- cbind(am.df[1:4], am.df[-(1:4)]/am.df[-(1:4)])
#   mn.now <- cbind(mn.df[1:4], mn.df[-(1:4)]/am.df[-(1:4)])
#   pm.now <- cbind(pm.df[1:4], pm.df[-(1:4)]/am.df[-(1:4)])
# 
#   plot_now <- rbind(plot_now, am.now, mn.now, pm.now)
#   
#   rm(q, am.df, mn.df, pm.df, n.sub, am.now, mn.now, pm.now)
# }
# 
# qtc_plot <- plot_now
# qtc_plot$Time <- as.character(qtc_plot$Time)
# ### generate the ones that with all PM and MN devided by AM

qtc_berry <- qtc_plot %>% group_by(Timing, Time) %>% group_split()

library("DescTools")

lineplot_list <- list()
mtbl <- colnames(qtc_plot[-(1:4)])

for(m in mtbl){
  
  classi <- time_min <- "temp"
  median_plotter <- data.frame("median")
  lwr_plotter <- data.frame("lwrCI")
  upr_plotter <- data.frame("uprCI")
  
  for(i in 1:length(qtc_berry)){
    df <- qtc_berry[[i]]
    classi <- classi %>% append(unique(df$Timing))
    time_min <- time_min %>% append(unique(df$Time))
    
    median_list <- vector()
    median_list <- median_list %>% append(median(unlist(dplyr::select(df, all_of(m))), na.rm = T))
    
    lwr_list <- vector()
    # lwr_list <- lwr_list %>% append(MedianCI(unlist(dplyr::select(df, all_of(m))), na.rm = T)[[2]])
    lwr_list <- lwr_list %>% append(quantile(unlist(dplyr::select(df, all_of(m))), probs = 0.25, na.rm = T))
    
    upr_list <- vector()
    # upr_list <- upr_list %>% append(MedianCI(unlist(dplyr::select(df, all_of(m))), na.rm = T)[[3]])
    upr_list <- upr_list %>% append(quantile(unlist(dplyr::select(df, all_of(m))), probs = 0.75, na.rm = T))
    
    median_plotter <- median_plotter %>% rbind(median_list)
    
    lwr_plotter <- lwr_plotter %>% rbind(lwr_list)
    upr_plotter <- upr_plotter %>% rbind(upr_list)
    
    rm(median_list, lwr_list, upr_list)
  }
  
  mtbl_sinplotting <- data.frame(cbind(classi, time_min, median_plotter, lwr_plotter, upr_plotter))
  
  colnames(mtbl_sinplotting)<- c('Category', 'Time.in.Minutes', m, "Median_lwr", "Median_upr")
  
  mtbl_sinplotting <- mtbl_sinplotting[-1,]
  rownames(mtbl_sinplotting) <- NULL
  mtbl_sinplotting[4:ncol(mtbl_sinplotting)] <- lapply(mtbl_sinplotting[4:ncol(mtbl_sinplotting)],
                                                       as.numeric)
  
  mtbl_sinplotting <- within(mtbl_sinplotting, c(Time.in.Minutes <-  factor(Time.in.Minutes, levels = c("0", "15", "30", "45", "60", "90", "120", "150", "180")),
                                                 Category <-  as.factor(Category)))
  names(mtbl_sinplotting)[2] <- 'Time'
  
  lineplot_list <- append(lineplot_list, list(mtbl_sinplotting))
  print(m)
  
  rm(m, classi, time_min, median_plotter, lwr_plotter, upr_plotter, mtbl_sinplotting)
}
rm(mtbl)

##plotting
if(!file.exists("otp/BL_devided_AUC_Line_plots_all_TMs_auto01")){dir.create("otp/BL_devided_AUC_Line_plots_all_TMs_auto01")}
for(df in lineplot_list){
  df <- within(df, Category <- factor(Category, levels = c("AM", "PM", "MN"), ordered = T))
  df[3:ncol(df)] <- lapply(df[3:ncol(df)], as.numeric)
  mt_name = colnames(df)[[3]]
  print(mt_name)

  p <- plot_graph_line2(df, 'Category', 'Time', mtbl_name = mt_name, ymin = "Median_lwr", ymax = "Median_upr")
  ggsave(filename = paste0("otp/BL_devided_AUC_Line_plots_all_TMs_auto01/", mt_name, ".png")
         , plot = p, units = "in",  width = 11, height = 8, dpi = 300
         , device = "png")

  rm(df, mt_name, p)
}

# ####### 2.5. Last time point contrast by wilcoxon signed rank test
# qtc_last <- qtc_info[which(qtc_info$Time == "180"), ]
# ## 30 can only work in paired test for MN vs. PM, cause 30A is missing...
# last_berry <- qtc_last %>% group_by(Timing) %>% group_split()
# 
# m.names <- fc.pa <- fc.ma <- fc.mp <- p.pa <- p.ma <- p.mp <- vector()
#   
# for (l in 5:ncol(qtc_last)){
#   n.mtbl <- colnames(qtc_last)[[l]]
#   last.am <- last_berry[[1]][l]
#   last.pm <- last_berry[[3]][l]
#   last.mn <- last_berry[[2]][l]
#   
#   p.pva <- wilcox.test(as.numeric(last.pm[-6, ][[1]]), as.numeric(last.am[[1]]), paired = T)$p.value
#   p.mva <- wilcox.test(as.numeric(last.mn[-6, ][[1]]), as.numeric(last.am[[1]]), paired = T)$p.value
#   p.mvp <- wilcox.test(as.numeric(last.mn[[1]]), as.numeric(last.pm[[1]]), paired = T)$p.value
#   
#   fc.pva <- median(as.numeric(last.pm[-6, ][[1]]))/median(as.numeric(last.am[[1]]))
#   fc.mva <- median(as.numeric(last.mn[-6, ][[1]]))/median(as.numeric(last.am[[1]]))
#   fc.mvp <- median(as.numeric(last.mn[[1]]))/median(as.numeric(last.pm[[1]]))
#   
#   m.names <- append(m.names, n.mtbl)
#   fc.pa <- append(fc.pa, fc.pva)
#   fc.ma <- append(fc.ma, fc.mva)
#   fc.mp <- append(fc.mp, fc.mvp)
#   p.pa <- append(p.pa, p.pva)
#   p.ma <- append(p.ma, p.mva)
#  p.mp <- append(p.mp, p.mvp)
# 
#  rm(l, n.mtbl, last.am, last.pm, last.mn, p.pva, p.mva, p.mvp, fc.pva, fc.mva, fc.mvp)
# }
# 
# last.df <- as.data.frame(cbind(m.names, fc.pa, fc.ma, fc.mp, p.pa, p.ma, p.mp))
# colnames(last.df) <- c("analytes", "FC_PvA", "FC_MvA", "FC_MvP", "p_PvA","p_MvA", "p_MvP")
# last.df$`p.adj_PvA` <- p.adjust(as.numeric(last.df$`p_PvA`), method = "BH")
# last.df$`p.adj_MvA` <- p.adjust(as.numeric(last.df$`p_MvA`), method = "BH")
# last.df$`p.adj_MvP` <- p.adjust(as.numeric(last.df$`p_MvP`), method = "BH")
# 
# last.df$`log2FC_PvA` <- log(as.numeric(last.df$`FC_PvA`), 2)
# last.df$`log2FC_MvA` <- log(as.numeric(last.df$`FC_MvA`), 2)
# last.df$`log2FC_MvP` <- log(as.numeric(last.df$`FC_MvP`), 2)
# 
# write.csv(last.df, file = paste0("otp/", "Baseline_devided_180min_differences_auto01", ".csv"), row.names = F)

###### 3. AUC Calculation ##### - AUC output table for suppelementary...

##sort qtc_info based on Time point
#convert to numeric first
qtc_info <- qtc_info %>% mutate(across(4:ncol(qtc_info), as.numeric))
qtc_info <- arrange(qtc_info, Time)

### Start calculating AUC, iAUC/δAUC and inAUC, deAUC...
AUC_berry <- qtc_info %>% group_by(Subject, Timing) %>% group_split()

AUC_df <- qtc_info[1, -c(1,4)] ## Just a container.
AUC_df <- AUC_df %>% add_column('AUC_Type' = NA, .after = "Timing")

master_df_ref_list <- list()

for (i in 1:length(AUC_berry)){ #L1 ## For every Patient, every RCT...
  this_patient_rct <- AUC_berry[[i]]
  unique_id <- paste(unique(this_patient_rct$Subject), unique(this_patient_rct$Timing), sep = '_') 
  print(unique_id)
  
  to_add_inAUC <- to_add_deAUC <- to_add_iAUC <- to_add_AUC <- 
    c(unique(this_patient_rct$Subject), unique(this_patient_rct$Timing))
  
  to_add_AUC <- to_add_AUC %>% append('AUC')
  to_add_iAUC <- to_add_iAUC %>% append('iAUC')
  to_add_inAUC <- to_add_inAUC %>% append('inAUC')
  to_add_deAUC <- to_add_deAUC %>% append('deAUC')
  
  time_series <- this_patient_rct$Time
  
  instance_df_ref_list <- list()
  
  for (j in 5:ncol(this_patient_rct)){ #L2 ## For every metabolite...
    mtbl_name <- colnames(this_patient_rct[j])
    mtbl_counts <- unlist(this_patient_rct[,j]) 
    
    if(any(is.na(mtbl_counts))){ #L3-1 ## NA/Inf Catcher - is this even legit... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      mtbl_counts[is.na(mtbl_counts)] <- mean(mtbl_counts, na.rm = T) 
      ##could try geometric mean, median also can, but due to the negative-binomial distribution of the data, that won't 
      #actually wouldn't affect much -- Hercules
    } #3-1
    if(any(is.infinite(mtbl_counts))){ #L3-2
      mtbl_counts[which(is.infinite(mtbl_counts))] <- mean(mtbl_counts[is.finite(mtbl_counts)], na.rm = T)
    } #L3-2
    
    inAUC <- deAUC <- iAUC <- AUC <- vector()
    
    for (k in 1:(length(mtbl_counts)-1)){ #L3-3 ## For each timepoint...
      bl_count <- mtbl_counts[1]
      initial_count <- mtbl_counts[k]
      next_count <- mtbl_counts[k+1]
      time_diff <- time_series[k+1] - time_series[k]
      
      ## inAUC/deAUC/iAUC (incremental and 'decremental' and iAUC) === THINK triangles and trapeziums!
      # Here, we think for each two timepoints, there're always deAUC, inAUC and iAUC (which equals to deAUC + inAUC), even sometimes deAUC/inAUC
      # equals to 0!!! -- Jiangwen
      if((initial_count >= bl_count & next_count >= bl_count)|
         (initial_count <= bl_count & next_count <= bl_count)){ #L4-1 ## if both points in question are above/below BL values
        initial_diff <- initial_count - bl_count
        next_diff <- next_count - bl_count
        my_area <- (initial_diff + next_diff)/2 * time_diff
        
        iAUC <- iAUC %>% append(my_area)
        if (my_area < 0){ #L5-1
          inAUC <- inAUC %>% append(0)
          deAUC <- deAUC %>% append(my_area)
        }else{ #L5-2
          inAUC <- inAUC %>% append(my_area)
          deAUC <- deAUC %>% append(0)
        }
        
        rm(initial_diff, next_diff, my_area)
        
      }else if((initial_count >= bl_count & next_count < bl_count)|
               (initial_count < bl_count & next_count >= bl_count)){ #L4-2 ## if both points aren't all above/below BL values (one above one below)
        first_time_segment <- (initial_count - bl_count)/(initial_count - next_count) * time_diff 
        second_time_segment <- (next_count - bl_count)/(next_count - initial_count) * time_diff ## these time segments should both be +ve
        my_first_area <- first_time_segment * (initial_count - bl_count) / 2
        my_second_area <- second_time_segment * (next_count - bl_count) / 2
        
        ## iAUC
        iAUC <- iAUC %>% append(my_first_area + my_second_area)
        
        ## inAUC and deAUC
        if (my_first_area > 0){ #L5-3 # and 2nd area is < 0 (and 2nd area must be < 0)
          inAUC <- inAUC %>% append(my_first_area)
          deAUC <- deAUC %>% append(my_second_area)
        }else if (my_second_area > 0){ #L5-4 # and 1st area is < 0
          inAUC <- inAUC %>% append(my_second_area)
          deAUC <- deAUC %>% append(my_first_area)
        }
        
        rm(first_time_segment, second_time_segment, my_first_area, my_second_area)
      }
      
      ## AUC
      AUC <- AUC %>% append((initial_count + next_count)/2 * time_diff) ## Whole area under the 2 points
      
      rm(initial_count, bl_count, next_count, time_diff, k)
    } #L3-3
    
    df_area_reference <- data.frame(metabolite = rep(mtbl_name, 4), inAUC, deAUC, iAUC, AUC)
    instance_df_ref_list <- instance_df_ref_list %>% append(list(df_area_reference))
    
    to_add_inAUC[mtbl_name] <- sum(inAUC)
    to_add_deAUC[mtbl_name] <- sum(deAUC)
    to_add_iAUC[mtbl_name] <- sum(iAUC)
    to_add_AUC[mtbl_name] <- sum(AUC) #mtbl_name are the names of the segement of the vector to_add_(i/d/n)AUC
    
    rm(mtbl_name, mtbl_counts, j, inAUC, deAUC, iAUC, AUC, df_area_reference)
  } #L2
  
  AUC_df <- AUC_df %>% rbind(to_add_deAUC, to_add_inAUC, to_add_iAUC, to_add_AUC) 
  master_df_ref_list <- master_df_ref_list %>% append(list(unique_id, instance_df_ref_list)) # Should probably add a label...
  
  rm(instance_df_ref_list, to_add_deAUC, to_add_inAUC, to_add_iAUC, to_add_AUC,
     time_series, this_patient_rct, unique_id, i)
} #L1

AUC_df <- AUC_df[-1,] ## Get rid of that filler first row

rm(AUC_berry, master_df_ref_list)

### Extract data-frame for AUC and iAUC
realAUC_df <- filter(AUC_df, AUC_Type == "AUC")
iAUC_df <- filter(AUC_df, AUC_Type == "iAUC")

# write.csv(realAUC_df, file = paste0("otp/BL_devided_AUC_auto01", ".csv"), row.names = T)
# write.csv(iAUC_df, file = paste0("otp/BL_devided_iAUC_auto01", ".csv"), row.names = T)
## Convert all iAUC to positive values, reference: http://younglab.wi.mit.edu/chromatin/foldchange.html

info_col <- ct_iAUC_df <- data.frame(matrix(ncol = 0, nrow = 26))##ct_iAUC_df is corrected iAUC
for (m in 4:ncol(iAUC_df)){
  crt_data <- iAUC_df[c(1:3, m)]
  
  AM <- subset(crt_data, Timing == "AM")
  PM <- subset(crt_data, Timing == "PM")
  MN <- subset(crt_data, Timing == "MN")
  
  MinAM <- min(as.numeric(unlist(AM[4])))
  MinPM <- min(as.numeric(unlist(PM[4])))
  MinMN <- min(as.numeric(unlist(MN[4])))
  Min <- min(MinAM, MinPM, MinMN)
  
  AM[4] <- as.numeric(unlist(AM[4])) - Min + 5
  PM[4] <- as.numeric(unlist(PM[4])) - Min + 5
  MN[4] <- as.numeric(unlist(MN[4])) - Min + 5
  
  corrected_data <- rbind(AM, PM, MN)
  info_col <- cbind(info_col, corrected_data[1:3])
  ct_iAUC_df <- cbind(ct_iAUC_df, corrected_data[4])
  
  rm(m, crt_data, AM, PM, MN, MinAM, MinPM, MinMN, Min, corrected_data)
}
ct_iAUC_df <- cbind(info_col[1:3], ct_iAUC_df)
# write.csv(ct_iAUC_df, file = paste0("otp/BL_devided_ct_iAUC_auto01", ".csv"), row.names = T)

## remove 30AM's row to check if increase the p-values
realAUC_df_rm <- filter(realAUC_df, !(Subject == 30 & Timing == "AM"))
iAUC_df_rm <- filter(iAUC_df, !(Subject == 30 & Timing == "AM"))
ct_iAUC_df_rm <- filter(ct_iAUC_df, !(Subject == 30 & Timing == "AM"))
###### 4. Separation analysis based on AUC value - Figure
### PLS-DA (try both raw and log-transformation AUC)
### reference: https://stats.stackexchange.com/questions/179891/how-to-pre-process-data-for-partial-least-square-pls-regression-in-r
###http://mixomics.org/faq/
### Package set-up
library(mixOmics)

### Input data generation
## log2-transformed data
##log2 AUC
lg_realAUC_df <- realAUC_df
lg_realAUC_df[4:ncol(lg_realAUC_df)] <- lapply(lg_realAUC_df[4:ncol(lg_realAUC_df)], as.numeric)
lg_realAUC_df <- lg_realAUC_df %>% mutate(across(where(is.numeric), ~ log2(.)))
infinite_indices <- which(sapply(lg_realAUC_df, is.infinite), arr.ind = TRUE)
View(infinite_indices)
lg_realAUC_df_rm <- filter(lg_realAUC_df, !(Subject == 30 & Timing == "AM"))

##log2 corrected iAUC
lg_ct_iAUC_df <- ct_iAUC_df
lg_ct_iAUC_df[4:ncol(lg_ct_iAUC_df)] <- lapply(lg_ct_iAUC_df[4:ncol(lg_ct_iAUC_df)], as.numeric)
lg_ct_iAUC_df <-lg_ct_iAUC_df %>% mutate(across(where(is.numeric), ~ log2(.)))
infinite_indices <- which(sapply(lg_ct_iAUC_df, is.infinite), arr.ind = TRUE)
View(infinite_indices)
lg_ct_iAUC_df_rm <- filter(lg_ct_iAUC_df, !(Subject == 30 & Timing == "AM"))


## different types of input data
AUC_pls_df <- 
realAUC_df
# realAUC_df_rm

# iAUC_df
# iAUC_df_rm

# ct_iAUC_df
# ct_iAUC_df_rm
  
# lg_realAUC_df
# lg_realAUC_df_rm
  
# lg_ct_iAUC_df
# lg_ct_iAUC_df_rm
  

AUC_pca_df <- AUC_pls_df <- within(AUC_pls_df, Timing <- factor(Timing, levels = c("AM","PM" ,"MN"), ordered = T))

# #pca
# pca.auc <- prcomp(as.matrix(sapply(AUC_pca_df[, -c(1:3)], as.numeric)), center = TRUE, scale. = TRUE)
# pca.scores <- as.data.frame(pca.auc$x)
# pca.scores$Timing <- AUC_pca_df$Timing
# 
# pp <- ggplot(pca.scores, aes(x = PC1, y = PC2, color = Timing)) +
#   geom_point(size = 3) +
#   theme_minimal() +
#   labs(title = "PCA Plot",
#        x = paste0("PC1 (", round(100*summary(pca.auc)$importance[2,1], 1), "%)"),
#        y = paste0("PC2 (", round(100*summary(pca.auc)$importance[2,2], 1), "%)")) +
#   scale_color_manual(values = c("AM" = "red", "PM" = "blue", "MN" = "green"))
# pp #kind of meaningless...

#pls.da
pls.auc <- plsda(X = as.matrix(sapply(AUC_pls_df[, -c(1, 3)], as.numeric)), Y = AUC_pls_df$Timing, ncomp = 2, scale = T)

plotIndiv(pls.auc, ind.names = F, ellipse = TRUE, legend = TRUE) # this is just to check the comp1 and comp2
          
p <- plotIndiv(pls.auc, ind.names = F, ellipse = TRUE, legend = TRUE, 
               title = "", 
               X.label = "Comp1 (8%)", Y.label = "Comp2 (6%)",
               col = c("skyblue2", "dodgerblue4", "plum"),
               point.lwd = 0.75, cex = 1.8, style = "ggplot2")$graph + 
  theme_classic() +
  # ggtitle(paste0("PLS-DA of Metabolome")) + 
  theme(plot.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(size = 12, face = 'bold'),
        axis.title.y = element_text(size = 12, face = 'bold'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 12, face = 'bold'),
        legend.position = 'right',
        strip.background = element_blank())

ggsave(filename = "otp/PLS-DA_BL_devided_AUC_auto01.pdf"
       , plot = p, units = "in", width = 11, height = 8
       , device = "pdf")

rm(AUC_pls_df, pls.auc, p)

###### . Comparison using linear model and using anova to extract p-values
library(lme4)
library(lmerTest)
library(stats)
library(emmeans)

### prepare the linear model dataframe: AM, PM and MN as separate factor
lgAUC_lmm <- lg_realAUC_df
iAUC_lmm <- iAUC_df

lgAUC_lmm <- lgAUC_lmm %>% mutate(across(4:ncol(lgAUC_lmm), as.numeric))
iAUC_lmm <- iAUC_lmm %>% mutate(across(4:ncol(iAUC_lmm), as.numeric))

lgAUC_lmm$Timing <- factor(lgAUC_lmm$Timing, levels = c("AM", "PM", "MN")) ### luckily this won't automaticlly reorder the Timing column... otherwise all data are miss-matched...
iAUC_lmm$Timing <- factor(iAUC_lmm$Timing, levels = c("AM", "PM", "MN"))

lgAUC_lmm$Subject <- as.character(lgAUC_lmm$Subject)
iAUC_lmm$Subject <- as.character(iAUC_lmm$Subject)

# ###run single test for ratio of lactate_b/pyruvate_a
# single.df <- cbind(lgAUC_lmm[1:2], lgAUC_lmm[which(colnames(lgAUC_lmm) %in% c("lactate_b", "pyruvate_1"))])
# single.df$ratio <- single.df$lactate_b/single.df$pyruvate_1
# single.lm <- lmer(ratio ~ Timing + (1|Subject), single.df)
# single.emm <- emmeans(single.lm, ~Timing)
# summary(pairs(single.emm, adjust = "none"))$p.value
# summary(pairs(single.emm, adjust = "none"))$estimate
# summary(pairs(single.emm, adjust = "none"))$contrast
# ###run single test for ratio of lactate_b/pyruvate_a

### Fitting model
mtbl.name <- overall.p <- AvP.p <- AvM.p <- PvM.p <- AvP.ES <- AvM.ES <- PvM.ES <- vector()
re.tab <- data.frame(matrix(ncol = 5, nrow = 0))
ES.tab <- p.tab <- data.frame(matrix(ncol = 0, nrow = (ncol(lgAUC_lmm)-3))) #nrow = amount of analytes

###LMM.lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
for (m in 4:ncol(lgAUC_lmm)){
ctr.name <- colnames(lgAUC_lmm)[[m]]
ctr.df <- lgAUC_lmm[, c(1:2, m)]
### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# ### LMM.iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ####
# for (m in 4:ncol(iAUC_lmm)){
# ctr.name <- colnames(iAUC_lmm)[[m]]
# ctr.df <- iAUC_lmm[, c(1:2, m)]
# ### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^

  mtbl.name <- append(mtbl.name, ctr.name)
  otc.name <- paste0('`', ctr.name, '`')
  re.tab.col <- rep(ctr.name, 2)
  
  ctr.lm <- lmer(eval(parse(text = otc.name)) ~ Timing + (1|Subject), ctr.df)
  
  ctr.re <- cbind(as.data.frame(VarCorr(ctr.lm)), re.tab.col)
  re.tab <- rbind(re.tab, ctr.re)
  
  ## overall p-values across AM, PM and MN
  overall.p <- append(overall.p, anova(ctr.lm)[[6]])
  
  ## pairwise p-values and estimated coefficients (log2FC for lgAUC and effect size for iAUC) for three contrast: AM vs. PM, AM vs. MN, PM vs. MN
  emm_options(disable.transform = TRUE)
  ctr.emm <- suppressWarnings(emmeans(ctr.lm, ~Timing, tran = "log2")) #adding "supressWarnings" and "tran = "log2" "as suggested by Claude, just to supress the warnings.
  AvP.p <- append(AvP.p, summary(pairs(ctr.emm, adjust = "none"))$p.value[[1]])
  AvM.p <- append(AvM.p, summary(pairs(ctr.emm, adjust = "none"))$p.value[[2]])
  PvM.p <- append(PvM.p, summary(pairs(ctr.emm, adjust = "none"))$p.value[[3]])
  
  AvP.ES <- append(AvP.ES, summary(pairs(ctr.emm, adjust = "none"))$estimate[[1]])
  AvM.ES <- append(AvM.ES, summary(pairs(ctr.emm, adjust = "none"))$estimate[[2]])
  PvM.ES <- append(PvM.ES, summary(pairs(ctr.emm, adjust = "none"))$estimate[[3]])
  
  print(ctr.name)
  
  rm(m, ctr.name, otc.name, re.tab.col, ctr.df, ctr.lm, ctr.re, ctr.emm)
}
p.tab <- cbind(mtbl.name, overall.p, AvP.p, AvM.p, PvM.p)

# reverse log2FC direction for easier report
PvA.ES <- AvP.ES * -1
MvA.ES <- AvM.ES * -1
MvP.ES <- PvM.ES * -1

ES.tab <- cbind(mtbl.name, PvA.ES, MvA.ES, MvP.ES)

colnames(p.tab) <- c("biomarker", "p-value_overall", "p-value_AvP", "p-value_AvM", "p-value_PvM")
colnames(ES.tab) <- c("biomarker", "ES_PvA", "ES_MvA", "ES_MvP")
colnames(re.tab)<- c("group", "var1", "var2", "variance", "Std Dev", "biomarker")

###1. [full dataset]LMM output with indicating random effect
###lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^####
write.csv(re.tab, file = paste0("otp/BL_devided_lgAUC_LMM_Random_effect_auto01", ".csv"), row.names = T)

ori.df <- cbind(as.data.frame(ES.tab), as.data.frame(p.tab)[-1])
write.csv(ori.df, file = paste0("otp/BL_devided_lgAUC_[Ori]pval_ES_auto01", ".csv"), row.names = F)
#### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# ### LMM.iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^######
# write.csv(p.tab, file = paste0("otp/iAUC_[original]pval_auto11", ".csv"), row.names = F)
# write.csv(ES.tab, file = paste0("otp/iAUC_EffectSize_auto11", ".csv"), row.names = F)
# write.csv(re.tab, file = paste0("otp/iAUC_LMM_Random_effect_auto11", ".csv"), row.names = T)
# ### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

rm(mtbl.name, re.tab, overall.p, AvP.p, AvM.p, PvM.p, AvP.ES, AvM.ES, PvM.ES, PvA.ES, MvA.ES, MvP.ES)

####### 5. p-value adjustments
library(qvalue)
if(!file.exists("otp/AUCs_pval_hist")){dir.create("otp/AUCs_pval_hist")}
if(!file.exists("otp/AUCs_pval_adj")){dir.create("otp/AUCs_pval_adj")}

p.tab <- as.data.frame(p.tab)
BH_pval_df <- qvalue_df <- p.tab[1]

for(i in 2:ncol(p.tab)){
  
  ##checking original p-value hist curve
  name <- colnames(p.tab)[[i]]
  #### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ####
  pdf(file = paste0("otp/AUCs_pval_hist/BL_devided_lgAUC_Histogram_", name, '.pdf'), width = 11, height = 8)
  #### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  # #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ####
  # pdf(file = paste0("otp/AUCs_pval_hist/iAUC_Histogram_", name, '.pdf'), width = 11, height = 8)
  # #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # 
  hist(x = as.numeric(p.tab[[i]]), main = paste0("Histogram of ", name), xlab = "p-values")
  dev.off()
  
  ##doing BH and q-value method for p-value adjustment
  col_suf <- str_remove(name, "^p-value_")
  BH_pval <- p.adjust(p = as.numeric(p.tab[[i]]), method = "BH") # or use qvalue(.., pi0 = 1)$qvalues
  qvalue <- qvalue(p = as.numeric(p.tab[[i]]))$qvalues
  
  BH_pval_df <- cbind(BH_pval_df, BH_pval)
  qvalue_df <- cbind(qvalue_df, qvalue)
  colnames(BH_pval_df)[[ncol(BH_pval_df)]] <- paste0("BH_pval_", col_suf)
  colnames(qvalue_df)[[ncol(qvalue_df)]] <- paste0("qvalue_", col_suf)
  
  ##graphic comparing two p-value adjustment method
  #### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
  pdf(file = paste0("otp/AUCs_pval_adj/BL_devided_lgAUC_BH_vs_qvalue", "_", name, ".pdf"), width = 8, height = 8)
  #### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  # #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ######
  # pdf(file = paste0("otp/AUCs_pval_adj/iAUC_BH_vs_qvalue", "_", name, ".pdf"), width = 8, height = 8)
  # #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  plot(x = as.numeric(p.tab[[i]]), y = qvalue, xlab = "p-values", ylab = "q-values/BH adjusted p-values",
       main = paste0("BH vs q-value of ", name), cex = .1, col = 2, ylim = c(0,1), xlim = c(0,1))
 
  
  points(x = as.numeric(p.tab[[i]]), y = BH_pval, cex = .1, col = 4)
  legend("bottomright", c("q-value", "BH adjusted p-value"), pch = 19, col = c(2, 4))
  abline(a = 0, b = 1, col = 1)
  dev.off()
  
  rm(col_suf, BH_pval, qvalue, name)
}
#### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
write.csv(x = BH_pval_df, file = paste0("otp/BL_devided_lgAUC_BH_adj_pvals_auto01", ".csv"), row.names = F)
#### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
# write.csv(x = BH_pval_df, file = paste0("otp/iAUC_BH_adj_pvals_auto11", ".csv"), row.names = F)
# #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

report.df <- cbind(ES.tab, BH_pval_df[-1])

#### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
write.csv(report.df, file = paste0("otp/BL_devided_lgAUC_BHpval_ES_auto01", ".csv"), row.names = F)
#### lgAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ #####
# write.csv(report.df, file = paste0("otp/iAUC_BHpval_ES_auto11", ".csv"), row.names = F)
# #### iAUC ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


##### n. basline contrasts using LMM - for thesis #####
###### . Comparison using linear model and using anova to extract p-values
library(lme4)
library(lmerTest)
library(stats)
library(emmeans)

### prepare the linear model dataframe: AM, PM and MN as separate factor
lg_qtc_bl <- qtc_bl
lg_qtc_bl[1:5] <- lapply(lg_qtc_bl[1:5], as.character)
lg_qtc_bl[5:ncol(lg_qtc_bl)] <- lapply(lg_qtc_bl[5:ncol(lg_qtc_bl)], as.numeric)
lg_qtc_bl <- lg_qtc_bl %>% mutate(across(where(is.numeric), ~ log2(.)))

infinite_indices <- which(sapply(lg_qtc_bl, is.infinite), arr.ind = TRUE)
View(infinite_indices) # just checking is there any NAs

lg_qtc_bl$Timing <- factor(lg_qtc_bl$Timing, levels = c("AM", "PM", "MN")) ### luckily this won't automaticlly reorder the Timing column... otherwise all data are miss-matched...

### Fitting model
mtbl.name.bl <- overall.p.bl <- AvP.p.bl <- AvM.p.bl <- PvM.p.bl <- AvP.ES.bl <- AvM.ES.bl <- PvM.ES.bl <- vector()
re.tab.bl <- data.frame(matrix(ncol = 5, nrow = 0))
ES.tab.bl <- p.tab.bl <- data.frame(matrix(ncol = 0, nrow = (ncol(lg_qtc_bl)-4))) #nrow = amount of analytes

for (m in 5:ncol(lg_qtc_bl)){
  ctr.name <- colnames(lg_qtc_bl)[[m]]
  ctr.df <- lg_qtc_bl[, c(2:3, m)]
  
  mtbl.name.bl <- append(mtbl.name.bl, ctr.name)
  otc.name <- paste0('`', ctr.name, '`')
  re.tab.col <- rep(ctr.name, 2)
  
  ctr.lm <- lmer(eval(parse(text = otc.name)) ~ Timing + (1|Subject), ctr.df)
  
  ctr.re <- cbind(as.data.frame(VarCorr(ctr.lm)), re.tab.col)
  re.tab.bl <- rbind(re.tab.bl, ctr.re)
  
  ## overall p-values across AM, PM and MN
  overall.p.bl <- append(overall.p.bl, anova(ctr.lm)[[6]])
  
  ## pairwise p-values and estimated coefficients (log2FC) for three contrast: AM vs. PM, AM vs. MN, PM vs. MN
  emm_options(disable.transform = TRUE)
  ctr.emm <- suppressWarnings(emmeans(ctr.lm, ~Timing, tran = "log2")) #adding "supressWarnings" and "tran = "log2" "as suggested by Claude, just to supress the warnings.
  AvP.p.bl <- append(AvP.p.bl, summary(pairs(ctr.emm, adjust = "none"))$p.value[[1]])
  AvM.p.bl <- append(AvM.p.bl, summary(pairs(ctr.emm, adjust = "none"))$p.value[[2]])
  PvM.p.bl <- append(PvM.p.bl, summary(pairs(ctr.emm, adjust = "none"))$p.value[[3]])
  
  AvP.ES.bl <- append(AvP.ES.bl, summary(pairs(ctr.emm, adjust = "none"))$estimate[[1]])
  AvM.ES.bl <- append(AvM.ES.bl, summary(pairs(ctr.emm, adjust = "none"))$estimate[[2]])
  PvM.ES.bl <- append(PvM.ES.bl, summary(pairs(ctr.emm, adjust = "none"))$estimate[[3]])
  
  print(ctr.name)
  
  rm(m, ctr.name, otc.name, re.tab.col, ctr.df, ctr.lm, ctr.re, ctr.emm)
}

p.tab.bl <- cbind(mtbl.name.bl, overall.p.bl, AvP.p.bl, AvM.p.bl, PvM.p.bl)

# reverse log2FC direction for easier report
PvA.ES.bl <- AvP.ES.bl * -1
MvA.ES.bl <- AvM.ES.bl * -1
MvP.ES.bl <- PvM.ES.bl * -1

ES.tab.bl <- cbind(mtbl.name.bl, PvA.ES.bl, MvA.ES.bl, MvP.ES.bl)

colnames(p.tab.bl) <- c("biomarker", "p-value_overall", "p-value_AvP", "p-value_AvM", "p-value_PvM")
colnames(ES.tab.bl) <- c("biomarker", "ES_PvA", "ES_MvA", "ES_MvP")
colnames(re.tab.bl)<- c("group", "var1", "var2", "variance", "Std Dev", "biomarker")


write.csv(re.tab.bl, file = paste0("otp/Baseline_LMM_Random_effect_auto01", ".csv"), row.names = T)

ori.df.bl <- cbind(as.data.frame(ES.tab.bl), as.data.frame(p.tab.bl)[-1])
write.csv(ori.df.bl, file = paste0("otp/Baseline_[Ori]pval_ES_auto01", ".csv"), row.names = F)

rm(mtbl.name.bl, re.tab.bl, overall.p.bl, AvP.p.bl, AvM.p.bl, PvM.p.bl, AvP.ES.bl, AvM.ES.bl, PvM.ES.bl, PvA.ES.bl, MvA.ES.bl, MvP.ES.bl)

#######p-value adjustments
library(qvalue)
if(!file.exists("otp/Baseline_pval_hist")){dir.create("otp/Baseline_pval_hist")}
if(!file.exists("otp/Baseline_pval_adj")){dir.create("otp/Baseline_pval_adj")}

p.tab.bl <- as.data.frame(p.tab.bl)
BH_pval_df.bl <- qvalue_df.bl <- p.tab.bl[1]

for(i in 2:ncol(p.tab.bl)){
  
  ##checking original p-value hist curve
  name <- colnames(p.tab.bl)[[i]]
  pdf(file = paste0("otp/Baseline_pval_hist/Baseline_Histogram_", name, '.pdf'), width = 11, height = 8)
  
  hist(x = as.numeric(p.tab.bl[[i]]), main = paste0("Histogram of ", name), xlab = "p-values")
  dev.off()
  
  ##doing BH and q-value method for p-value adjustment
  col_suf <- str_remove(name, "^p-value_")
  BH_pval <- p.adjust(p = as.numeric(p.tab.bl[[i]]), method = "BH") # or use qvalue(.., pi0 = 1)$qvalues
  qvalue <- qvalue(p = as.numeric(p.tab.bl[[i]]))$qvalues
  
  BH_pval_df.bl <- cbind(BH_pval_df.bl, BH_pval)
  qvalue_df.bl <- cbind(qvalue_df.bl, qvalue)
  colnames(BH_pval_df.bl)[[ncol(BH_pval_df.bl)]] <- paste0("BH_pval_", col_suf)
  colnames(qvalue_df.bl)[[ncol(qvalue_df.bl)]] <- paste0("qvalue_", col_suf)
  
  ##graphic comparing two p-value adjustment method
  pdf(file = paste0("otp/Baseline_pval_adj/Baseline_BH_vs_qvalue", "_", name, ".pdf"), width = 8, height = 8)
  
  plot(x = as.numeric(p.tab.bl[[i]]), y = qvalue, xlab = "p-values", ylab = "q-values/BH adjusted p-values",
       main = paste0("BH vs q-value of ", name), cex = .1, col = 2, ylim = c(0,1), xlim = c(0,1))
  
  
  points(x = as.numeric(p.tab.bl[[i]]), y = BH_pval, cex = .1, col = 4)
  legend("bottomright", c("q-value", "BH adjusted p-value"), pch = 19, col = c(2, 4))
  abline(a = 0, b = 1, col = 1)
  dev.off()
  
  rm(col_suf, BH_pval, qvalue, name)
}

write.csv(x = BH_pval_df.bl, file = paste0("otp/Baseline_BH_adj_pvals_auto01", ".csv"), row.names = F)

report.df.bl <- cbind(ES.tab.bl, BH_pval_df.bl[-1])

write.csv(report.df.bl, file = paste0("otp/Baseline_BHpval_ES_auto01", ".csv"), row.names = F)

