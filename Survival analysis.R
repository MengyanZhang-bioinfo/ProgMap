
rm(list = ls()) 
library(survival)
library(survminer)

dt <- read.csv("./os.csv", header = TRUE)
dt$OS.time <- dt$OS.time / 30.42

dt$cluster <- factor(dt$cluster, levels = c(1, 2), labels = c("ISPC", "ACHC"))

dt$stage <- trimws(as.character(dt$stage))
dt$Stage_Group <- NA
dt$Stage_Group[grepl("I|1", dt$stage, ignore.case = TRUE)] <- "Stage I"
dt$Stage_Group[grepl("II|2|III|3", dt$stage, ignore.case = TRUE)] <- "Stage II/III"
dt <- dt[!is.na(dt$Stage_Group), ]
dt$Four_Groups <- paste(dt$cluster, dt$Stage_Group, sep = " + ")

target_levels <- c(
  "ISPC + Stage I",       
  "ISPC + Stage II/III",  
  "ACHC + Stage I",       
  "ACHC + Stage II/III"   
)

target_colors <- c(
  "#C0C0C0",  
  "#4F4F4F",  
  "#7FFFD4",  
  "#008B8B"   
)

dt$Four_Groups <- factor(dt$Four_Groups, levels = target_levels)
dt$Four_Groups <- droplevels(dt$Four_Groups)
final_levels <- levels(dt$Four_Groups) 
match_indices <- match(final_levels, target_levels) 
final_palette <- target_colors[match_indices] 


fit <- survfit(Surv(OS.time, OS) ~ Four_Groups, data = dt)

ggsurvplot(
  fit,                     
  data = dt,
  palette = final_palette, 
  linetype = "solid",   
  size = 1,                
  pval = TRUE,             
  pval.method = TRUE,
  conf.int = FALSE,        
  
  risk.table = TRUE,       
  risk.table.col = "strata",
  risk.table.height = 0.30,
  xlab = "Time (Months)",  
  ylab = "Overall Survival Probability",
  legend.title = "Subgroups",

  legend.labs = final_levels,
  ggtheme = theme_classic()
)


