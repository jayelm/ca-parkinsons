# GLOBAL LIBRARIES ====
library(plyr)
# For PCA
library(FactoMineR)
library(psych)

# GLOBAL CONSTANTS ====
SAVE.PREPROCESSING.PLOTS <- FALSE
DB_FILE <- "../data/DATABASE_NMS Burden levels_15-4-2012.csv"
MCP_FILE <- "../data/Database 951_Subtyping_30-5-2016.csv"
ALL.SYMPTOMS <- c("nms_d1",
                 "nms_d2",
                 "nms_d3",
                 "nms_d4",
                 "nms_d5",
                 "nms_d6",
                 "nms_d7",
                 "nms_d8",
                 "nms_d9",
                 "tremor",
                 "bradykin",
                 "rigidity",
                 "axial",
                 "scmmotcp")
MOTOR.SYMPTOMS <- c("tremor",
                   "bradykin",
                   "rigidity",
                   "axial",
                   "scmmotcp")
MOTOR.PUB <- c("Tremor",
                   "Bradykinesia",
                   "Rigidity",
                   "Axial",
                  "Motor_comp")
NMS.30.NAMES <- c(
  'd1-1-lightheaded',
  'd1-2-fainting',
  'd2-3-drowsiness',
  'd2-4-fatigue',
  'd2-5-insomnia',
  'd2-6-rls',
  'd3-7-loss_interest',
  'd3-8-loss_activities',
  'd3-9-anxiety',
  'd3-10-depression',
  'd3-11-flat_affect',
  'd3-12-loss_pleasure',
  'd4-13-hallucination',
  'd4-14-delusion',
  'd4-15-diplopia',
  'd5-16-loss_concentration',
  'd5-17-forget_explicit',
  'd5-18-forget_implicit',
  'd6-19-drooling',
  'd6-20-swallowing',
  'd6-21-constipation',
  'd7-22-urinary_urgency',
  'd7-23-urinary_frequency',
  'd7-24-nocturia',
  'd8-25-sex_drive',
  'd8-26-sex_dysfunction',
  'd9-27-unexplained_pain',
  'd9-28-gust_olfact',
  'd9-29-weight_change',
  'd9-30-sweating'
)
NMS.30.NAMES.PUB <- c(
  'Lightheadedness',
  'Fainting',
  'Drowsiness',
  'Fatigue',
  'Insomnia',
  'RLS',
  'Loss_interest',
  'Loss_activities',
  'Anxiety',
  'Depression',
  'Flat_affect',
  'Loss_pleasure',
  'Hallucination',
  'Delusion',
  'Diplopia',
  'Loss_concentration',
  'Forget_explicit',
  'Forget_implicit',
  'Drooling',
  'Swallowing',
  'Constipation',
  'Urinary_urgency',
  'Urinary_frequency',
  'Nocturia',
  'Sex_drive',
  'Sex_dysfunction',
  'Unexplained_pain',
  'Gustation_olfaction',
  'Weight_change',
  'Sweating'
)
NMS.30 <- c(
  "nms1", "nms2", "nms3", "nms4", "nms5", "nms6",
  "nms7", "nms8", "nms9", "nms10", "nms11", "nms12",
  "nms13", "nms14", "nms15", "nms16", "nms17", "nms18",
  "nms19", "nms20", "nms21", "nms22", "nms23", "nms24",
  "nms25", "nms26", "nms27", "nms28", "nms29", "nms30"
)
NMS.D <- c(
  "nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5", "nms_d6", "nms_d7", "nms_d8", "nms_d9"
)
NMS.D.NAMES <- c(
  "Cardiovascular", "Sleep/fatigue", "Mood/apathy", "Perception/hallucination", "Attention/memory",
  "Gastrointestinal", "Urinary", "Sexual", "Miscellaneous"
)
NMS.NUM.TO.PUB <- setNames(sapply(NMS.30.NAMES, function(s) paste(strsplit(s, "-")[[1]][c(1, 3)], collapse = "-")),
                           NMS.30)
NMS.D.TO.NUM <- setNames(sapply(1:length(NMS.D), function(i) paste(i, "-", NMS.D[i], sep = "")), NMS.D)
NMS.D.MAP.PUB <- setNames(NMS.D.NAMES, NMS.D)
NMS.D.MAP.PUB.N <- setNames(sapply(1:length(NMS.D),
                                   function(i) paste(i, NMS.D.NAMES[i], sep = "-")), NMS.D)
NMS.30.MAP.PUB <- setNames(NMS.30.NAMES.PUB, NMS.30)
NMS.30.LONG.SHORT.MAP <- setNames(NMS.30.NAMES.PUB, NMS.30.NAMES)
MISC.MAP <- c(
  'age' = 'Age',
  'sex' = 'Sex',
  'pdonset' = 'PD_onset',
  'durat_pd' = 'PD_duration',
  'dyskinesia' = 'Dyskinesia',
  'fluctuat' = 'Fluctuations',
  'scmmotcp' = 'Motor_comp',
  'cisitot' = 'CISI_PD_total',
  'tremor' = 'Tremor',
  'bradykin' = 'Bradykinesia',
  'rigidity' = 'Rigidity',
  'axial' = 'Axial',
  'ldopa' = 'L-Dopa',
  'surgery' = 'Surgery',
  'cluster' = 'Cluster'
)
PUB.MAP <- c(NMS.D.MAP.PUB, NMS.30.MAP.PUB, MISC.MAP)
PUB.MAP.N <- c(NMS.D.MAP.PUB.N, NMS.30.MAP.PUB, MISC.MAP)  # This one just has n
ALL.BUT.NMS <- c('age', 'sex', 'pdonset', 'durat_pd', 'cisitot',
                 'tremor', 'bradykin', 'rigidity', 'axial', 'scmmotcp')
INTERPRETED <- c("age", "sex", "pdonset", "durat_pd", "cisitot",
                 "nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5",
                 "nms_d6", "nms_d7", "nms_d8", "nms_d9",
                 "tremor", "bradykin", "rigidity", "axial", 'scmmotcp',
                 NMS.30)
# Change this variable to change symptoms
SYMPTOMS.TO.USE <- ALL.SYMPTOMS

# IMPORT DATA, PREPROCESSING ====
raw <- read.csv(DB_FILE)
mcp <- read.csv(MCP_FILE)
raw$dyskinesia <- mcp$dyskinesia
raw$fluctuat <- mcp$fluctuat
# Study is irrelevant
# raw$study <- NULL
# Get rid of NAs (later may use some kind of missing value compensation method
# This only brings it down to 904, rather than ~700 something!
raw.omitted <- raw[, INTERPRETED]
raw.omitted <- na.omit(raw.omitted)
raw.omitted.full <- cbind(
  raw.omitted,
  # Now get the rownames from the raw
  raw[rownames(raw.omitted), -which(names(raw) %in% INTERPRETED)]
)
# To get identity of excluded patients: study/nident/country
raw[!rownames(raw) %in% rownames(raw.omitted), c("study", "nident", "country")]
# Verify they all have NAs in the things we care about
all(apply(raw[!rownames(raw) %in% rownames(raw.omitted), INTERPRETED], MARGIN = 1, anyNA))

# raw with only select nms and motor symptoms
raw.filtered <- raw.omitted[, SYMPTOMS.TO.USE]

# DESCRIPTIVE STATISTICS (BEFORE STANDARDIZATION) ====
# TODO - change raw.filtered to raw.omitted.filtered, or something like that
raw.omitted.stats <- describe(raw.omitted)
raw.filtered.stats <- describe(raw.filtered)

# STANDARDIZATION ====
raw.filtered.unscaled <- raw.filtered
raw.filtered <- as.data.frame(scale(raw.filtered))

# PCA (Note: not helpful for identifying specific factors) ====
par(mfrow = c(1, 2))
pca <- PCA(raw.filtered)
if (SAVE.PREPROCESSING.PLOTS) {
  dev.copy(pdf, "../figures/pca.pdf")
  dev.off()
}
# pca$eig
# pca$var$coord
# head(pca$ind$coord)
par(mfrow = c(1, 1))
plot(x = 1:length(rownames(pca$eig)), y = (pca$eig$eigenvalue), type="b",
     xlab="Factor", ylab="Eigenvalue", xaxp=c(0, 19, 19))
if (SAVE.PREPROCESSING.PLOTS) {
  dev.copy(pdf, "../figures/pca-eigenvalues.pdf")
  dev.off()
}

# Elbow appears around 3 factors
# Reset par for the rest of the script

