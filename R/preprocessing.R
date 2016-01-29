# GLOBAL LIBRARIES ====
library(plyr)
# For PCA
library(FactoMineR)
library(psych)

# GLOBAL CONSTANTS ====
SAVE.PREPROCESSING.PLOTS <- FALSE
DB_FILE <- "../data/DATABASE_NMS Burden levels_15-4-2012.csv"
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
                 "axial")
MOTOR.SYMPTOMS <- c("tremor",
                   "bradykin",
                   "rigidity",
                   "axial")
NMS.30 <- c(
  "nms1", "nms2", "nms3", "nms4", "nms5", "nms6",
  "nms7", "nms8", "nms9", "nms10", "nms11", "nms12",
  "nms13", "nms14", "nms15", "nms16", "nms17", "nms18",
  "nms19", "nms20", "nms21", "nms22", "nms23", "nms24",
  "nms25", "nms26", "nms27", "nms28", "nms29", "nms30"
)
INTERPRETED <- c("age", "sex", "pdonset", "durat_pd", "cisitot",
                 "nms_d1", "nms_d2", "nms_d3", "nms_d4", "nms_d5",
                 "nms_d6", "nms_d7", "nms_d8", "nms_d9",
                 "tremor", "bradykin", "rigidity", "axial",
                 NMS.30)
# Change this variable to change symptoms
SYMPTOMS.TO.USE <- ALL.SYMPTOMS

# UTILITY FUNCTIONS ====
splitdf <- function(dataframe, seed=NULL, trainfrac=0.7) {
  if (trainfrac<=0 | trainfrac>=1) stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

# IMPORT DATA, PREPROCESSING ====
raw <- read.csv(DB_FILE)
# Study is irrelevant
raw$study <- NULL
# Get rid of NAs (later may use some kind of missing value compensation method
# This only brings it down to 904, rather than ~700 something!
raw.omitted <- raw[, INTERPRETED]
raw.omitted <- na.omit(raw.omitted)

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
