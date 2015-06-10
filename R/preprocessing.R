# ==== GLOBAL LIBRARIES ====
library(plyr)
# For PCA
library(FactoMineR)
library(psych)

# ==== GLOBAL CONSTANTS ====
DB_FILE = "../data/DATABASE_NMS Burden levels_15-4-2012.csv"
ALL.SYMPTOMS = c("nms_d1",
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
                 "pigd")
MOTOR.SYMPTOMS = c("tremor",
                   "bradykin",
                   "rigidity",
                   "axial",
                   "pigd")
# Change this variable to change symptoms
SYMPTOMS.TO.USE = ALL.SYMPTOMS

# ==== UTILITY FUNCTIONS ====
splitdf <- function(dataframe, seed=NULL, trainfrac=0.7) {
  if (trainfrac<=0 | trainfrac>=1) stop("Training fraction must be between 0 and 1, not inclusive")
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/(1/trainfrac)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

# ==== IMPORT DATA, PREPROCESSING ====
raw <- read.csv(DB_FILE)
# Study is irrelevant
raw$study <- NULL
# Get rid of NAs (later may use some kind of missing value compensation method
raw.omitted <- na.omit(raw)

# raw with only select nms and motor symptoms
raw.filtered <- raw.omitted[, SYMPTOMS.TO.USE]

# ==== DESCRIPTIVE STATISTICS (BEFORE STANDARDIZATION) ====
# TODO - change raw.filtered to raw.omitted.filtered, or something like that
raw.omitted.stats <- describe(raw.omitted)
raw.filtered.stats <- describe(raw.filtered)

# ==== STANDARDIZATION ====
raw.filtered <- as.data.frame(scale(raw.filtered))

# ==== PCA (Note: not helpful for identifying specific factors) ====
par(mfrow = c(1, 2))
pca <- PCA(raw.filtered)
# pca$eig
# pca$var$coord
# head(pca$ind$coord)
plot(x = 1:length(rownames(pca$eig)), y = (pca$eig$eigenvalue), type="b",
     xlab="Factor", ylab="Eigenvalue", xaxp=c(0, 14, 14))
# Elbow appears around 3 factors
# Reset par for the rest of the script
par(mfrow = c(1, 1))
