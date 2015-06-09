# ==== GLOBAL LIBRARIES ====
library(plyr)

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