# LIBRARIES ====
library(plyr)
library(aplpack)

# PREPROCESSING ====
source('./preprocessing.R')

# STARS ====
stars(head(raw.filtered))

# CHERNOFF FACES ====
faces(head(raw.filtered, n = 100))

# SOURCE KMEANS (TAKES A WHILE!) ====
source('./kmeans-dtree.R')

# CHERNOFF KMEANS ====
faces(head(clusters.raw[["4"]], n = 100))

# Plot by ascending age!
faces(head(arrange(clusters.raw[["4"]], age), n = 100))
stars(head(arrange(clusters.raw[["4"]], age), n = 100))
summaries[["4"]]
# 

# PLOT KMEANS CLUSTER CENTERS ====
centers4 <- trees$clusters4$clustering$centers
par(mfrow=c(1, 2))
plot(faces(centers4))
stars(centers4, col.stars = 4:8)
par(mfrow=c(1, 1))