# Parkinson's disease subtypes identified from cluster analysis of motor and non-motor symptoms

This repository contains data, analyses, figures, and manuscripts for

> Jesse Mu, Kallol Ray Chaudhuri, Concha Bielza, Jesús de Pedro Cuesta, Pedro Larrañaga, and Pablo Martinez-Martin. Parkinson's disease subtypes identified from cluster analysis of motor and non-motor symptoms To appear in *Frontiers in Aging Neuroscience*

The code to reproduce figures is mostly in `publication.R`
which depends first on `source`ing `preprocessing.R` and (maybe??)
`kmeans-dtree.R` and `nms30.R`.

## TODO

- Organize and figure out procedure for exactly reproducing figures. Might be
 able to put everything in preprocessing
- Rename files to something more sensible (kmeans-dtree -> domains, nms30 -> symptoms)
- Continue to get rid of unnecessary figures

![PD analysis](./figures/png/analysis.png)
![Hierarchical clustering](./figures/png/hc.png)
![Longitudinal analysis](./figures/png/long.png)
