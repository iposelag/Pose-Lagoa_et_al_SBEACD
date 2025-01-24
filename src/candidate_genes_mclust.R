#########################################################################################
# Selection of candidate genes mclust method
#########################################################################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load required libraries
library(pheatmap)
library(mclust)
load("../COPD/raw_data/expression.Rda")
load("../COPD/raw_data/phenotypic.Rda")

## ----------------------------------------------------------------------------------------------------------------------------------------
# Load data
norm <- read.table("../COPD/norm_results.csv", header=T, sep=",", stringsAsFactors = F)
exp <- t(combat_copd_ctrl[, 2:ncol(combat_copd_ctrl)])
sel <- c("dis_condition", "GOLD_stage", "sex", "age", "platform_id")
meta <- phenotypic_ctics[colnames(exp), sel]
meta$age <- as.numeric(meta$age)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Obtain clusters using the model distribution of aggregated scores
X <- norm$aggregated_max
BIC <- mclustBIC(X)
plot(BIC)
mod1 <- Mclust(X, x = BIC, G=1:9)
summary(mod1, parameters = TRUE)
table(mod1$classification)
ICL <- mclustICL(X)
summary(ICL)

## ----------------------------------------------------------------------------------------------------------------------------------------
# Define threshold

# Plot selected clusters and treshold
sel_groups <- c(8,9)
pdf("../COPD/plots/cluster_aggregated_SHAP.pdf")
par(mfrow=c(2,1))
plot(mod1, what = "classification", xlim=c(0, max(X)))
hist(norm$aggregated_max, 100, xlim=c(0, max(X)))
abline(v=min(X[mod1$classification %in% sel_groups]), col="red")
dev.off()
summary(X[mod1$classification==sel_groups])

# Generate candidate list
th <- min(X[mod1$classification %in% sel_groups])
cand <- norm$variable_name[norm$aggregated_max>=th]
length(cand)
write.table(cand, file="candidate_ml_clust.tsv", sep="\t", row.names=F, col.names=F, quote=F)

# cand_exp <- exp[cand, ]
# pheatmap(cand_exp, annotation_col = meta, scale="none")
# pheatmap(cor(cand_exp, method = "pearson"), annotation_col = meta, scale="none")


