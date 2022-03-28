## LOAD LIBRARIES
library(dplyr)
library(tidyverse)
library(mlr)
library(gridExtra)
library(ggfortify)
library(edgeR)
library(org.Hs.eg.db)
library(RColorBrewer)
library(gplots)
library(tidyr)
library(pheatmap)
library(reshape2)
library(survminer)
library(rms)
library(pheatmap)
library(ggplot2)

## SET THEME FOR SOME PLOTS
my_theme <- function(base_size =9, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text = element_text(size=8),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      axis.title = element_text(size =8),
      panel.grid.major = element_line(color = "gray"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "#ffffff"),
      strip.background = element_rect(fill = "black", color = "black", size =0.5),
      strip.text = element_text(face = "bold", size = 8, color = "white"),
      legend.position = "bottom",
      legend.justification = "center",
      legend.background = element_blank(),
      panel.border = element_rect(color = "grey5", fill = NA, size = 0.5)
    )
}

theme_set(my_theme())

mycolors=c("darkred","black","#004431","#000c44","#2d1600","purple","#202302")
myfillcolors=c("#d10c0c","grey10","#1e6651","#213989","#5b431a","purple","#393f03")

##########################
## CLINICAL DATA
##########################

clin_data <- read.csv("sc3_Training_ClinAnnotations.csv")

## group age into bins for visualization purposes
labs <- c(paste(seq(20, 80, by=15), seq(20+15-1, 95-1, by=15), sep="-"))
clin_data$agegrp <- cut(clin_data$D_Age, breaks=c(seq(20,80, by=15), Inf), labels=labs, right=F)


## Pre-processing
table(clin_data$HR_FLAG)

## assumption: events where HR_FLAG it is not TRUE (event did not occur) are assumed to be censored
## 131 out of 583 events occurred and rest are censored
clin_data$HR_FLAG <- clin_data$HR_FLAG==TRUE

## function to check if any NAs present in the data
nacount <- apply(clin_data, 2, function(x){
  sum(is.na(x))
})
nacount <- as.data.frame(nacount, row.names=names(clin_data))
nacount <- add_rownames(nacount) %>% dplyr::mutate(percentNA=round(nacount/583*100,2))%>%
  dplyr::filter(percentNA>=0.00)

## remove columns where NAs >50%
cols_to_remove <- nacount[nacount$percentNA>50,]$rowname
clin_data <- clin_data[, !(names(clin_data) %in% cols_to_remove)]

## remove unnecessary columns -> keep only D_PFS as time column,  remove file info columns, 
## study and patient status columns also removed (all same)
clin_data <- clin_data[, -c(1, 5:6, 8, 10:11,13:17)]

## remove NAs (only 20 in D_ISS column)
clin_data <- clin_data[complete.cases(clin_data),]

## Clinical data exploration

## age by gender distribution plot
age_sex <- clin_data %>% ggplot(aes(x=agegrp, color=D_Gender, fill=D_Gender))+geom_bar()+
  ggtitle("Age distribution by gender")+xlab("Age group [years]")+ylab("Number of subjects")+theme_minimal()

## Kaplan-Meier curve for events in patients stratified by age and ISS
surv_plots <- list()
fit_iss=survival::survfit(survival::Surv(D_PFS,HR_FLAG)~D_ISS,data=clin_data)
iss_plot <- ggsurvplot(fit_iss, pval = TRUE, conf.int = TRUE)

#autoplot(fit_iss)+scale_color_manual(values=mycolors)+scale_fill_manual(values=myfillcolors)

fit_age = survival::survfit(survival::Surv(D_PFS,HR_FLAG)~agegrp,data=clin_data)
age_plot <- ggsurvplot(fit_age, pval = TRUE, conf.int = TRUE)
#autoplot(fit_age)+scale_color_manual(values=mycolors)+scale_fill_manual(values=myfillcolors)

## simple 2-predictor model with coxph using only age and ISS 
age_iss_model <- survival::coxph(survival::Surv(D_PFS,HR_FLAG) ~ 
                          D_Age + D_ISS , data=clin_data)


summary(age_iss_model)

## both age and ISS have significant effect, so keep both features for further analysis

## simple 13-predictor model using only cyto features with coxph - to detect which cyto features could be relevant as significant predictors
cyto_feat_model <- survival::coxph(survival::Surv(D_PFS,HR_FLAG) ~ 
                                     CYTO_predicted_feature_01+CYTO_predicted_feature_02+CYTO_predicted_feature_03+
                                     CYTO_predicted_feature_05+CYTO_predicted_feature_06+CYTO_predicted_feature_08+
                                     CYTO_predicted_feature_12+CYTO_predicted_feature_13+CYTO_predicted_feature_14+
                                     CYTO_predicted_feature_15+CYTO_predicted_feature_16+CYTO_predicted_feature_17+
                                     CYTO_predicted_feature_18, data=clin_data)
summary(cyto_feat_model)

## only feature 01, 03 and 14 have a significant effect (based on p-values) - plot Kaplan meier curves for these features
fit_cyto1 = survival::survfit(survival::Surv(D_PFS,HR_FLAG)~CYTO_predicted_feature_01,data=clin_data)
fit_cyto3 = survival::survfit(survival::Surv(D_PFS,HR_FLAG)~CYTO_predicted_feature_03,data=clin_data)
fit_cyto14 = survival::survfit(survival::Surv(D_PFS,HR_FLAG)~CYTO_predicted_feature_14,data=clin_data)
cyto1_plot <- ggsurvplot(fit_cyto1,
           pval = TRUE, conf.int = TRUE)
cyto3_plot <- ggsurvplot(fit_cyto3,
           pval = TRUE, conf.int = TRUE)
cyto14_plot <- ggsurvplot(fit_cyto14,
           pval = TRUE, conf.int = TRUE)

## include only these 3 cyto features as they seem to have significant effect on survival 
clin_data <- clin_data %>% dplyr::select(-c("CYTO_predicted_feature_02", 
                                            "CYTO_predicted_feature_05", "CYTO_predicted_feature_06", "CYTO_predicted_feature_08", 
                                            "CYTO_predicted_feature_12", "CYTO_predicted_feature_13",  
                                            "CYTO_predicted_feature_15", "CYTO_predicted_feature_16", "CYTO_predicted_feature_17", 
                                            "CYTO_predicted_feature_18"))

## Descriptive statistics

psych::describeBy(clin_data$D_PFS,group=clin_data$CYTO_predicted_feature_01)
psych::describeBy(clin_data$D_PFS,group=clin_data$CYTO_predicted_feature_03)
psych::describeBy(clin_data$D_PFS,group=clin_data$CYTO_predicted_feature_14)
psych::describeBy(clin_data$D_PFS,group=clin_data$agegrp)
psych::describeBy(clin_data$D_PFS,group=clin_data$D_ISS)

#arrange_ggsurvplots(surv_plots, print = TRUE, ncol = 2, nrow = 1)

############################
## EXPRESSION DATA
############################

exp_data <- read.csv("MMRF_CoMMpass_IA9_E74GTF_Salmon_entrezID_TPM_hg19.csv")

## list the risk genes - based on literature research

riskgenes <- c("KRAS", "NRAS", "DIS3", "TENT5C",
               "BRAF","HUWE1", "TP53", "TRAF3",
               "EGR1", "ATM", "H1-4", "FGFR3", "UBR5", "PRKD2", "CYLD", "ACTG1", 
               "IRF4", "MAX", "KMT2C", "CREBBP", "CCND1", "ARID1A", "EPAS1", 
               "ERC2", "PRC1", "CSGALNACT1", "PHF19", "NSD2")

# Remove first column from expdata - countdata, contains only the counts for all samples.
countdata_initial <- exp_data[,-1]

# Store EntrezGeneID as rownames
rownames(countdata_initial) <- exp_data[,1]

# make dgelist object
dge_initial <- DGEList(countdata_initial)

## get annotations for entrez ids
ann <- select(org.Hs.eg.db,keys=rownames(dge_initial$counts),columns=c("ENTREZID","SYMBOL","GENENAME"))

##double check that the ENTREZID column matches exactly to our dge_inital$counts rownames.
table(ann$ENTREZID==rownames(dge_initial$counts))

## add this information to expression dataset
exp_data$genes <- ann$SYMBOL

## filter expression dataset for risk genes
mm_risk <- exp_data %>% dplyr::filter(genes %in% riskgenes)

## filter for samples where clinical info also available 
mm_risk_expdata <- mm_risk[, colnames(mm_risk) %in% clin_data$RNASeq_transLevelExpFileSamplId]

## bring columns in expression data in the correct order as in the clinical dataset
mm_risk_expdata <- mm_risk_expdata[clin_data$RNASeq_transLevelExpFileSamplId]

##rename entrez id and gene columns
mm_risk_expdata$Entrez_Id <- mm_risk$X
mm_risk_expdata$Gene <- mm_risk$genes

## make new countdata using the filtered dataset
countdata <- mm_risk_expdata[,-(564:565)]
# Store EntrezGeneID as rownames
rownames(countdata) <- mm_risk_expdata[,564]

##Make sure that the column names are now the same as sample name in the clinical data file. This is good because it means our sample information in clinical data is in the same order as the columns in countdata
table(colnames(countdata)==clin_data$RNASeq_transLevelExpFileSamplId)

## make a dgelist object for visualizations, qc checks
dge <- DGEList(countdata)

## add clinical covariates (age+sex) to DGEList object
group_age <- factor(clin_data$agegrp)
group_sex <- factor(clin_data$D_Gender)

# Add the group information into the DGEList
dge$samples$age <- group_age
dge$samples$gender <- group_sex

## qc checks
## library size and distribution
barplot(dge$samples$lib.size/1e06, names=colnames(dge), las=2, ann=FALSE, cex.names=0.2)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# examine the distributions of the raw counts using log of the counts.
# Get log2 counts per million
logcounts <- cpm(dge,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, cex.axis=0.2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

## multi-dimensional scaling plots

# set up colour schemes for age and gender

col.sex <- c("purple","orange")[as.factor(clin_data$D_Gender)]
col.age <- c("blue","red","black","green", "purple", "orange")[clin_data$agegrp]

#MDS with disease stage colouring
plotMDS(dge,col=col.sex)
# Let's add a legend to the plot so we know which colours correspond to which stage
legend("topleft",fill=c("purple","orange"),legend=levels(as.factor(clin_data$D_Gender)))
# Add a title
title("Sex")

# Similarly for age

plotMDS(dge,col=col.age)
legend("topleft",fill=c("blue","red","black","green", "purple", "orange"),legend=levels(clin_data$agegrp),cex=0.5)
title("Age")

#############################################
## MERGE EXP AND CLIN DATASETS
#############################################

## reshape the expression dataset so that genes become column names
final_exp_data <- mm_risk_expdata[,-564] %>%
  gather(key = key, value = value, 1:563) %>%
  spread(key = names(mm_risk_expdata)[565], value = "value")

## rename the sample id column to merge with clinical data in the next step
final_exp_data <- final_exp_data %>% dplyr::rename(RNASeq_transLevelExpFileSamplId=key)

## merge clinical and expression data
clin_exp_merge <- merge(x=clin_data, y=final_exp_data, by="RNASeq_transLevelExpFileSamplId", all.x=T)

## change column names that are not suitable for R
clin_exp_merge <- clin_exp_merge %>% dplyr::rename(H1_4 = "H1-4")

## correlation of genes with event time
corres_genes <- cor(clin_exp_merge[,c(5, 12:39)], method="pearson", use="pairwise.complete.obs")

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(corres_genes)

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# plot Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  geom_text(aes(label = round(value, 2)), size=1.5)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()



#########################
## MODELLING
#########################

## remove irrelevant columns from the dataset
train_data <- clin_exp_merge %>% dplyr::select(-c("RNASeq_transLevelExpFileSamplId", "Patient", "agegrp", "D_Gender"))

## Feature exploration: Attribution of each predictor to the predicted  risk
rfi=generateFilterValuesData(surv.task,method=c("randomForestSRC_importance"))%>%.$data%>%ggplot(aes(x=reorder(name,value),y=value,fill=reorder(name,value)))+geom_bar(alpha=0.8,stat="identity",color="black",show.legend=F)+scale_x_discrete("Features")+coord_flip()


uv=generateFilterValuesData(surv.task,imp.learner=makeLearner("surv.coxph"),method="univariate.model.score")%>%.$data%>%ggplot(aes(x=reorder(name,value),y=value,fill=reorder(name,value)))+geom_bar(alpha=0.8,stat="identity",color="black",show.legend=F)+scale_x_discrete("Features")+coord_flip()


grid.arrange(rfi,uv)

## Results show some features are not important predictors, however, they still might contribute to the mortality risk (based on the univariate model). 
## Thus, all features from the dataset will be included in our model.

## make survival task
surv.task <- makeSurvTask(data=train_data, target = c("D_PFS", "HR_FLAG"))


## z-score normalize all numeric features 
surv.task <- normalizeFeatures(surv.task, method="standardize", cols=c("D_Age", "ACTG1", "ARID1A", "ATM", 
                                                                       "BRAF", "CCND1", "CREBBP", "CSGALNACT1", "CYLD", "DIS3", "EGR1", 
                                                                       "EPAS1", "ERC2", "FGFR3", "H1_4", "HUWE1", "IRF4", "KMT2C", "KRAS", 
                                                                       "MAX", "NRAS", "NSD2", "PHF19", "PRC1", "PRKD2", "TENT5C", "TP53", 
                                                                       "TRAF3", "UBR5"))

#Benchmark study : comparing 2 different algorithms for survival task

# Creating a list of 2 learners
svlearners=list(
  makeLearner("surv.coxph"),
  makeLearner("surv.randomForestSRC")
)

rdesc=makeResampleDesc("CV", iters=5, predict = "both")

# Initalising the Benchmark study
set.seed(123)
svbnmrk=benchmark(svlearners,surv.task,rdesc)

bmrkdata=getBMRPerformances(svbnmrk, as.df = TRUE)%>%.[,-1]

## plot and check the performance of the two learners
## Distribution of C-index of 2 survival algorithms, based on a 5-fold cross-validation
plotBMRRanksAsBarChart(svbnmrk)+scale_fill_manual(values=myfillcolors)
plotBMRBoxplots(svbnmrk)+aes(fill=learner.id)+coord_flip()+scale_fill_manual(values=myfillcolors,name="Learners")
ggplot(bmrkdata)+geom_path(aes(x=iter,y=cindex,color=learner.id),size=1,alpha=0.8)+facet_wrap(~learner.id,scales="free",ncol=2)+scale_color_manual(values=myfillcolors)

## both learners perform equally well, proceed with the random forest SRC model

## make final learner (random forest SRC)
svlearner <- makeLearner("surv.randomForestSRC")

## set parameter space
params <- makeParamSet(makeDiscreteParam("mtry", c(2,4,6,8,10)), makeDiscreteParam("ntree", c(3,10,50,100,300,1000)))

## define type of search (GridsearchCV)
ctrl <- makeTuneControlGrid()

## define CV scheme (5-fold CV)
rdesc=makeResampleDesc("CV", iters=5, predict = "both")

set.seed(2356)

## tune hyperparameters
tune <- tuneParams(learner=svlearner, task=surv.task, resampling=rdesc, par.set=params, control = ctrl,
                   show.info=T)

## change params (take only optimal params)
svlearner_tune <- setHyperPars(makeLearner("surv.randomForestSRC"), mtry=tune$x$mtry, ntree=tune$x$ntree)

svlearner_tune$par.vals <- list(importance=T)

## resampling - to anticipate how the model might possibly perform on a new dataset
## the CV scheme in mlr by default makes 80/20 split of the training dataset
cv <- resample(svlearner_tune, surv.task, rdesc, measures=cindex, models=T, extract=getFeatureImportance)

#train a model
rforestsrc <- train(svlearner_tune, surv.task)
getLearnerModel(rforestsrc)

## performance evaluation - check c-index across 5 iterations
c_index_df <- cv$measures.test

##plot
ggplot(c_index_df)+geom_path(aes(x=iter,y=cindex),size=1,alpha=0.8)+scale_color_manual(values=myfillcolors)

## get feature importance using the cross validation set
imp_df_cv <- rbind(cv$extract[[1]]$res,
                   cv$extract[[2]]$res,
                   cv$extract[[3]]$res,
                   cv$extract[[4]]$res,
                   cv$extract[[5]]$res)%>% dplyr::group_by(variable) %>% dplyr::summarise(mean=mean(importance), sd=sd(importance)) %>% arrange(desc(mean))

## plot
ggplot(imp_df_cv, aes(x=reorder(variable, mean), y=mean, fill=mean))+
  geom_bar(stat="identity", position="dodge")+coord_flip()+
  ylab("Variable importance")+xlab("")+ggtitle("Information Value Summary")+
  guides(fill=F)+scale_fill_gradient(low="red", high="blue")+theme_minimal()

## Kaplan-Meier curves for the top 5 genes that contribute as important predictors of survival risk

## first:divide patients into two groups using gene expression median as a cutoff
clin_exp_merge = clin_exp_merge %>%
  mutate(PHF19_exp = case_when(
    PHF19 > quantile(PHF19, 0.5) ~ 'PHF19_High',
    PHF19 < quantile(PHF19, 0.5) ~ 'PHF19_Low',
    TRUE ~ NA_character_
  ))

clin_exp_merge = clin_exp_merge %>%
  mutate(PRC1_exp = case_when(
    PRC1 > quantile(PRC1, 0.5) ~ 'PRC1_High',
    PRC1< quantile(PRC1, 0.5) ~ 'PRC1_Low',
    TRUE ~ NA_character_
  ))

clin_exp_merge = clin_exp_merge %>%
  mutate(MAX_exp = case_when(
    MAX > quantile(MAX, 0.5) ~ 'MAX_High',
    MAX< quantile(MAX, 0.5) ~ 'MAX_Low',
    TRUE ~ NA_character_
  ))

clin_exp_merge = clin_exp_merge %>%
  mutate(CCND1_exp = case_when(
    CCND1 > quantile(CCND1, 0.5) ~ 'CCND1_High',
    CCND1< quantile(CCND1, 0.5) ~ 'CCND1_Low',
    TRUE ~ NA_character_
  ))

clin_exp_merge = clin_exp_merge %>%
  mutate(CSGALNACT1_exp = case_when(
    CSGALNACT1> quantile(CSGALNACT1, 0.5) ~ 'CSGALNACT1_High',
    CSGALNACT1< quantile(CSGALNACT1, 0.5) ~ 'CSGALNACT1_Low',
    TRUE ~ NA_character_
  ))

fit_phf19 = survival::survfit(Surv(D_PFS, HR_FLAG) ~ PHF19_exp, data = clin_exp_merge)
fit_prc1 = survival::survfit(Surv(D_PFS, HR_FLAG) ~ PRC1_exp, data = clin_exp_merge)
fit_max = survival::survfit(Surv(D_PFS, HR_FLAG) ~ MAX_exp, data = clin_exp_merge)
fit_ccnd1 = survival::survfit(Surv(D_PFS, HR_FLAG) ~ CCND1_exp, data = clin_exp_merge)
fit_csg = survival::survfit(Surv(D_PFS, HR_FLAG) ~ CSGALNACT1_exp, data = clin_exp_merge)

ggsurvplot(fit_phf19, pval = T, conf.int = T)
ggsurvplot(fit_prc1, pval = T, conf.int = T)
ggsurvplot(fit_max, pval = T, conf.int = T)
ggsurvplot(fit_ccnd1, pval = T, conf.int = T)
ggsurvplot(fit_csg, pval = T, conf.int = T)


## save the model
saveRDS(rforestsrc, "pred_model.rds")


#########################
## DO NOT RUN
############################

## guide to test it on an external dataset
# load the model
pred_model <- readRDS("pred_model.rds")
print(super_model)
# make a predictions on "new data" using the final model
final_predictions <- predict(pred_model, validation_data)###

