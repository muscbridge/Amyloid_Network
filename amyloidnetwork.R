
#load and clean data
library(readxl)
Network_mSUVr_results_240808 <- read_excel("~/Desktop/Network_mSUVr_results_240808.xlsx", 
                                           col_types = c("text", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", 
                                                         "numeric", "text", "text", "text", 
                                                         "text", "text", "text", "text", "text", 
                                                         "text", "text", "text", "text", "text", 
                                                         "text", "text"))
myvars<- c("studyid", "msuvr", "DAT_msuvr...14", "DMN_msuvr...15", "FPN_msuvr...16", "LIM_msuvr...17", "SOM_msuvr...18", "VAT_msuvr...19", "VIS_msuvr...20", "Nets_msuvr_avg");
amyloidbynetwork<-Network_mSUVr_results_240808[myvars];
amyloidbynetwork <- amyloidbynetwork[apply(amyloidbynetwork!=0, 1, all),] #remove obs with any zeros
amyloidbynetwork <- na.omit(amyloidbynetwork); #remove na
library(dplyr)
amyloidbynetwork <- amyloidbynetwork %>% 
  rename(DAT = DAT_msuvr...14, FPN = FPN_msuvr...16, LIM = LIM_msuvr...17, SOM = SOM_msuvr...18, VAT= VAT_msuvr...19, VIS = VIS_msuvr...20, DMN = DMN_msuvr...15, Network_Average = Nets_msuvr_avg)

# Find tertiles
vTert = quantile(amyloidbynetwork$msuvr, probs=2/3)
# classify values (errors but still works)
amyloidbynetwork$msuvr_groups <- ifelse(amyloidbynetwork$msuvr > vTert, "AB+", "AB-")

#subset data to plot
networks<-c("VIS", "SOM", "LIM", "DAT", "VAT", "FPN", "DMN", "Network_Average", "msuvr_groups")
plotdata<-amyloidbynetwork[networks]
library(reshape2)
plotdata_long<-melt(plotdata, id="msuvr_groups")
plotdata_long$value = as.double(plotdata_long$value);

library(ggplot2)
ggplot(plotdata_long, aes(x=variable, y=value, color=msuvr_groups))+
  geom_boxplot() +
  geom_point(position=position_jitterdodge(seed=4), alpha=.25) +
  ggtitle("Amyloid Deposition by Network") +
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(),legend.position = "bottom")+
  xlab("Network")+
  ylab("Mean Network Amyloid")+ 
  labs(color='Level of Global Amyloid') 

#medians
library(dplyr)
medians<- amyloidbynetwork %>%
group_by(msuvr_groups) %>%
  summarise(across(c(VIS, SOM, LIM, DAT, VAT, FPN, DMN), median))

#t-tests
amyloidbynetwork$FPNDMN_average <- rowMeans(amyloidbynetwork[, c("DMN", "FPN")])
amyloidbynetwork$VISSOM_average <- rowMeans(amyloidbynetwork[, c("VIS", "SOM")])
t.test(amyloidbynetwork$FPNDMN_average, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")
t.test(amyloidbynetwork$FPN, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")
t.test(amyloidbynetwork$DMN, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")
t.test(amyloidbynetwork$DAT, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")
t.test(amyloidbynetwork$VAT, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")
t.test(amyloidbynetwork$LIM, amyloidbynetwork$VISSOM_average, paired = TRUE, alternative = "two.sided")

#subset data to plot
networks<-c("VIS", "SOM", "LIM", "DAT", "VAT", "FPN", "DMN", "studyid", "msuvr_groups")
subdata<-amyloidbynetwork[networks]
subdata$cognet <- rowMeans(subdata[, c("LIM", "DAT", "VAT", "FPN", "DMN")])
subdata$smnet <- rowMeans(subdata[, c("VIS", "SOM")])
subsubdata<-subdata[c("studyid", "msuvr_groups", "cognet", "smnet")]
library(tidyr)
data_long<-gather(subdata, network, value, VIS:DMN, factor_key = TRUE)
data_long<-gather(subsubdata, network, value, cognet:smnet, factor_key = TRUE)
library(rstatix)
data_long$msuvr_groups<-as.factor(data_long$msuvr_groups)

#Full sample
res.aov1 <- anova_test(data = data_long, dv = value, wid =studyid, within = network) 
get_anova_table(res.aov1)

#which networks diff regardless of group
pwc1 <- data_long %>%   
  pairwise_t_test(value ~ network, paired = TRUE, p.adjust.method = "fdr") 
pwc1

#Does effect differ by group?
#repeated measures anova to test for interaction btw network X group
#shows main effects of network and group
#in general the AB+ are higher than AB-
#since there is a sig X, followup with pairwise comp
res.aov2 <- anova_test(data = data_long, dv = value, wid =studyid, between = msuvr_groups, within = network) 
get_anova_table(res.aov2)

#one-way anova to test if each network differs by group
#in general the AB+ are higher than AB- in every network
one.way <- data_long %>%   
  group_by(network) %>%   
  anova_test(dv = value, wid = studyid, between = msuvr_groups) %>%   
  get_anova_table() %>%   adjust_pvalue(method = "fdr") 
one.way

#pairwise comp, within each group which networks differ from each other
#little diff in AB-, in AB+ FPN and DMN were higher than SOM
pwc2 <- data_long %>%   
  group_by(msuvr_groups) %>%   
  pairwise_t_test(value ~ network, p.adjust.method = "fdr") 
pwc2





