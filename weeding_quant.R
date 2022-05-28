################################################################################
### R code for the analysis and graphing of ant weeding bioassay data for pub: 
##          Trachymyrmex septentrionalis ants promote fungus garden hygiene 
##          using Trichoderma-derived metabolite cues
################################################################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggsignif)
library(ggpubr)

# create colorblind friendly pallete
cborange="#E69F00";cbblue="#0072B2";cbdorange="#D55E00";cbgreen="#009E73"
cbyellow="#F0E442";cblblue="#56B4E9";cbpink="#CC79A7";cbblack="#000000";cbgray="#999999";
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # can also add black or gray

### CRUDE DATA ###
# files for crude analysis (slightly diff csv format for crude (more spread))
crude.df=read.csv("master_Su2020antExperiments-alldata.tidy.crude.ncdups.csv")  # reminder code

# careful bc filter by trashing success !n may not be objective
# removes half exp. 17, exp. 18, but keeps 19
filter(crude.df, humidity.perc == '100') %>% 
  filter(exp.hours == 24) %>% filter(trashing.success != 'n') %>% 
  filter(extract.id.fill != "MB1079") %>% 
  filter(extract.id.fill != "MB1082") -> fullfilt.df
# get crude mean
fullfilt.trash.mean.nodup$trash.g=c(0.3950425, 0.080814, 0.0368370) # numbers calculated in gsheet
###################

### FRACTION DATA ###
dat=read.csv("alldata.tidy.csv")

fullfilt.fracdat.df <- filter(dat,exp.num %in% c(11,12,15,19)) %>% 
  filter(!grepl("d crude",treat2)) %>%  # get rid of 3d7d14d
  filter(treat2 != "G") %>%             # G is a column wash not true frac
  mutate(parent.extract = case_when(
    grepl("MB1081", extract.id) ~ "MB1081",
    grepl("MB1084", extract.id) ~ "MB1084",
    extract.id == "" ~ "DMSO/NT"
    )
  )
# get fracs mean
fullfilt.fracdat.df %>% group_by(treat2) %>% 
  summarise(trash.g=mean(trash.g)) -> fullfilt.fracdat.trash.mean
#######################

### PURE COMPOUNDS (Oberlies) DATA ###
oberdat.df <- filter(dat,grepl("Oberlies",exp.id)) %>% 
  filter(treat1 != "transEKODE") %>%    # not a peptaibol
  filter(colony.id != 392)              # failed rep, PC fail

# get mean for treatments
oberdat.df %>% group_by(treat1) %>% 
  summarise(trash.g=mean(trash.g)) -> oberdat.trash.mean

######################################

##################
## Figure 3A graph
##################
p=ggplot(fullfilt.df, aes(x=treat1, y=trash.g)) +
  geom_line(aes(group=exp.col.ext.hum),size=0.5,col='grey90') + 
  geom_point(aes( fill = as.factor(colony.id)), shape=21, size=3.5)+
  geom_errorbar(data=fullfilt.trash.mean.nodup, aes(ymax=trash.g,ymin=trash.g,width=0.4),alpha=0.8)+
  xlab("") +ylab("Waste (g)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        panel.background = element_rect(fill = 'white',color='black'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=12,color = 'black'),
        axis.text.y=element_text(size=10,color='black')) + 
  scale_x_discrete(labels = c("NT", "DMSO","Crude")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  theme(legend.key=element_blank()) + 
  scale_fill_manual(values = cbPalette)+ 
  labs(fill = "Colony ID") + 
  theme(legend.background = element_blank(),
        legend.title = element_text(size=12),
        legend.justification=c(0,1),
        legend.position=c(0.025, 0.965),
        legend.text = element_text(size=12)) 
p$data$treat1 = factor(p$data$treat1, levels = c("na","dmso","crude"))
p


##################
## Figure 3B graph
##################
p=ggplot(fullfilt.fracdat.df, aes(x=treat2, y=trash.g)) + 
  geom_point(aes(fill=as.factor(colony.id)), shape=21,size=3.5) +
  geom_errorbar(data=fullfilt.fracdat.trash.mean, 
                aes(ymax=trash.g,ymin=trash.g,width=0.5),alpha=0.8) + 
  xlab("") + ylab("Waste (g)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        panel.background = element_rect(fill = 'white',color='black'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=12,color = 'black'),
        axis.text.y= element_text(size=12,color='black'), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        legend.key=element_blank(), 
        legend.background = element_blank(), 
        legend.justification=c(0,1),
        legend.position=c(0.01, 0.985)) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  scale_fill_manual(values = c(cblblue,cbgreen,cbblue,cbdorange)) + 
  labs(fill = "Colony ID", shape= "Parent Extract") +
  scale_x_discrete(labels = c("A","B","C","D","E","F","NT","DMSO","Crude"))
p$data$treat2 = factor(p$data$treat2, levels = c("A","B","C","D","E","F","na","dmso","crude"))
p

##################
## Figure S12 graph
##################
p=ggplot(fullfilt.df, aes(x=reorder(treat1, trash.g),y = trash.g)) + 
  geom_line(aes(group=exp.col.ext.hum),size=0.5,col='grey90') + 
  geom_point(aes(shape=extract.id.fill, fill=as.factor(colony.id)), size=3)+
  scale_shape_manual(values=c(22,23,24))+
  geom_errorbar(data=fullfilt.trash.mean.nodup, 
                aes(ymax=trash.g,ymin=trash.g,width=0.4),alpha=0.8) +
  xlab("") +ylab("Waste (g)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        panel.background = element_rect(fill = 'white',color='grey'), 
        axis.text.x = element_text(size=12,color = 'black'), 
        axis.text.y=element_text(size=12,color='black')) + 
  scale_x_discrete(labels = c("NT", "DMSO","Crude")) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  scale_fill_manual(values = cbPalette)+ 
  labs(fill = "Colony ID", shape= "Crude ID")
p$data$treat1 = factor(p$data$treat1, levels = c("na","dmso","crude"))
p


###################
## Figure S13 graph
###################
p=ggplot(fullfilt.fracdat.df, aes(x=treat2, y=trash.g)) + 
  geom_point(aes(shape=parent.extract, fill=as.factor(colony.id)), size=3.5) +
  scale_shape_manual(values=c(21,23,24)) + 
  geom_errorbar(data=fullfilt.fracdat.trash.mean, 
                aes(ymax=trash.g,ymin=trash.g,width=0.5),alpha=0.8) + 
  xlab("") + ylab("Waste (g)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        panel.background = element_rect(fill = 'white',color='black'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=12,color = 'black'),
        axis.text.y= element_text(size=12,color='black'), 
        legend.title = element_text(size=12),
        legend.text = element_text(size=12), 
        legend.key=element_blank(), 
        legend.background = element_blank(), 
        legend.justification=c(0,1),
        legend.position=c(0.01, 0.985)) + 
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  scale_fill_manual(values = c(cblblue,cbgreen,cbblue,cbdorange)) + 
  labs(fill = "Colony ID", shape= "Parent Extract") +
  scale_x_discrete(labels = c("A","B","C","D","E","F","NT","DMSO","Crude"))
p$data$treat2 = factor(p$data$treat2, levels = c("A","B","C","D","E","F","na","dmso","crude"))
p


###################
## Figure S14 graph
###################
ggplot(oberdat.df, aes(x=treat1, y=trash.g)) + 
  geom_point(aes(fill=as.factor(colony.id)), shape=21,size=4) +
  geom_errorbar(data=oberdat.trash.mean, 
                aes(ymax=trash.g,ymin=trash.g,width=0.5),alpha=0.8)  +
  theme(panel.background = element_rect(fill = 'white',color='black'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=12,color = 'black',angle = 30,vjust = 1, hjust=1), 
        axis.text.y=element_text(size=12,color='black'), 
        axis.title.y = element_text(size = 12), # don't need axis.title.x bc xlab is blank
        legend.key=element_blank(), 
        legend.background = element_blank(), 
        legend.justification=c(0,1),
        legend.position=c(0.01, 0.985),
        legend.title = element_text(size=12), 
        legend.text = element_text(size=12)) +
  scale_fill_manual(values = c(cbgreen,cbyellow,cbpink)) +
  xlab("") + ylab("Waste (g)") + labs(fill = "Colony ID") +
  guides(fill = guide_legend(override.aes = list(shape = 21)))
