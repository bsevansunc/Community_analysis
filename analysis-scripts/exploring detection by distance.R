#---------------------------------------------------------------------------------------
# Set-up
#---------------------------------------------------------------------------------------

# Libraries:

library(plyr)
library(ggplot2)

# Read in the file:

setwd('C:/Users/Brian/Documents/community_analysis/derived_data')
pc = read.csv('pc.csv')

# Remove unbi, water birds (e.g., Canada Goose), and non-birds (e.g., squirrels, cats, etc.):

source('C:/Users/Brian/Documents/community_analysis/r_scripts/fix_names_and_subset.R')

pc = fix.birds.only(pc)

pc = fix.no.unbi(pc)

pc = fix.landbirds.only(pc)

#---------------------------------------------------------------------------------------
# Explore distance data: Raw count proportions
#---------------------------------------------------------------------------------------

d.class = colSums(pc[,11:15], na.rm = T)/sum(colSums(pc[,11:15], na.rm = T))

barplot(d.class, ylim = c(0,.3), 
        ylab = 'Proportion of count', xlab = 'Distance class',
        cex.lab = 1.5)

# How about for a given impervious class?

dist.imp = function(imp1, imp2){
  df = pc[pc$imp>imp1&pc$imp<imp2,] # Subset dataframe
  n = dim(df)[1]
  sites = length(unique(df$site))
  d.class = colSums(df[,11:15], na.rm = T)/sum(colSums(df[,11:15], na.rm = T)) #props
  barplot(d.class, ylim = c(0,max(d.class)+.05), 
          ylab = 'Proportion of count', xlab = 'Distance class',
          main = paste('Count by distance', imp1, 'to', imp2, '% impervious surface'),
          cex.lab = 1.5)
  legend('topright',c(paste('n, birds =', n),paste('n, sites =', sites)), 
          bty = 'n', cex =1.25)
  }

dist.imp(-Inf, 1)
dist.imp(1, 10.01)
dist.imp(10, 30.01)
dist.imp(30,50.01)
dist.imp(50,70.01)
dist.imp(70,100)

#---------------------------------------------------------------------------------------
# Explore distance data: Average proportions
#---------------------------------------------------------------------------------------

sites = unique(pc[order(pc$site),c('site','imp')])

d.props = ddply(pc, .(site),d.func)

d.props = merge(sites, d.props)

tapply(d.props,)

mn.narm = function(x) mean(x, na.rm = T)

se = function(x) {
  x = x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

propt.func = function(imp1, imp2){
  df = d.props[d.props$imp>=imp1&d.props$imp<=imp2,] # Subset dataframe
  df2 = data.frame(matrix(ncol = 1))
  d.classes = seq(10,50, by = 10)
  means = c(mn.narm(df$d10),mn.narm(df$d20),
               mn.narm(df$d30),mn.narm(df$d40),
               mn.narm(df$d50))
  ses = c(se(df$d10),se(df$d20),se(df$d30),se(df$d40),se(df$d50))
  data.frame(d.classes,means,ses)
}


# The basic plot:

all.sites = propt.func(0, 100)

t2 = qplot(d.classes, y = means, data = all.sites,geom = c('point','line'))

# Define top and bottom of error bars:

eb = aes(ymax = all.sites$means + all.sites$ses, ymin = all.sites$means - all.sites$ses)
#eb90 = aes(ymax = all.sites$means + all.sites$ses*1.645, ymin = all.sites$means - all.sites$ses*1.645)
#eb95 = aes(ymax = all.sites$means + all.sites$ses*1.96, ymin = all.sites$means - all.sites$ses*1.96)

# Add error bars to the plot: 

t3  = t2 + geom_errorbar(data = all.sites, eb, color = 1, width = .5, size = .5,position=position_dodge(.9))# +
  #geom_errorbar(data = all.sites, eb90, color = eb.col, width = 0, size = 1.2,position=position_dodge(.9)) +
  #geom_errorbar(data = all.sites, eb95, color = eb.col, width = 0, size = .6,position=position_dodge(.9))
t3
# Add large points on top:

t3 = t3 + geom_point(size = 4) 

# Add plot labels:

plot.labs = list(x = 'Distance class', y = 'Mean proportion', title = 'Proportional counts across sites')

t3 = t3 + labs(plot.labs)

# Adjust the theme (size of labels)

t3+theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=20,face="bold", vjust = -.5),
        axis.title.y=element_text(size=20,face="bold", vjust = .3),
        plot.title = element_text(lineheight=.8, face="bold", size = 25, vjust = 1))

##########################################################################################
#

imp.dist = rbind(propt.func(0,100),propt.func(0,10),propt.func(10,40), propt.func(40,70),propt.func(70,100))
imp.dist$imp.class = rep(c('Across sites','IMP 0-10','IMP 10-40','IMP 40-70','IMP 70-100'), each = 5, times = 1)

imp.dist = imp.dist[,c(4,1:3)]
imp.dist = imp.dist[imp.dist!='Across sites',]
imp.dist = imp.dist[!is.na(imp.dist$ses),]

imp.class.plot = qplot(d.classes, y = means,data = imp.dist,
        geom = c('point','line'), color = imp.class)
imp.class.plot

eb = aes(ymax = means + ses, ymin = means - ses)
#imp.class.plot = imp.class.plot + geom_errorbar(data = imp.dist, eb, color = imp.class, width = .5, size = .5,
#                                                position=position_dodge(.9))
imp.class.plot  


plot.labs = list(x = 'Distance class', y = 'Mean proportion', title = 'Proportional counts by distance\nand impervious surface')

t3 = imp.class.plot + labs(plot.labs)

# Adjust the theme (size of labels)

t3 = t3+theme(axis.text=element_text(size=12),
         axis.title.x=element_text(size=20,face="bold", vjust = -.5),
         axis.title.y=element_text(size=20,face="bold", vjust = .3),
         plot.title = element_text(lineheight=.8, face="bold", size = 25, vjust = 1))

t3 + guides(fill = guide_legend(title  ="Impervious\nclass"))
t3 + theme(legend.text = element_text, size = 3)

# Adjust legend

t3 = t3  + theme(legend.title=element_blank())

t3 + scale_fill_continuous(guide = guide_legend(keywidth  = 40))

###########################################################################

pc$total = rowSums(pc[,11:15], na.rm = T)

d2 = ddply(pc, .(site, time.int,year),summarize, sum(total))

colnames(d2)[4] = 'count'

d2$rate.count = d2$count/d2$time.int

d2 = d2[!is.na(d2$time.int)&d2$time.int!='-9999',]

t2 = qplot(factor(time.int), y = rate.count, data = d2, geom = 'boxplot')
t2

###########################################################################

d2 = ddply(pc, .(site, year, jdate),summarize, sum(total))

hist(d2$jdate)

tapply(factor(d2$jdate),factor(d2$jdate),length)
colnames(d2)[4] = 'count'

d2$rate.count = d2$count/d2$time.int

d2 = d2[!is.na(d2$time.int)&d2$time.int!='-9999',]

t2 = qplot(factor(time.int), y = rate.count, data = d2, geom = 'boxplot')
t2

