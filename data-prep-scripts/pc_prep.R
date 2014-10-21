#================================================================================*
# ---- BIRD PRE-PROCESSING ----
#================================================================================*
# The goal of this script is to make the bird data ready for analysis.

# Libraries:

library(lubridate)
library(plyr)

# Load data:

setwd('/Users/bsevans/Desktop/community')

p1 = read.csv('sampling_data/pc_2009.csv')      # Records 2009
p2 = read.csv('sampling_data/pc_2010.csv')      # Records 2010   
p3 = read.csv('sampling_data/2012 Nestwatch Data.csv')  # Records 2012
p4 = read.csv('sampling_data/pc_supplement_2010.csv')   # Supplemental point counts

lc = read.csv('sampling_data/derived_data/pts_lc.csv')

#---------------------------------------------------------------------------------*
# ---- Bind data ----
#=================================================================================*

# Fixing site codes for 2009 records, remove and rename columns:

p1 = merge(p1,lc,by = c('X','Y'), all = F)

p1 = p1[,-c(1:6,20:22,24)]

p1 = p1[,c(14,1:13)]

names(p1)[1] = 'Site'

# Bind p2-p4 (similar column names):

p2 = p2[,-21]

pc = rbind.fill(p2,p3,p4)

# Remove start time (not available in p1) and unnecessary fields:

pc = pc[,-c(4,16:20)]

# Bind with p1:

pc = rbind(p1,pc)

#---------------------------------------------------------------------------------*
# ---- Format dates ----
#=================================================================================*

# Extract year of sample and provide decimal date:

pc$year = year(as.Date(pc$Date, format = '%m/%d/%Y'))
month = month(as.Date(pc$Date, format = '%m/%d/%Y'))
day = day(as.Date(pc$Date, format = '%m/%d/%Y'))
pc$ddate = ifelse(month == 5, 5 + day/31, ifelse(month == 6, 6 + day/30, 
              ifelse(month == 7, 7 + day/31,8 + day/31)))

# Remove Date column and reorder:

pc = pc[,c(1,15:16,3:14)]

#---------------------------------------------------------------------------------*
# ---- Technician names ----
#=================================================================================*

unique(pc$Observer)

pc$Observer = tolower(pc$Observer)

pc$Observer = gsub('reitsma','rr',pc$Observer) # Bob
pc$Observer = gsub('jw','jlw',pc$Observer) # Josie

pc = pc[pc$Observer != '',] # Remove blanks

# Collapse multiple observer records:

pc$Observer = gsub('cw,jsm,jlw','multiple',pc$Observer) 
pc$Observer = gsub('all interns and brian','multiple',pc$Observer)
pc$Observer = gsub('kam, jag, jlw, amc','multiple',pc$Observer)
pc$Observer = gsub('kam, jag, jlw, amc bse','multiple',pc$Observer)
pc$Observer = gsub('jag, kam, bse','multiple',pc$Observer)
pc$Observer = gsub('amc, kam','multiple',pc$Observer)
pc$Observer = gsub('jag, bse','multiple',pc$Observer)
pc$Observer = gsub('amc, bse','multiple',pc$Observer)
pc$Observer = gsub('jag, kam, lac','multiple',pc$Observer)
pc$Observer = gsub('jag, kam','multiple',pc$Observer)
pc$Observer = gsub('multiple bse','multiple',pc$Observer)

#---------------------------------------------------------------------------------*
# ---- Species names and removals ----
#=================================================================================*

pc$Species = tolower(pc$Species)

pc = pc[pc$Species!='',] # Remove blanks

pc = pc[pc$Species!='-',] # Remove dashes

pc = pc[pc$Species!='--',] # Remove double-dashes

pc = pc[pc$Species!='-9999',] # Remove -9999 NA

pc$Species = gsub('amco','amro',pc$Species) # American robin typo

pc$Species = gsub('amrc','amro',pc$Species) # American robin typo

pc$Species = gsub('armo','amro',pc$Species) # American robin typo

pc$Species = gsub('amgf','amgo',pc$Species) # American goldfinch typo

pc$Species = gsub('basw','bars',pc$Species) # Barn swallow typo

pc$Species = gsub('bhci','bhco',pc$Species) # Brown-headed cowbird typo

pc$Species = gsub('brhc','bhco',pc$Species) # Brown-headed cowbird typo

pc$Species = gsub('btbl','btbw',pc$Species) # Assuming Black-throated blue typo

pc$Species = gsub('carw\n','carw',pc$Species) # I have no idea what this notation means

pc$Species = gsub('cawr','carw',pc$Species) # CARW typo

pc$Species = gsub('cewa','cedw',pc$Species) # Cedar waxwing

pc = pc[pc$Species!='chwv',] # I have no idea what species this refers to

pc$Species = gsub('cogr\n','cogr',pc$Species) # Common grackle

pc = pc[pc$Species!='cona',] # I have no idea what species this refers to

pc$Species = gsub('crga','grca',pc$Species) # Gray catbird

pc$Species = gsub('eape','eawp',pc$Species) # Eastern bluebird

pc$Species = gsub('ebbl','eabl',pc$Species) # Eastern bluebird

pc$Species = gsub('etti','tuti',pc$Species) # Tufted titmouse

pc$Species = gsub('ewpe','eawp',pc$Species) # Eastern Wood-Pewee

pc$Species = gsub('grca\n','grca',pc$Species) # I have no idea what this notation means

pc$Species = gsub('grfl','gcfl',pc$Species) # great-crested flycatcher

pc$Species = gsub('hofi\n','hofi',pc$Species) # I have no idea what this notation means

pc$Species = gsub('hosp ','hosp',pc$Species) # An extra space in the name

pc$Species = gsub('hosp\n','hosp',pc$Species) # I have no idea what this notation means

pc$Species = gsub('hsop','hosp',pc$Species) # House sparrow

pc$Species = gsub(' noca','noca',pc$Species) # An extra space before the name

pc$Species = gsub('noca ','noca',pc$Species) # An extra space after the name

pc$Species = gsub('noca\n','noca',pc$Species) # I have no idea what this notation means

pc$Species = gsub('nofl','ysfl',pc$Species) # Yellow-shafted flicker

pc = pc[pc$Species!='osfl',] # SOOO doubtful that an olive-sided flycatcher was seen

pc = pc[pc$Species!='rewo',] # I have no idea what species this refers to

pc$Species = gsub('rbth','rthu',pc$Species) # Rock pigeon

pc$Species = gsub('rodo','ropi',pc$Species) # Rock pigeon

pc$Species = gsub('sosp ','sosp',pc$Species) # An extra space after the name

pc$Species = gsub('sosp\n','sosp',pc$Species) # I have no idea what this notation means

pc$Species = gsub('wbu','wbnu',pc$Species) # White-breated Nuthatch

pc = pc[pc$Species!='ybsa',] # Migrants

# Removing some values (e.g., UNID, water birds (canada goose), cats, chipmunks, etc.):

pc = pc[pc$Species!='cang',] # Remove Canada Goose records

pc = pc[pc$Species!='cat',] # Remove cat records

pc = pc[pc$Species!='chipmunk',] # Remove chipmunk records

pc = pc[pc$Species!='lagu',] # Remove Laughing Gull records

pc = pc[pc$Species!='ospr',] # Remove Osprey records

pc = pc[pc$Species!='rabbit',] # Remove rabbit

pc = pc[pc$Species!='squirrel',] # Remove squirrel

pc = pc[pc$Species!='Squirrel',] # Remove squirrel

pc = pc[pc$Species!='squirrell',] # And now introducing ... Ali's many spellings of "squrill"

pc = pc[pc$Species!='squriel',] # Remove squirrel

pc = pc[pc$Species!='squril',] # Remove squirrel

pc = pc[pc$Species!='squrill',] # Remove squirrel

pc = pc[pc$Species!='squrril',] # Remove squirrel

pc = pc[pc$Species!='ukn. woodpecker',] # Remove UNID's

pc = pc[pc$Species!='unbi',] # Remove UNID's

pc = pc[pc$Species!='undu',] # Remove UNID's

pc = pc[pc$Species!='unid',] # Remove UNID's

pc = pc[pc$Species!='unid-x',] # Remove UNID's

pc = pc[pc$Species!='unid (woodpecker knocking)',] # Remove UNID's

pc = pc[pc$Species!='unid gull',] # Remove UNID's

pc = pc[pc$Species!='unk',] # Remove UNID's

pc = pc[pc$Species!='unk young',] # Remove UNID's

pc = pc[pc$Species!='unkn',] # Remove UNID's

pc = pc[pc$Species!='unwo',] # Remove UNID's

pc = pc[pc$Species!='unwa',] # Remove UNID's

pc = pc[pc$Species!='wodu',] # Remove Wood Duck

# REMOVALS VERIFIED WITH EBIRD RECORDS FOR THE REGION:

pc = pc[pc$Species != 'brcr',] # Keep? 1 data point, but mid-June ### EBIRD: REMOVE ###

pc = pc[pc$Species != 'btbw',] # Keep? 1 data point, but mid-June ### EBIRD: REMOVE ###

pc = pc[pc$Species != 'cswa',] # Remove, 1 record early June ### EBIRD: REMOVE ###

pc = pc[pc$Species != 'swth',] # Remove, 1 record mid June ### EBIRD: REMOVE ###

pc = pc[pc$Species != 'ssha',] # Remove? 1 record, but in mid July ### EBIRD: REMOVE ###

# Remove the Species that are blatantly just passing through (eBird verified):

pc = pc[pc$Species!='blpw'&pc$Species!='brcr'&pc$Species!='btbw'
         &pc$Species!='cswa'&pc$Species!='swth'&pc$Species!='ssha',]

# BIRDS WITH JUST ONE OR TWO RECORDS BUT DO BREED IN DC (eBird sightings):

pc = pc[pc$Species != 'amre',] # Remove amre, one sighting, early June

pc = pc[pc$Species != 'coha',] # Hard to say, 1 record but in early July

#pc[pc$Species != 'coye',]  # Interestingly, this Species has been only observed at 1 site, but 3 years apart!

pc = pc[pc$Species != 'howa',] # Remove? 1 record, but mid-July

pc = pc[pc$Species != 'kill',] # Keep? Probably not. Just one record in early June

pc = pc[pc$Species != 'prow',] # Keep? Probably not. Just one record in mid June

pc = pc[pc$Species != 'wavi',] # Keep? Probably not, though the one record was in late July

pc = pc[pc$Species != 'wewa',] # Keep? Probably not, one record, early June

pc = pc[pc$Species != 'ybch',] # Remove? 1 record late June

pc = pc[pc$Species != 'ybcu',] # Remove? 1 record, but in mid July

pc = pc[pc$Species != 'puma',] # These are just purple martin boxes!

unique(pc[order(pc$Species),]$Species)

#---------------------------------------------------------------------------------*
# ---- Standardization and NA hunting ----
#=================================================================================*
# Simplify names:

names(pc)[12:15] = c('detection','temp','sky','wind')

# Change to lower case:

names(pc) = tolower(names(pc))

pc$detection = tolower(pc$detection)

# summary(pc$year)
# 
# summary(pc$wind)
# 
# pc[is.na(pc$temp),] 

# pc[pc$wind<0,'wind']
# pc[is.na(pc$d10),'d10'] = 0
# head(pc)
#---------------------------------------------------------------------------------*
# ---- Names and writing file ----
#=================================================================================*

# Set distance names:

names(pc)[7:11] = c('d10','d20','d30','d40','d50')

# Set NA's in distance columns to 0's:

pc[is.na(pc$d10),'d10'] = 0
pc[is.na(pc$d20),'d20'] = 0
pc[is.na(pc$d30),'d30'] = 0
pc[is.na(pc$d40),'d40'] = 0
pc[is.na(pc$d50),'d50'] = 0



# Write file:

write.csv(pc,'sampling_data/derived_data/pc_10_14.csv', row.names = F)
