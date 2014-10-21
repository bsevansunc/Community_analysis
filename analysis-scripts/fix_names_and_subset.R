#---------------------------------------------------------------------------------------
# Function to fix all spec names in counts (all mispellings included):
#---------------------------------------------------------------------------------------

fix.sp.names = function(df){
  df$spec = gsub('amco','amro',df$spec) # American robin typo
  df$spec = gsub('amrc','amro',df$spec) # American robin typo
  df$spec = gsub('amgf','amgo',df$spec) # American goldfinch typo
  df$spec = gsub('armo','amro',df$spec) # American Robin typo
  df$spec = gsub('basw','bars',df$spec) # Barn swallow typo
  df$spec = gsub('bhci','bhco',df$spec) # Brown-headed cowbird typo
  df$spec = gsub('brhc','bhco',df$spec) # Brown-headed cowbird typo
  df$spec = gsub('btbl','btbw',df$spec) # Assuming Black-throated blue typo
  df$spec = gsub('carw\n','carw',df$spec) # I have no idea what this notation means
  df$spec = gsub('cawr','carw',df$spec) # CARW typo
  df = df[df$spec!='chwv',] # I have no idea what spec this refers to
  df$spec = gsub('cogr\n','cogr',df$spec) # Common grackle
  df = df[df$spec!='cona',] # I have no idea what spec this refers to
  df$spec = gsub('crga','grca',df$spec) # Gray Catbird typo
  df$spec = gsub('eape','eawp',df$spec) # Eastern Wood-Peewee
  df$spec = gsub('ebbl','eabl',df$spec) # Eastern bluebird
  df$spec = gsub('etti','tuti',df$spec) # Tufted titmouse
  df$spec = gsub('EUST','eust',df$spec) # European starling
  df$spec = gsub('ewpe','eawp',df$spec) # Eastern Wood-Pewee
  df$spec = gsub('grca\n','grca',df$spec) # I have no idea what this notation means
  df = df[df$spec!="grfl",] # I highly doubt a gray flycatcher was observed
  df$spec = gsub('hofi\n','hofi',df$spec) # I have no idea what this notation means
  df$spec = gsub('hosp ','hosp',df$spec) # An extra space in the name
  df$spec = gsub('hosp\n','hosp',df$spec) # I have no idea what this notation means
  df$spec = gsub('hsop','hosp',df$spec) # House Sparrow
  df$spec = gsub(' noca','noca',df$spec) # An extra space before the name
  df$spec = gsub('noca ','noca',df$spec) # An extra space after the name
  df$spec = gsub('noca\n','noca',df$spec) # I have no idea what this notation means
  df = df[df$spec!='osfl',] # SOOO doubtful that an olive-sided flycatcher was seen
  df = df[df$spec!='rewo',] # I have no idea what spec this refers to
  df$spec = gsub('rbth','rthu',df$spec) # Ruby-throated Hummingbird
  df$spec = gsub('rodo','ropi',df$spec) # Rock pigeon
  df$spec = gsub('sosp ','sosp',df$spec) # An extra space after the name
  df$spec = gsub('sosp\n','sosp',df$spec) # I have no idea what this notation means
  df$spec = gsub('squirrell','squirrel',df$spec) # Ali's many spellings of "squrill"
  df$spec = gsub('squriel','squirrel',df$spec) # Ali's many spellings of "squrill"
  df$spec = gsub('squril','squirrel',df$spec) # Ali's many spellings of "squrill"
  df$spec = gsub('squrill','squirrel',df$spec) # Ali's many spellings of "squrill"
  df$spec = gsub('squrril','squirrel',df$spec) # Ali's many spellings of "squrill"
  df$spec = gsub('ukn. woodpecker','unbi',df$spec) # UNBI
  df$spec = gsub('undo','unbi.w',df$spec) # UNBI
  df$spec = gsub('unid','unbi',df$spec) # UNBI
  df$spec = gsub('unbi-x','unbi',df$spec) # UNBI
  df$spec = gsub('unbi (woodpecker knocking)','unbi',df$spec) # UNBI
  df$spec = gsub('unbi gull','unbi.w',df$spec) # UNBI
  df$spec = gsub('unk','unbi',df$spec) # UNBI
  df$spec = gsub('unk young','unbi',df$spec) # UNBI
  df$spec = gsub('unbi young','unbi',df$spec) # UNBI
  df$spec = gsub('unkn','unbi',df$spec) # UNBI
  df$spec = gsub('unwa','unbi',df$spec) # UNBI
  df$spec = gsub('unwo','unbi',df$spec) # UNBI
  df$spec = gsub('unbin','unbi',df$spec) # UNBI
  df$spec = gsub('wbu','wbnu',df$spec) # White-breated Nuthatch
  df$spec = gsub('ysfl','nofl',df$spec) # Northern flicker
  df = df[df$spec!="",]
  df = df[df$spec!="-",]
  df = df[df$spec!="--",]
  df
}

#---------------------------------------------------------------------------------------
# Function to remove non-birds:
#---------------------------------------------------------------------------------------

fix.birds.only = function(df){

df = df[df$spec!='cat',] # Remove cat records

df = df[df$spec!='chipmunk',] # Remove chipmunk records

df = df[df$spec!='rabbit',] # Remove rabbit

df = df[df$spec!='squirrel',] # Remove squirrel

df

}

#---------------------------------------------------------------------------------------
# Function to include only identified birds:
#---------------------------------------------------------------------------------------

fix.no.unbi = function(df) {
  df = df[df$spec!="unbi",] 
  df[df$spec!='unbi.w',]
  }

# Function to only include land-birds:

fix.landbirds.only = function(df){
  df = df[df$spec!='cang',] # Remove Canada Goose records
  df = df[df$spec!='lagu',] # Remove Laughing Gull records
  df = df[df$spec!="ospr",] # Remove osprey
  df = df[df$spec!="puma",] # Remove purple martin
  df = df[df$spec!='wodu',] # Remove Wood Duck
}








