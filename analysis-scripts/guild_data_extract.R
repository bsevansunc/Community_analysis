# Copy-and-paste data from http://www.esapubs.org/archive/ecol/E095/178/BirdFuncDat.txt
# Note: See metadata at: http://www.esapubs.org/archive/ecol/E095/178/metadata.php and study
# at http://www.esajournals.org/doi/abs/10.1890/13-1917.1

#d = read.table('clipboard', header=T, sep = '\t', quote = "\", fill=T)

d = read.table('clipboard', header=T, sep = '\t', quote = "\"", fill=T)

d = read.table('guild_data/diet_foraging', header=T, sep = '\t', quote = "\"", fill=T)
