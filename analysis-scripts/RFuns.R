# A function for matching or not matching records in a data frame
# Note: Modified from Brian Ripley

dfMatcher <- function(df1,df2, match.or.no){
  v1 <- do.call("paste", df1)
  v2 <- do.call("paste", df2)
  if (match.or.no == 'match')
  {df1[v1 %in% v2,]} else df1[!v1 %in% v2,]
}

# A function for determining to return julian dates from a vector of
# dates.

jDate = function(date){
  jdate = numeric()
  for(i in 1:length(date)){
  jdate[i] = julian(date[i])[1]
  }
  return(jdate)
  }

# A function to determine the day-of-year for a given date.
# Input date is a character string of format ('yyyy/mm/dd')

dayOfYear = function(date){
  y = format(as.Date(date), "%Y")
  d0 = as.Date(paste(y, '1/1', sep = '/')) 
  d1 = as.Date(date)
  as.numeric(julian(d1, origin = d0))+1
  #return(j1 - j0 + 1)[1]
}

# Function to change a factor to a proper number:

as.numericF = function(x) as.numeric(as.character(x))

t = as.factor(c(3,2,5))
