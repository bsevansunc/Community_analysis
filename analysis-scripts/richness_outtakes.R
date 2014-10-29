

# Create a dataframe of expected richness values by potentential combinations of
# impervious and canopy:

expected.frame = function(response,imp, can){
  # Create modeled values of imp and can:
  imp = imp
  can = can
  # Make grid:
  grid = mesh(imp,can)
  # Exclude values for which the percent combination exceeds 100
  g1 = grid$x + grid$y
  g1 = ifelse(g1>100, NA, 1)
  # Project values:
  z = with(grid, expected.values(sr,imp,can))
  # Limit to potential combinations of imp and can:
  z = z*g1
  # Create a dataframe for output
  out.df = melt(z)
  out.df = na.omit(out.df)
  names(out.df) = c('imp','can','richness')
  return(out.df)
}

expected.frame(sr, lc$imp,lc$can)

R2 <- cor(sr,e1)^2

plot(lc$imp, sr)
plot(lc$can, sr)

plot(lc$imp, e1)
plot(lc$can, e1)

# Function that takes the model average frame, imp and canopy and
# calculates the expected values (global model):

expected.values = function(response, imp, can){
  mod = summary.outs(response)[[2]][,2]
  mod[1]+mod[2]*imp+mod[3]*can+mod[4]*imp*can+mod[5]*(imp^2)+mod[6]*(can^2)
}

# Function to determine the adjusted R2 averaged across the model:

adj.r2 = function(response, imp, can){
  aic.table = summary.outs(sr)[[1]]
  observed = response
  expected = expected.values(sr, imp, can)
  r2 = cor(observed, expected)^2
  k.ave = sum(aic.table$k*aic.table$w)
  ((length(response)-1)*r2-k.ave)/(length(response)-1-k.ave)
}
