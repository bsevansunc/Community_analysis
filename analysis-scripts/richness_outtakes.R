

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