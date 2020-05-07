###############################
# R helpers for NHC_2019_metabolism project

# run length encoder: 
# for finding start and end points of NA sequences for interpolating discharge

rle_custom = function(x){
  r = rle(x)
  ends = cumsum(r$lengths)
  r = as.data.frame(cbind(values=r$values,
                         starts=c(1, ends[-length(ends)] + 1),
                         stops=ends, lengths=r$lengths, deparse.level=1))
  return(r)
}
