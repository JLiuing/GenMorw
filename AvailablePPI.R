AvailablePPI <- function()
{
  ppis <- list.files("Data/networks_rds", pattern = "\\.rds$")
  ppis <- unlist(strsplit(ppis,split = ".rds"))
  return(ppis)
}
