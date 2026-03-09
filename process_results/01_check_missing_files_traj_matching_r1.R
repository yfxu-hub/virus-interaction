# ---------------------------------------------------------------------------------------------------------------------
# Code to determine whether any results are missing after running cluster code
# ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# Load libraries:
library(tidyverse)
library(testthat)
library(readxl)
vir1_name<-'h1'
vir2_name<-'rsv'

vir1<-vir1_name
vir2<-vir2_name
# Get cluster environmental variables:
#fit_canada <- as.logical(Sys.getenv("FITCANADA")); print(fit_canada)
fit_canada<-FALSE
# Set size of parameter start space:
sobol_size <- 500

# Get list of completed runs for each virus:
res_files <- list.files(path = '/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/round1', 
                        pattern = sprintf('res_%s_%s',vir1, vir2), full.names = TRUE)

check_list <- list.files(path = '/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/round1', 
                         pattern = sprintf('res_%s_%s',vir1, vir2), full.names = TRUE)

# Get vector of seasons for which fitting was done:
# if (fit_canada) {
#   all_yrs <- unique(str_sub(check_list, 13, 18))
# } else {
#   all_yrs <- unique(str_sub(check_list, 23, 28))
# }
# print(all_yrs)
# print(length(all_yrs))

# Set parameters for run:
#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_allperiod.xlsx")%>%as.data.frame()
outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1_with_flursv.xlsx")%>%as.data.frame()
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]

all_yrs<-outbreak$season
#all_yrs<-c('s22-23','s23-24','s24-25')

# Check for complete results:
if (length(check_list) != sobol_size * length(all_yrs) & length(check_list) > 0) {
  
  for (yr in all_yrs) {
    #yr<-"s23"
    temp_list <- check_list[str_detect(check_list, pattern = yr)]#
    
    if (length(temp_list) != sobol_size) {
      
      if (fit_canada) {
        
        completed_runs <- str_split(temp_list, '_') %>%
          lapply(., function(ix) {ix[12]}) %>%
          unlist() %>%
          str_split(fixed('.')) %>%
          lapply(., function(ix) {ix[1]}) %>%
          unlist() %>%
          as.numeric()
        
      } else {
        
        completed_runs <- str_split(temp_list, '_') %>%
          lapply(., function(ix) {ix[11]}) %>%
          unlist() %>%
          str_split(fixed('.')) %>%
          lapply(., function(ix) {ix[1]}) %>%
          unlist() %>%
          as.numeric()
        
      }
      
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_h1_plus_b/RSV in ', yr, ', missing run: ', run))
      }
      
    }
    
  }
}

# Clean up:
#rm(list = ls())
print('Done.')

