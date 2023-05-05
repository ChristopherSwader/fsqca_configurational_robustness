
# We check for robustness to 
# (a) case selection, 
# (b) calibration decisions, 
# (c) necessity decisions, and 
# (d) truth table row inclusion score cutoffs. We check each of these simultaneously by checking all possible theoretically relevant combinations of the above.

testing_these_solution_paths <-  "CI*I*A*~PI + I*A*CP*P + I*A*P*PI"
necessary_conditions <- c("I", "A")


#Warmup####

  setwd()

library(BBmisc)
library(QCA)
library(SetMethods)
library(dplyr)
library(data.table)
library(stringr)

#load uncalibrated data
rawdt <- read.csv("uncalibratedreadytouse v9.csv") 
'Data format:
GEO	LL	CI	social_protection	rec_sport_culture	I	A	CP	P
Netherlands	1.22	946.31	6289	353	96	0.3221	2.48	1.3
Denmark	1.23	794.79	8629	372	93	0.2051	2.29	0.95
Finland	1.24	936.64	7914	304	90	0.1063	2.12	0.455
'


#load calibrated data
caldt <- read.csv("CALIBRATED v12.csv")
'Data format:
variable	e	c	i	valence
LL	1.49	1.39	1.25	-1
CI	550	1100	1600	1
social_protection	2700	5200	8300	1
'


#1. Robustness to case selection. Will randomly rotate through each of these to remove along with removing none.
cases <- rawdt$GEO

#2. Robustness to calibration choices. Will rotate through each of these.
#data in following format
'variable	e	c	i	valence
1	LL	1.49	1.39	1.25	-1
2	LL	1.49	1.34	1.25	-1
3	LL	1.6	1.39	1.25	-1
4	LL	1.6	1.34	1.25	-1
1	CI	550	1100	1600	1
2	CI	350	1100	1600	1
'

calibration_alternatives <- read.csv("calibration_alternatives v2.csv")


#3. Robustness to necessity decisions##


#4. Robustness to TT inclusion cutoffs...will generate this on the fly

#The robustness Loop.

grand_results_table <- data.table(matrix(nrow=0, ncol=9))
colnames(grand_results_table) <- c( "iteration"  ,    "case_excluded" , "LL_thresh"   ,   "other_thresh"  , "incl_N" , "necessary", "TT_inclusion_thresh"    ,   
                                    "solutions"  ,    "solution_stats")
  


iteration <- 0

#suppress warnings
defaultW <- getOption("warn") 
options(warn = -1) 


  
  #now rotate through each version of each condition, nested after one another
  LL_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="LL",]
  CI_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="CI",]
  social_protection_loop <- calibration_alternatives[calibration_alternatives$variable=="social_protection",]
  rec_sport_culture_loop <- calibration_alternatives[calibration_alternatives$variable=="rec_sport_culture",]
  I_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="I",]
  A_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="A",]
  CP_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="CP",]
  P_thresholds_loop <- calibration_alternatives[calibration_alternatives$variable=="P",]
  
  iterations_to_run <- nrow(LL_thresholds_loop)*nrow(CI_thresholds_loop)*nrow(social_protection_loop)*nrow(rec_sport_culture_loop)*nrow(I_thresholds_loop)*nrow(A_thresholds_loop)*nrow(CP_thresholds_loop)*nrow(P_thresholds_loop)*2 #final factor 2 is for the 2 different inclusion cutoffs
  cat("\nIterations to run:", iterations_to_run, "\n")
  
  
  #main calibration loop####
  for (ll_this in 1:nrow(LL_thresholds_loop)){
    for (ci_this in 1:nrow(CI_thresholds_loop)){
      for (soc_prot_this in 1:nrow(social_protection_loop)){
        for (rec_sport_culture_this in 1:nrow(rec_sport_culture_loop)){
          for (i_this in 1:nrow(I_thresholds_loop)){
            for (a_this in 1:nrow(A_thresholds_loop)){
              for (cp_this in 1:nrow(CP_thresholds_loop)){
                for (p_this in 1:nrow(P_thresholds_loop)){
           
          
                remove_case <- sample(1:(length(cases)+length(cases)), 1, replace=T )  #50% chance to remove one case or keep all
                raw_loop <- rawdt[-remove_case,] #remove_case was a lottery number to be removed. Only 50% are acctually possible to remove
                
            
                this_threshold_table <- rbind(LL_thresholds_loop[ll_this,], CI_thresholds_loop[ci_this,], social_protection_loop[soc_prot_this,], rec_sport_culture_loop[rec_sport_culture_this,],I_thresholds_loop[i_this,], A_thresholds_loop[a_this,], CP_thresholds_loop[cp_this,], P_thresholds_loop[p_this,])
                
                #now calibrate by combining this with raw_loop
                
                
                calibrated_loop <- data.frame(matrix(nrow=nrow(raw_loop), ncol=9))
                colnames(calibrated_loop) <- colnames(raw_loop)
                calibrated_loop$GEO <- raw_loop$GEO
              
                for (r in 1:nrow(this_threshold_table)){
                  
                  one <- this_threshold_table[r, 'i']
                  pt_five <- this_threshold_table[r, 'c']
                  zero <- this_threshold_table[r, 'e']
                  thresh <- ifelse(this_threshold_table[r,'valence']==1, paste("e=", zero, ", c=", pt_five, ",i=", one), paste("i=", one, ", c=", pt_five, ",e=", zero) )
                  this_variable <- this_threshold_table[r, 'variable']
                  data_to_calibrate <- unlist(unname(raw_loop[,this_variable]))
                  new_data <-  calibrate(data_to_calibrate, thresholds =thresh)
                  
                  calibrated_loop[,this_variable] <- new_data
                  
                 
                  
                }
                
                #calculate now the new PI
                
                calibrated_loop$PI <-  with(calibrated_loop, fuzzyand(social_protection,rec_sport_culture))

                #now remove the components used to build PI
                calibrated_loop <- calibrated_loop[,!colnames(calibrated_loop) %in% c('rec_sport_culture', 'social_protection')]
                
                dt <- calibrated_loop
                
                
                necessity <-superSubset(dt, 
                                         outcome = "LL", 
                                         neg.out = FALSE, 
                                         relation = "necessity", 
                                         conditions = names(dt)[3:8],
                                         incl.cut = 0,  #should be .89, but I manually cut it to ensure no error
                                         use.tilde = FALSE, 
                                         use.letters = FALSE,
                                      depth=1  #this prevents any conjunctions or disjunctions
                                         )
                

                inclusion_necessity <- ifelse(max(necessity$incl.cov$inclN)>=.89, list(necessity$incl.cov$inclN[necessity$incl.cov$inclN>=.89]), NA) #modify the one to put in actual inclusion scores
                if (is.na(inclusion_necessity)){
                  conditions_necessary <- NA
                }else{
                  conditions_necessary <-  list(rownames( necessity$incl.cov)[necessity$incl.cov$inclN>=.89])  
                }
                              
                
                #now TT analysis
                TT  <- truthTable(dt, outcome = "LL", neg.out = FALSE, conditions = names(dt)[3:8], n.cut = 1, incl.cut = 0.75,complete = FALSE, show.cases = TRUE,sort.by = c("incl", "n"), decreasing = TRUE, use.letters = FALSE)
                #now 1. look for contradictory rows
                # 2. look for TT gaps
                #3. compare the .75 cutoff to the other gaps. 
                cases_above_pt75 <- TT$tt$cases[TT$tt$incl>=.75 & TT$tt$PRI>=.50]
                #now convert these to integers... report which ones are not having the outcome, set their rows to 0.

                cases_above_pt75 <-  as.integer( unlist(str_split(cases_above_pt75, ",", n = Inf, simplify = FALSE)))
                which_not_LL <- cases_above_pt75[which(dt$LL[cases_above_pt75]<.5)]
                which_tt_rows_to_zero  <- as.integer( rownames(TT$tt)[TT$tt$cases %in% which_not_LL])
               TT$tt$OUT[which_tt_rows_to_zero] <- 0
               
               #now get the gaps.
               inclusion_vector <- sort(as.numeric(  TT$tt$incl), decreasing = T)
               inclusion_vector <- inclusion_vector[inclusion_vector>=.75]
               gaps <- inclusion_vector[ -length(inclusion_vector)]-inclusion_vector[-1]
               big_gaps <- sort(gaps, decreasing = T)
               gap2 <- which(gaps==big_gaps[1]) #gap1_threshold is .75
               #gap3 <- which(gaps==big_gaps[2]) #gap1_threshold is .75
                gap2_threshold <- round(as.numeric( TT$tt$incl[ TT$rowsorder[gap2]] ), digits=3)-.001 #will be immediately beneath that row number in the displayed ordered TT
                #gap3_threshold <- round(as.numeric( TT$tt$incl[ TT$rowsorder[gap3]] ), digits=3)-.001
                
                #now, go through minimization two times, once with gap1=.75 and once with the other two thresholds      
                # finish with this first gap1. TT already built. Then cycle through 2 . can apply the row cancellation for contradictory cases each time
                
#                inclusion_thresholds <- c(.75, gap2_threshold, gap3_threshold)
                inclusion_thresholds <- c(.75, gap2_threshold)
                
                for (threshold in inclusion_thresholds){
                  
                  iteration <- iteration+1
                  
                  cat("\nIteration", iteration)  
                  
                 
                  
                  TT_loop  <- truthTable(dt, outcome = "LL", neg.out = FALSE, conditions = names(dt)[3:8], n.cut = 1, incl.cut = threshold,complete = FALSE, show.cases = TRUE,sort.by = c("incl", "n"), decreasing = TRUE, use.letters = FALSE)
                  
                  TT_loop$tt$OUT[which_tt_rows_to_zero] <- 0
                  
                  #now run the minimization... include the necessary conditions then directional expectations
                  #record the combinations and their inclusion and coverage.

              
              if (is.na( conditions_necessary)){
                ET_loop <- esa(TT_loop, nec_cond = NULL) 
                
              }else{
                ET_loop <- esa(TT_loop, nec_cond = unlist(conditions_necessary)) 
                
              }
                  
                   
                    intermediate_loop <- try(minimize(ET_loop, details = TRUE, include = "?",show.cases = TRUE, dir.exp = c(I, CI, P, CP, PI, A)))

                    if (is.character( intermediate_loop)){  #if an error popped up...
                      intermediate_solutions <- NA
                      intermediate_solution_stats <- NA
                    }else{
                    
                      
                 
                      
                      intermediate_solutions <- list()
                   intermediate_solution_stats <- list()
                   
                   for (i in 1:length(intermediate_loop$i.sol)){
                   intermediate_solutions <- append(intermediate_solutions, intermediate_loop[["i.sol"]][[i]][["solution"]] )
                   
                   if (length(intermediate_loop[["i.sol"]][[i]][["solution"]] )>1 ){
                     my_length <- length(intermediate_loop[["i.sol"]][[i]][["solution"]] )
                     these_stats <- sapply(intermediate_loop$i.sol[[i]]$IC$individual, "[[", "incl.cov", simplify = F)
                     intermediate_solution_stats <- append(intermediate_solution_stats, these_stats)
                     
                   }else{ #if there is only one solution
                  intermediate_solution_stats <- append(intermediate_solution_stats, list( intermediate_loop$i.sol[[i]][['IC']][['incl.cov']]))
                   }
                   } 
                   
                   which_to_keep <- match(unique(intermediate_solutions),intermediate_solutions )
                   intermediate_solutions <- unique(intermediate_solutions)
                   
                   intermediate_solution_stats <- intermediate_solution_stats[which_to_keep]

                   }
                   
                 

          
                  grand_results_table_new <- data.table(iteration=iteration,case_excluded=cases[remove_case], LL_thresh=list(unname(this_threshold_table[1,3:5])), other_thresh=list(unname(this_threshold_table[2:8,3:5])) ,  incl_N=inclusion_necessity, necessary=conditions_necessary,TT_inclusion_thresh=list( inclusion_thresholds), solutions=intermediate_solutions, solution_stats=intermediate_solution_stats)
                  
               grand_results_table <- rbind(grand_results_table,grand_results_table_new )  
               
               #reset those lists
               intermediate_solutions <- intermediate_solutions[0]
               intermediate_solution_stats <- intermediate_solution_stats[0]
               
               if (iteration %% 2000==0){ #change to 2000. That is just how often the results will be saved
                 
               

                 saveRDS(grand_results_table, paste0(iteration," logged robustness checks.RDS"))
               
               grand_results_table <- data.table(matrix(nrow=0, ncol=9))
               colnames(grand_results_table) <- c( "iteration"  ,    "case_excluded" , "LL_thresh"   ,   "other_thresh"  , "incl_N" , "necessary", "TT_inclusion_thresh"    ,   
                                                   "solutions"  ,    "solution_stats")    

                     }
                }
                
                          }
              
              
          }
          }
      }
      }
      }
    }
    
  
  }

#after main loop####
#options(warn = defaultW)




saveRDS(grand_results_table, paste0(iteration,"new logged robustness checks.RDS"))

  
found <-list.files(pattern = "logged robustness",  full.names=TRUE) 

robust_dt<-NULL
for (i in found){
  robust_dt<-rbind(robust_dt,data.table(readRDS(paste0(i)) )  ) 
}

saveRDS(robust_dt, "main robustness table.RDS")
#Load data####
#robust_dt <-data.frame(  readRDS("main robustness table new debugged.RDS"))

robust_dt$necessary <- replace(robust_dt$necessary, is.na(robust_dt$necessary), 0)
robust_dt$necessary  <- lapply(robust_dt$necessary, FUN=function(x)  paste(x, collapse = "*"))
robust_dt$necessary <- factor(unlist(robust_dt$necessary))

#and how make it summarized by iteration! take only the first iteration
mini_robust_dt <- robust_dt[match(1:max(robust_dt$iteration),robust_dt$iteration ),]

#and now I need to take only the ODD lines, because the inclusion threshold variation is irrelevant for necessity testing
mini_robust_dt <- mini_robust_dt[as.logical(mini_robust_dt$iteration %% 2),] 

necessity_table <- data.frame(table(mini_robust_dt$necessary))

necessity_table <- necessity_table[order(necessity_table$Freq, decreasing = T),]


#TABLE. Necessary conditions robustness table####

#we do not want this done for model ambiguity... because the necessity analysis is independent of that!

necessity_table$percent <- round(necessity_table$Freq/sum(necessity_table$Freq), digits = 3)*100
  
#write.csv(necessity_table, "necessity.csv")

#also check if necessary conditions are found within individual solution terms####

#Then sufficiency
check_sufficiency_robustness <- function(solutions=robust_dt$solutions, 
                                         solution_stats=robust_dt$solution_stats, 
                                         combination=testing_these_solution_paths,
                                         calc_stats=F){


combination_chunks <- unlist(str_split(combination, " \\+ "))
original_combination_chunks <- combination_chunks

#negations become X. There was a problem with the tilde

combination_chunks <- gsub("\\~", "X", combination_chunks)
  
#exact match. and subsets found#
#add another column for direct presence and another for subsets of this presence. Calculate stats only for direct presence.
#FIX this to add also subsets of the chunk

#now replace solution term negations .. the tilde with an X
solutions <- lapply(solutions, FUN=function(x) gsub("\\~", "X", x))

solutions_atomized <- lapply(solutions, FUN=function(x) str_split(x, "\\*"))

combination_chunks_atomized <-str_split(combination_chunks, "\\*")


is_this_combo_found <- function(solution_list_mini){
  
  lapply(solution_list_mini, FUN=function(x) lapply(combination_chunks_atomized, FUN=function(y) as.logical(prod(y %in% x ))))
}


  is_it_here <- lapply(solutions_atomized, FUN=is_this_combo_found)

#relists these results by combination chunk, showing true or false applied to each solution term for each chunk 
by_combination <-   lapply(is_it_here, FUN=function(x) split(unlist(x), rep(1:length(combination_chunks), times=length(x))))

#for each iteration, each list object contains a vector of whether or not the corresponding combination chunk is found as a solution term or as a subset of the solution term
by_iteration <- lapply(by_combination, FUN=function(x) unlist(lapply(x, FUN=function(x) T %in% x)))

variable_vectors <- split(unlist(by_iteration), f=names(by_iteration[[1]]))

robust_dt <- data.table(robust_dt)
set(robust_dt, j=combination_chunks, value=as.logical(NA))


for (v in 1:length(variable_vectors)){
  set(robust_dt, j=combination_chunks[v], value=variable_vectors[[v]])
  
}
robust_dt <- data.frame(robust_dt)



#now pull solution stats for exact matches only. 


exact_is_it_here <- lapply(solutions, FUN=function(x) combination_chunks %in% x)

  
  if (length( unlist(exact_is_it_here))>0){
  
    #first get the row indices!
    #make my own loop. mapply inappropriate.
    suff_stats_list <- vector("list", length = 0)
    for (r in 1:length(solution_stats)){
  
which_stats <-       solution_stats[[r]][exact_is_it_here[[r]],]
suff_stats_list <- append(suff_stats_list, list(which_stats))
    }
    

    
    #now add these data to robust_dt
    combination_variables <- paste0(rep(combination_chunks, each=3), c(".consistency", ".PRI", ".coverage"))
    combination_variables <- gsub("\\*", ".",combination_variables )
    #creates these variables
    robust_dt <- data.table(robust_dt)
    
    set(robust_dt, j=combination_variables, value=as.numeric(NA))
    robust_dt <- data.frame(robust_dt)
    
    options(warn = 2)    
    
    if (calc_stats==T){
    #now use is_it_here to create a vector of values for each row in robust_dt, mapply to columns is 9-11, 12-14, 15-17, etc..
    for (suff in 1:length(suff_stats_list)){
   
      
       x <-  suff_stats_list[[suff]] 
     x <-   x[rownames(x) %in% original_combination_chunks,]
     
     if (!is.null(x) & !nrow(x)==0){
 

       
       
sufficiency_value_vector <-      unname(as.vector(t(x[,1:3] )))
#column_indices <- c(12+(1:(length(combination_chunks)*3)))*rep(exact_is_it_here[[s]]*1, each =3)

correct_variable_names <- combination_variables[   as.logical( c(rep(1, length(combination_chunks)*3)*rep(exact_is_it_here[[suff]]*1, each =3)))]
column_indices <- which(colnames(robust_dt) %in% correct_variable_names)


 robust_dt[suff,column_indices ] <- sufficiency_value_vector
     }

    }
    
    } #end if for calc_stats

  
  }
  
  return(list( dt=robust_dt, chunks=combination_chunks, variables=combination_variables))
 
}

results <- check_sufficiency_robustness()
robust_dt <- results[[1]]
combination_chunks <- results[[2]]
combination_variables <- results[[3]]


nrow(robust_dt)
nrow(distinct(robust_dt))

#TABLES. exact terms or their subsets in ALL solutions####
for (combo in 1:length(combination_chunks)){
  #FIX!
  converted_chunks <- gsub("\\*", "\\.", combination_chunks)

this_table <-   data.frame(table(robust_dt[, colnames(robust_dt)== converted_chunks [combo]]))
this_table <- rename(this_table, present=Var1)  
cat( "\n\n",combination_chunks[combo], round(this_table$Freq[this_table$present==TRUE]/sum(this_table$Freq), digits=3)*100,  "%\n")
print(this_table)

#get for each combo the right columns
correct_variable_names <- combination_variables[3*(combo-1)+c(1:3)  ]

these_columns <- which(colnames(robust_dt) %in% correct_variable_names)

par(mfrow=n2mfrow(3, asp = 1)) #the three corresponds to the number of stats displayed
for (s in these_columns){
if (!is.null(length( robust_dt[,s][!is.na(robust_dt[,s])]))){

  hist(robust_dt[,s][!is.na(robust_dt[,s])], breaks=20, main=combination_chunks[combo], xlab=colnames(robust_dt)[s])
}
  }

}



#TABLES. by robustness choice combinations####
#now frequencies based on unique robustness parameter combinations####
these_combo_columns <- which(colnames(robust_dt) %in% gsub("\\*", ".", combination_chunks) )

#now running these by ITERATION of robustness parameters, not by solution
for (cc in 1:length( these_combo_columns)){
xx <-   aggregate(robust_dt[,these_combo_columns[cc]], by=list(iteration=robust_dt$iteration),  FUN=function(x) T %in% x  )
xx2 <-data.frame(  table(xx$x))
colnames(xx2)[1] <- "present"
cat( "\n\n",combination_chunks[cc], round(xx2$Freq[this_table$present==TRUE]/sum(xx2$Freq), digits=2)*100,  "%\n")
print(xx2)
}

#TABLE.Case exclusion table####
these_rows <- match(unique(robust_dt$iteration), robust_dt$iteration)
case_table <-data.frame( table(robust_dt$case_excluded[these_rows], useNA = "ifany"))
print(case_table, row.names = F)

saveRDS(robust_dt, "robustness results.RDS")

#necessity robustness####
dt <- readRDS("main robustness table newest.RDS")


#now frequencies based on unique robustness parameter combinations#
#"How many parameter combinations result in at least one model where the necessary conditions are present in all paths?##

necessity_chunks <- unlist(str_split(necessary_conditions, " \\+ "))

#negations become X

necessity_chunks <- gsub("\\~", "X", necessity_chunks)

#now cycle through all iterations in dt... searching for one solution that matches the necessity claim.

necessity_vector <- vector("logical", max(dt$iteration))

for (i in sort(unique(dt$iteration))){
  these_rows <- which(dt$iteration==i)
  
  this_iteration_has_nec <- FALSE
  
  for(i2 in these_rows){
   this_solution <- unlist( dt$solutions[i2])
   
   #now replace solution term negations .. the tilde with an X
   this_solution <- gsub("\\~", "X", this_solution)
   
   solutions_atomized <- lapply(this_solution, FUN=function(x) str_split(x, "\\*"))
   
each_path_has_nec <-    lapply(solutions_atomized, FUN=function(x)  as.logical(prod(necessity_chunks %in% unlist( x) )))
   
if (prod(unlist(each_path_has_nec))==1){
  this_iteration_has_nec <- TRUE
  break
}

  }
  
  necessity_vector[i] <- this_iteration_has_nec
}

data.frame(table(necessity_vector))

#now running these by ITERATION of robustness parameters, not by solution
for (cc in 1:length( these_combo_columns)){
  xx <-   aggregate(robust_dt[,these_combo_columns[cc]], by=list(iteration=robust_dt$iteration),  FUN=function(x) T %in% x  )
  xx2 <-data.frame(  table(xx$x))
  colnames(xx2)[1] <- "present"
  cat( "\n\n",combination_chunks[cc], round(xx2$Freq[this_table$present==TRUE]/sum(xx2$Freq), digits=2)*100,  "%\n")
  print(xx2)
}





#TABLES. crucial pairs####
critical_pairs <-   "CI*~PI + I*A + P*PI + CP*P"
# critical_pairs <-   "~CI*PI"
critical_pairs <-   "I + A + CI + ~PI + PI + CP + P"

results <- check_sufficiency_robustness(combination = critical_pairs)
robust_dt <- results[[1]]
combination_chunks <- results[[2]]
combination_variables <- results[[3]]

robust_dt <- data.frame(robust_dt)

for (combo in 1:length(combination_chunks)){

  converted_chunks <- gsub("\\*", "\\.", combination_chunks)
  
  this_table <-   data.frame(table(robust_dt[, colnames(robust_dt)== converted_chunks [combo]]))
  this_table <- rename(this_table, present=Var1)  
  cat( "\n\n",combination_chunks[combo], round(this_table$Freq[this_table$present==TRUE]/sum(this_table$Freq), digits=3)*100,  "%\n")
  print(this_table)
  
  #get for each combo the right columns
  correct_variable_names <- combination_variables[3*(combo-1)+c(1:3)  ]
  
  these_columns <- which(colnames(robust_dt) %in% correct_variable_names)
  

  
}


# #now running these by ITERATION of robustness parameters, not by solution
# for (cc in 1:length( these_combo_columns)){
#   xx <-   aggregate(robust_dt[,these_combo_columns[cc]], by=list(iteration=robust_dt$iteration),  FUN=function(x) T %in% x  )
#   xx2 <-data.frame(  table(xx$x))
#   colnames(xx2)[1] <- "present"
#   cat( "\n\n",combination_chunks[cc], round(xx2$Freq[this_table$present==TRUE]/sum(xx2$Freq), digits=2)*100,  "%\n")
#   print(xx2)
# }
