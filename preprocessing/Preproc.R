
rm(list= ls())


# if('devtools' %in% rownames(installed.packages())==FALSE){
#   install.packages('devtools')
#   library(devtools)
# }else{
#   library(devtools)
# }
# install_github('martin-vasilev/EMreading')

library(EMreading)


##Question accuracy:
quest<- Question(data_list = 'C:/Data/corr_sacc/padded',
                maxtrial = 90)

write.csv(quest, file= 'data/quest.csv')


trial_time<- trialTime(data_list = 'C:/Data/corr_sacc/padded',
                       maxtrial = 90)

write.csv(trial_time, file = 'data/trial_time.csv')


# Get fixation data from raw file and combine with Eyedoctor:

raw_fix<- preprocFromDA1(data_dir = 'C:/Data/corr_sacc/asc',
                      maxtrial = 90, tBlink = 150, padding = 5)

write.csv(raw_fix, file = 'data/raw_fixations.csv')

# DC<- ExtractMessages(data_list = 'C:/Data/corr_sacc/padded', 
#                      maxtrial = 90,   message_name =  
#                       c('DISPLAY CHANGE STARTED',
#                         'DISPLAY CHANGE COMPLETED'))
# 
# #DC$flag_time<- NA
# DC$flag<- NA
# 
# for(i in 1:nrow(DC)){
#   
#   msg<- unlist(strsplit(DC$whole_message[i], ' '))
#   time<- msg[1]
#   change<- msg[4] 
#   # DC$flag_time[i]<- as.numeric(unlist(strsplit(time, '\t'))[2])
#    DC$flag[i]<- change
# }
# 
# DC$message<- NULL
# DC$whole_message<- NULL


DC<- Boundary(data_list = 'C:/Data/corr_sacc/padded', boundary_loc = 618, maxtrial = 90)

# save display change data frame:
write.csv(DC, file = 'data/display_changes.csv')


##############################
# Preprocessing of raw data: #
##############################

# See pre-registered protocol for data preprocessing criteria:
# https://osf.io/9sngw#analysis-plan.data-exclusion

#load("preproc/raw_fix.Rda")

raw_fix <- read.csv("~/R/Jitter/data/raw_fixations.csv")

raw_fix$hasText<- NULL


#######################################
# first, let's code some new variables:

raw_fix_new<- NULL

raw_fix$prev_RS<- NA
raw_fix$next_RS<- NA
raw_fix$prevChar<-NA
raw_fix$nextChar<- NA
raw_fix$prevX<- NA
raw_fix$nextX<- NA
raw_fix$prevY<- NA
raw_fix$prev_max_char_line<- NA
raw_fix$prevEFIX<- NA

nsubs<- unique(raw_fix$sub)

for(i in 1:length(nsubs)){
  n<- subset(raw_fix, sub==nsubs[i])
  nitems<- unique(n$item)
  cat(i); cat(" ")
  
  for(j in 1:length(nitems)){
    m<- subset(n, item== nitems[j])
    
    l1<- subset(m, line==1)
    max_l1<- l1$max_char_line[1]
    
    for(k in 1:nrow(m)){
      if(k==1){
        m$prev_RS[k]<- 0
        m$next_RS[k]<- 0
        m$next_RS[k+1]<- 0
        
        ####
        m$nextChar[k]<- m$char_line[k+1] # next char
        m$nextX[k]<- m$xPos[k+1]
        
      }else{
        
        if(is.na(m$SFIX[k])){
          next
        }
        
        if(m$Rtn_sweep[k]==1){
          m$prev_RS[k-1]<- 1
          
          if(k+1 <= nrow(m)){
            m$next_RS[k+1]<- 1
          }
          
        }else{
          m$prev_RS[k-1]<- 0
          
          if(k+1 <= nrow(m)){
            m$next_RS[k+1]<- 0
          }
        }
        ###
        m$prevChar[k]<- m$char_line[k-1] # prev char
        m$prevX[k] <- m$xPos[k-1] # prev x
        m$prevY[k]<- m$yPos[k-1]
        m$prevEFIX[k]<- m$EFIX[k-1]
        
        if(k+1<= nrow(m)){
          m$nextChar[k]<- m$char_line[k+1] # next char
          m$nextX[k]<- m$xPos[k+1] # next x
        }
        
        
      }
      
      if(k== nrow(m)){
        m$prev_RS[k]<- 0
      }
      
      ## map previous line length (for launch site calculation):
      if(!is.na(m$line[k])){
        if(m$line[k]==2){
          m$prev_max_char_line[k]<- max_l1
        }else{
          m$prev_max_char_line[k]<- NA
        }
      }else{
        if(m$Rtn_sweep[k]==1){
          m$prev_max_char_line[k]<- max_l1
        }
      }
      
    } # end of m
    raw_fix_new<- rbind(raw_fix_new, m)
  } # end of j
  
  
}

raw_fix<- raw_fix_new;
rm(raw_fix_new)



nAllTrials<- length(nsubs)*90

########################################
# check number of trials per subject
nTrials<- NULL

for(i in 1:length(nsubs)){
  n<- subset(raw_fix, sub== nsubs[i])
  nTrials[i]<- length(unique(n$item))
}
nTrials

# 3 trials were discarded during manual processing (due to track losses, etc.)

nDiscardedTrials<- nAllTrials- sum(nTrials)


##############################
# check if there is only 1 return sweep per trial (expected, due to the two lines):

noRS_sub<- NULL # no return sweeps this trial
noRS_item<- NULL # no return sweeps this trial
for(i in 1:length(nsubs)){
  n<- subset(raw_fix, sub==nsubs[i])
  nitems<- unique(n$item)
  
  for(j in 1:length(nitems)){
    m<- subset(n, item== nitems[j] & !is.na(SFIX))
    
    if(sum(m$Rtn_sweep)>1){ # >1 RS detected..
      cat(sprintf("Subject %g, item %g, %g return sweeps\n", i, nitems[j], sum(m$Rtn_sweep)))
    }
    
    if(sum(m$Rtn_sweep)==0){ # no return sweeps detected
      cat(sprintf("Subject %g, item %g, %g return sweeps\n", i, nitems[j], sum(m$Rtn_sweep)))
      noRS_sub<- c(noRS_sub, i)
      noRS_item<- c(noRS_item, nitems[j])
      
    }
  }
}


if(length(noRS_sub)>0){
  
  # remove trials with no return sweeps:
  for(i in 1:length(noRS_sub)){
    out<- which(raw_fix$sub== noRS_sub[i]& raw_fix$item== noRS_item[i])
    raw_fix<- raw_fix[-out,]
    
  }
  
}


nNoRS<- length(noRS_sub)

#########################
# let's merge fixations smaller than 80 ms

raw_fix_new<- cleanData(raw_fix, removeOutsideText = F, removeBlinks = F, combineNearbySmallFix = T, 
                        combineMethod = "char", combineDist = 1, removeSmallFix = F, 
                        removeOutliers = F, keepRS = T)
raw_fix<- raw_fix_new
rm(raw_fix_new)

less80<- raw_fix[which(raw_fix$fix_dur<80 & raw_fix$Rtn_sweep==1), ]


# remove remaining fixations less than 80 only if they are not return sweeps:
#raw_fix<- raw_fix[which(raw_fix$fix_dur<80 & raw_fix$Rtn_sweep==1), ]
raw_fix<- raw_fix[-which(raw_fix$fix_dur<80), ]

# Some of the discarded <80 ms fixations may be return sweep ones. Therefore, we will remap fixations for such trials,
# by taking the next fixation as the return sweep one

nsubs<- unique(raw_fix$sub)

new<- NULL
for(i in 1:length(nsubs)){ # for each subject..
  n<- subset(raw_fix, sub== nsubs[i])
  nitems<- unique(n$item)
  
  for(j in 1:length(nitems)){ # for each item..
    m<- subset(n, item== nitems[j])
    
    if(sum(m$Rtn_sweep[which(!is.na(m$SFIX))])>0){
      new<- rbind(new, m)
    }else{
      line= 1
      for(k in 1:nrow(m)){
        if(!is.na(m$line[k])){
          if(m$line[k]> line){
            line= line+1
            m$Rtn_sweep[k]<- 1
            if(m$xPos[k+1]> m$xPos[k]){
              m$Rtn_sweep_type[k]<- "accurate" 
            }else{
              m$Rtn_sweep_type[k]<- "undersweep"
            }
            
            a<- which(less80$sub== m$sub[1] & less80$item== m$item[1])
            if(length(a)>0){
              if(!is.na(less80$char_line[a])){
                m$prevChar[k]<- less80$prevChar[a]
                m$prevX[k]<- less80$prevX[a]
                m$prevY[k]<- less80$prevY[a] 
              }else{
                m$prevChar[k]<- m$char_trial[k-1]
                m$prevX[k]<- m$xPos[k-1]
                m$prevY[k]<- m$yPos[k-1] 
              }
              
            }else{
              m$prevChar[k]<- m$char_trial[k-1]
              m$prevX[k]<- m$xPos[k-1]
              m$prevY[k]<- m$yPos[k-1] 
            }
            
          }
        }
      } 
      new<- rbind(new, m)
      cat(sprintf('changed sub %g item %g \n', m$sub[1], m$item[1]))
    }
    
  }
  cat(i); cat("\n")
}

raw_fix<- new; rm(new)


### double-check to make sure there are no trials with 0 RS..
noRS_sub<- NULL # no return sweeps this trial
noRS_item<- NULL # no return sweeps this trial
for(i in 1:length(nsubs)){
  n<- subset(raw_fix, sub==nsubs[i])
  nitems<- unique(n$item)
  
  for(j in 1:length(nitems)){
    m<- subset(n, item== nitems[j]& !is.na(SFIX))
    
    if(sum(m$Rtn_sweep)>1){ # >1 RS detected..
      cat(sprintf("Subject %g, item %g, %g return sweeps\n", i, nitems[j], sum(m$Rtn_sweep)))
    }
    
    if(sum(m$Rtn_sweep)==0){ # no return sweeps detected
      cat(sprintf("Subject %g, item %g, %g return sweeps\n", i, nitems[j], sum(m$Rtn_sweep)))
      noRS_sub<- c(noRS_sub, i)
      noRS_item<- c(noRS_item, nitems[j])
      
    }
  }
}

# # one more trial is removed as there are no more fixations on the second line
# for(i in 1:length(noRS_sub)){
#   out<- which(raw_fix$sub== noRS_sub[i]& raw_fix$item== noRS_item[i])
#   raw_fix<- raw_fix[-out,]
#   
# }
# nNoRS<- nNoRS+ length(noRS_item)


##########################
# check for blinks occuring on return sweep saccade/fixation:
RS_blinks<- raw_fix[which(raw_fix$Rtn_sweep==1), ]
RS_blinks<- subset(RS_blinks, blink==1)

nBlinks<- nrow(RS_blinks)

if(nBlinks> 0){
  for(i in 1:nrow(RS_blinks)){
    raw_fix<- raw_fix[-which(raw_fix$sub== RS_blinks$sub[i]& raw_fix$item== RS_blinks$item[i]),]
  }
  
}


# remove also blinks that did not occur next to return sweeps:
raw_fix<- raw_fix[-which(raw_fix$blink==1),]


# confirm we got the correct num of trials:
RS_blinks<- raw_fix[which(raw_fix$Rtn_sweep==1), ]
nrow(RS_blinks)== nAllTrials- nBlinks- nNoRS- nDiscardedTrials

RS_blinks<- subset(RS_blinks, blink==1)

if(nrow(RS_blinks)==0){
  cat("GOOD!")
  rm(RS_blinks)
}else{
  cat(":(")
}

table(raw_fix$blink)

# remove blink columns
raw_fix$blink<- NULL
raw_fix$prev_blink<- NULL
raw_fix$after_blink<- NULL


####################
# check for outliers:

outliers<- raw_fix[which(raw_fix$fix_dur>1000),]

outliers<- subset(outliers, Rtn_sweep==1 | prev_RS==1 | next_RS==1)

nOutliers<- nrow(outliers)

if(nOutliers>0){
  
  for(i in 1:nrow(outliers)){
    raw_fix<- raw_fix[-which(raw_fix$sub== outliers$sub[i] & raw_fix$item== outliers$item[i]), ]
  }
  
}


remove_blink<- which(raw_fix$fix_dur>1000)

if(length(remove_blink)>0){
  # remove remaining outlier fixations (not next to return sweeps)
  raw_fix<- raw_fix[-which(raw_fix$fix_dur>1000),]
  
}

outOfBnds<- which(raw_fix$outOfBnds==1 & raw_fix$Rtn_sweep==0 & raw_fix$prev_RS==0 & raw_fix$next_RS==0)

if(length(outOfBnds)){
  # remove fixations outside screen bounds:
  raw_fix<- raw_fix[-outOfBnds,]
}


# not really sure why we still have these, but we remove them here:
outsideLeft<- which(raw_fix$prevX<0 |raw_fix$nextX<0)

if(length(outsideLeft)){
  raw_fix<- raw_fix[-outsideLeft,]
}


###################################################
# let's verify we have the correct number of trials:
RS<- subset(raw_fix, Rtn_sweep==1)

nrow(RS)+ nOutliers+ nBlinks+ nNoRS+ nDiscardedTrials == nAllTrials


### Check display changes:

library(readr)
library(tidyverse)
DC <- read_csv("data/display_changes.csv")

RS<- RS %>%
  inner_join(DC, by= c('sub', 'item', 'cond'))

# check display change occured within return-sweep:
RS$DS_in_RS<- ifelse(RS$tStarted>= RS$prevEFIX & RS$tCompleted<= (RS$SFIX+10),1 ,0)





# now let's print a summary:

fileConn<- file("preprocessing/Preproc_summary.txt", "w")

writeLines(sprintf("- %#.2f percent of trials manually discarded due to tracking loss, etc.",  
                   round((nDiscardedTrials/nAllTrials)*100, 2)), fileConn)

writeLines(sprintf("- %#.2f percent of trials discarded due to the lack of a return sweep in the trial",  
                   round((nNoRS/nAllTrials)*100, 2)), fileConn)

writeLines(sprintf("- %#.2f percent of trials discarded due to blinks on or around return sweeps",  
                   round((nBlinks/nAllTrials)*100, 2)), fileConn)

writeLines(sprintf("- %#.2f percent of trials discarded due to outliers",  
                   round((nOutliers/nAllTrials)*100, 2)), fileConn)

writeLines(sprintf("- %#.2f percent of trials remaining for analysis",  
                   round((nrow(RS)/nAllTrials)*100, 2)), fileConn)

close(fileConn)


###################
# re-organise file:
raw_fix$sent<- NULL # all items are 1 sentence
raw_fix$outOfBnds<- NULL # remove above
raw_fix$time_since_start<- NULL # no longer necessary
raw_fix$prev_RS<- NULL
raw_fix$next_RS<- NULL





