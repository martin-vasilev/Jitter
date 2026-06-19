
##### Old data- add corrective saccade landing:
rm(list= ls())

load("~/R/Jitter/Pilot/data/raw_fix.Rda")

loc<- which(raw_fix$Rtn_sweep_type=='undersweep')

raw_fix$corr_sacc_land<- NA
raw_fix$corr_sacc_land[loc]<- raw_fix$xPos[loc+1]


# Extract messages with line margin:
msg<- ExtractMessages(data_list = 'C:/Users/Martin/Documents/R/Jitter/Pilot/raw', maxtrial = 105,   message_name =  
                      'NEXT OFFSET')

msg$offset<- NA

for(i in 1:nrow(msg)){
  msg$offset[i]<- as.numeric(unlist(strsplit(msg$whole_message[i], ' '))[4])
}

msg<- msg[,c(1:3,8)]

library(tidyverse)
raw_fix<- inner_join(raw_fix, msg, by= c('sub', 'item'))


Old<- subset(raw_fix, Rtn_sweep==1)

Old$Rtnsweep_land_char<- ceiling((Old$xPos - Old$offset)/13)
Old$corr_sacc_land_char<- ceiling((Old$corr_sacc_land - Old$offset)/13)

save(Old, file = 'Pilot/data/Old_RS_data.Rda')


library(EMreading)
#EyeDoctor_PadLines(data_dir = 'D:/Data/JITTER/new', paddingSize = 5)

raw_fix<- preprocFromDA1(data_dir = 'C:/Data/corr_sacc',
                         maxtrial = 90, tBlink = 150, padding = 0)

t<- ExtractMessages(data_list = 'C:/Data/corr_sacc', maxtrial = 105,   message_name =  
                      c('DISPLAY CHANGE STARTED', 'DISPLAY CHANGE COMPLETED'))


# save(t, file= 'Pilot/data/t.Rda')
# write.csv(t, 'Pilot/data/t.csv')
# save(raw_fix, file= 'Pilot/data/raw_fix.Rda')
# write.csv(raw_fix, 'Pilot/data/raw_fix.csv')


#RS_new<- subset(raw_fix, Rtn_sweep==1)
# save(RS_new, file= "Pilot/data/RS_new.Rda")
# write.csv(RS_new, "Pilot/data/RS_new.csv")


loc<- which(raw_fix$Rtn_sweep_type=='undersweep')

raw_fix$corr_sacc_land<- NA
raw_fix$corr_sacc_land[loc]<- raw_fix$xPos[loc+1]

# Extract messages with line margin:
msg<- ExtractMessages(data_list = 'C:/Data/corr_sacc', maxtrial = 90,   message_name =  
                        'NEXT OFFSET')

msg$offset<- NA

for(i in 1:nrow(msg)){
  msg$offset[i]<- as.numeric(unlist(strsplit(msg$whole_message[i], ' '))[4])
}

msg<- msg[,c(1:3,8)]

library(tidyverse)
raw_fix<- inner_join(raw_fix, msg, by= c('sub', 'item'))


New<- subset(raw_fix, Rtn_sweep==1)

New$Rtnsweep_land_char<- ceiling((New$xPos - New$offset)/18)
New$corr_sacc_land_char<- ceiling((New$corr_sacc_land - New$offset)/18)



New<- New[, c('sub', 'item', 'cond.x', 'Rtn_sweep_type', 'char_line', 'Rtnsweep_land_char', 'corr_sacc_land_char')]
Old<- Old[, c('sub', 'item', 'cond.x', 'Rtn_sweep_type', 'char_line', 'Rtnsweep_land_char', 'corr_sacc_land_char')]

New$sample<- "New"
Old$sample<- "Old"

New$sub<- New$sub+4

final<- rbind(Old, New)


final$undersweep_prob<- NA
final$undersweep_prob[which(final$Rtn_sweep_type=='accurate')]<- 0
final$undersweep_prob[which(final$Rtn_sweep_type=='undersweep')]<- 1


write.csv(final, 'Pilot/data/final.csv')


library(tidyverse)

RS$bnd<- NA

for(i in 1:nrow(RS)){
  loc<- which(t$sub== RS$sub[i]& t$item== RS$item[i])
  
  if(length(loc)>0){
    RS$bnd[i]<- t$MSG2[loc]
  }

    
}

RS$bnd_time_fix<- RS$SFIX - (RS$bnd +8)


library(reshape)
DesRS<- melt(RS, id=c('sub', 'item', 'cond'), 
                measure=c("undersweep_prob", 'char_line'), na.rm=TRUE)
m<- cast(DesRS, cond ~ variable
              ,function(x) c(M=signif(mean(x),3)
                             , SD= sd(x) ))


DesFD<- melt(RS, id=c('sub', 'item', 'cond', 'undersweep_prob'), 
             measure=c("fix_dur"), na.rm=TRUE)
m2<- cast(DesFD, cond + undersweep_prob ~ variable
         ,function(x) c(M=signif(mean(x),3)
                        , SD= sd(x) ))

write.csv(m, 'prob.csv')
