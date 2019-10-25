
library(EMreading)

#EyeDoctor_PadLines(data_dir = 'D:/Data/JITTER', paddingSize = 5)

raw_fix<- preprocFromDA1(data_dir = 'D:/Data/JITTER', maxtrial = 105, tBlink = 150, padding = 5)

t<- ExtractMessages(data_list = 'D:/Data/JITTER', maxtrial = 105,   message_name =  
                      c('DISPLAY CHANGE STARTED', 'DISPLAY CHANGE COMPLETED'))


save(t, file= 'Pilot/t.Rda')
save(raw_fix, file= 'Pilot/raw_fix.Rda')

RS<- subset(raw_fix, Rtn_sweep==1)
save(RS, file= "Pilot/RS.Rda")


RS$undersweep_prob<- NA
RS$undersweep_prob[which(RS$Rtn_sweep_type=='accurate')]<- 0
RS$undersweep_prob[which(RS$Rtn_sweep_type=='undersweep')]<- 1


library(reshape)
DesRS<- melt(RS, id=c('sub', 'item', 'cond'), 
                measure=c("undersweep_prob", 'char_line'), na.rm=TRUE)
m<- cast(DesRS, cond ~ variable
              ,function(x) c(M=signif(mean(x),3)
                             , SD= sd(x) ))
