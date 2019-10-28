
library(EMreading)

#EyeDoctor_PadLines(data_dir = 'D:/Data/JITTER/new', paddingSize = 5)

raw_fix<- preprocFromDA1(data_dir = 'D:/Data/JITTER', maxtrial = 105, tBlink = 150, padding = 5)

t<- ExtractMessages(data_list = 'D:/Data/JITTER', maxtrial = 105,   message_name =  
                      c('DISPLAY CHANGE STARTED', 'DISPLAY CHANGE COMPLETED'))


save(t, file= 'Pilot/data/t.Rda')
write.csv(t, 'Pilot/data/t.csv')
save(raw_fix, file= 'Pilot/data/raw_fix.Rda')
write.csv(raw_fix, 'Pilot/data/raw_fix.csv')


RS<- subset(raw_fix, Rtn_sweep==1)
save(RS, file= "Pilot/data/RS.Rda")
write.csv(RS, "Pilot/data/RS.csv")

RS$undersweep_prob<- NA
RS$undersweep_prob[which(RS$Rtn_sweep_type=='accurate')]<- 0
RS$undersweep_prob[which(RS$Rtn_sweep_type=='undersweep')]<- 1


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
