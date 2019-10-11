
library(EMreading)

EyeDoctor_PadLines(data_dir = 'D:/Data/JITTER', paddingSize = 5)

raw_fix<- preprocFromDA1(data_dir = 'D:/Data/JITTER', maxtrial = 105, tBlink = 150, padding = 5)

t<- ExtractMessages(data_list = 'D:/Data/JITTER', maxtrial = 105,   message_name =  
                      c('DISPLAY CHANGE STARTED', 'DISPLAY CHANGE COMPLETED'))


