#!/usr/bin/env Rscript

#frequency to count

args = commandArgs(trailingOnly=TRUE)

file = args[1]

freq_table<-read.csv(file, sep = '', header = FALSE, col.names = c('frequency','length'), colClasses = c("numeric","numeric"))
freq_table<-freq_table[freq_table$length>-1 & freq_table$length < 501,]
freq_table<-freq_table[complete.cases(freq_table$frequency),]

freq_table$percent<-NA
for(i in 1:nrow(freq_table)){ 
  freq_table$percent[i] <- 100 * (freq_table$frequency[i] / sum(freq_table$frequency))
}

plot<-ggplot(freq_table, aes( x=as.numeric(as.character(freq_table$length)), y = as.numeric(as.character(freq_table$percent)))) +
  geom_line() + 
  xlim(0,500) + 
  ylim(0,2.5) + 
  scale_color_manual(values=c("#ED6464", "#05188B"), labels = c("CMV", "Human")) +
  geom_vline(xintercept = 167) + 
  theme(legend.position = c(0.8, 0.6)) + 
  ylab('percent')+
  xlab('Fragment size')+
  theme_classic() +
  theme(text = element_text(size=8)) +
  theme(legend.position = c(0.8, 0.6)) + 
  theme(legend.title=element_blank())
plot
ggsave(plot = plot, 'fragment_length.pdf')