rm(list = ls())
f=read.table("AIMS_selection/from_suyash/fig1c_AFR_vs_EUR.txt",header=T)
#head(f)
ob=lm( Estimate ~ Truth,f)
print(ob)
print(summary(ob))
print(attributes(ob))

library(ggplot2)

# jpeg("fig1c_AFR_vs_EUR.png")
# with(f,plot(Truth,Estimate,xlim=c(0,1),ylim=c(0,1)))
# abline(b=1,a=0,col="red")
# abline(ob)
# legend("topleft",legend=c("Ideal fit","Regression fit"),col=c("red","black"),lwd=1)
# dev.off()

rsquared=round(cor(f$Truth,f$Estimate)^2,digits=2)
x="r"
rsquaredstr=expression(paste(x^2,"=",rsquared))

p<-ggplot(f,aes(x=Truth,y=Estimate))+geom_point()+geom_smooth(method=lm,se=T)
p<-p+geom_text(data=NULL,x=0.125,y=1.0,size=7,label=paste("R^2==",rsquared),parse=T)#+geom_abline(intercept=0,slope=1,colour="red")
p<-p+theme_bw()+xlim(0,1)+ylim(0,1)
p<-p+theme(axis.title.x = element_text(size=20),
	axis.text.x  = element_text(size=20,angle=0,colour='black'))
p<-p+theme(axis.title.y = element_text(size=20),
	axis.text.y  = element_text(size=20,colour='black'))
p<-p+theme(legend.title = element_text(size=20),
	legend.text = element_text(size = 20))
p<-p+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),panel.background=element_blank()) + theme(strip.background=element_blank(), strip.text.y = element_blank())
p<-p+xlab("Genome-wide African ancestry")+ylab("Estimated African ancestry")	

ggsave('AIMS_selection/figures/fig1c_AFR_vs_EUR.png',p, width = 6, height = 6)
ggsave('AIMS_selection/figures/fig1c_AFR_vs_EUR.pdf',p, width = 6, height = 6)