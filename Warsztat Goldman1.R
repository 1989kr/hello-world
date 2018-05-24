setwd("C:/Users/Krzysztof/Desktop/Goldman")
setwd("C:/Users/krroman/Desktop/Goldman")
mydata<-read.table("quantworkshop.csv", sep=";", header = TRUE)
input<-mydata

debug(utils:::unpackPkgZip)
trace(utils:::unpackPkgZip, edit=TRUE)

install.packages("ggplot2")
library(ggplot2)
install.packages("dplyr")
library("dplyr")
install.packages("magrittr")
library(magrittr)

library(ggplot2)
library("dplyr")
library(magrittr)

#"Black Scholes

BSAnalyticalPricer<-function (S,K, Tenor, sigma, r){
relationSK=S/K
d1=(log(relationSK)+(r+sigma^2/2)*Tenor)/(sigma*sqrt(Tenor)) 
d2=d1-sigma*sqrt(Tenor)  
price=S*pnorm(d1)-exp(-r*Tenor)*K*pnorm(d2)
return(price)
}


BSAnalyticalPricer(1,1,1,0.1,0.01)
# S = 1
# K=0
# Tenor=1
# sigma=1
# r=0.1

to_plot=input%>%mutate(Price=BsAnalyticalPricer(S,K, Tenor, sigma, r))
to_plot$Tenor=as.factor(to_plot$Tenor)
to_plot




plot<-ggplot(to_plot,aes(x = K,y =Price, col=Tenor))+geom_line()+ geom_point()+
  xlab("Strike")+ylab("Price")+ ggtitle("Dependence of prices on Tenor and Strike")+
  theme_light()+theme(legend.position = "bottom",legend.title=element_blank())
plot

plot1<-ggplot(to_plot, aes(x=K, y=sigma, col=Tenor))+geom_point()+geom_line()+
  xlab("Strike")+ylab("Implied vol")+ ggtitle(paste("Dependence of sigma on Tenor and Strike"))+
  theme_light()+theme(legend.position = "bottom",legend.title=element_blank())
plot1

#Monte Carlo  

BSPricerMC<-function(S,K,Tenor,sigma,r,repetition){
St=S*exp((r-sigma^2/2)*Tenor+sigma*rnorm(repetition))  
SK=St-K
payoff=pmax(SK,0)
return(exp(-Tenor*r)*mean(payoff))
}

BSPricerMC1<-function(S,K,Tenor,sigma,r,repetition){
  St=S*exp((r-sigma^2/2)*Tenor+sqrt(Tenor)*sigma*rnorm(repetition))  
  tresholding=sapply(St-K,function(x) max(x,0))
  return(exp(-Tenor*r)*mean(tresholding))
}

BSPricerMC(1,1,1,0.1,0.01,1000)
BSPricerMC1(1,1,1,0.1,0.01,1000)


inputMC=input%>% filter(Tenor==2)%>%inner_join(data.frame(Tenor=2,repetitions=seq(1,5,by=0.2)))

y=data.frame(Tenor=2,repetitions=seq(1,5,by=0.2))
x=filter(input,Tenor==2)
inputMC=inner_join(x,y)

to_plotMC=inputMC%>%group_by(K,repetitions)%>%mutate(Price=BSAnalyticalPricer(S,K,Tenor,sigma,r), PriceMC=BSPricerMC1(S,K,Tenor,sigma,r,10^repetitions))

yMC<-group_by(inputMC,K,repetitions)
yMC1<-mutate(yMC,Price=BSAnalyticalPricer(S,K,Tenor,sigma,r), PriceMC=BSPricerMC1(S,K,Tenor,sigma,r,10^repetitions))
  
plotMC<-ggplot(to_plotMC,aes(x=repetitions, y= PriceMC))+geom_point()+
  geom_line(aes(x=repetitions, y=Price))+
  xlab("MC sample size, logscale")+ ylab("Price")+
  facet_wrap(~paste("Strike",K), scale="free")+theme_classic()+theme(text=element_text(size=100))+
  ggtitle(paste("G³upi tytu³"))+theme_light()+theme(legend.position="bottom", legend.title = element_blank())

plotMC


#Implied Vol

ImpliedVol<-function(Price,K2,Tenor2){
res=uniroot(
  function(sigma1) return(BSAnalyticalPricer(1,K2,Tenor2,sigma1,0.01)-Price)
  ,c(0,1), tol=0.0001)  
return(res$root)  
}

sigmain=0.1
Price=BsAnalyticalPricer(1,1,1,sigmain,0.01)
Price
  
sigma=ImpliedVol(Price,1,1)
abs(sigma-sigmain)<.00001
  
#Gamma MonteCarlo

BSPricerMCGamma<-function(S,K,Tenor, sigmamean,r,repetitions,sdgamma){  
  #generate different states of nature
rate=sigmamean/sdgamma^2
shape=sigmamean*rate
vol=rgamma(repetitions, shape=shape, rate=rate)
Price=sapply(vol,function(x) BSAnalyticalPricer(S,K,Tenor,x,r))
return(mean(Price))
}

sigmain=0.1
PriceMCGamma=BSPricerMCGamma(1,1,1,sigmain,0.01,10000,1)
PriceMCGamma

xG<-filter(input,Tenor==1, K==1)
xG1<-data.frame(Tenor=1, K=seq(0.5,1.5, by=0.1))
xG3<-data.frame(Tenor=1,sigmamean=0.1,sdgamma=c(0.1,0.2,0.5)*0,1)
xG2<-select( xG3,-K)
inputMCGamma=input%>%filter(Tenor==1,K==1)%>%inner_join(data.frame(Tenor=1,sigmamean=0.1,sdgamma=c(0.1,0.2,0.5)*0.1))%>% 
  select(-K)%>%inner_join(data.frame(Tenor=1, K=seq(0.5,1.5, by=0.1)))
  
to_plotMCG<-inputMCGamma%>%group_by(Tenor, K, sigmamean,sdgamma)%>%mutate(Price=BSPricerMCGamma(S,K,Tenor, sigmamean, r,10000,sdgamma),PriceBS=BSAnalyticalPricer(S,K,Tenor,sigmamean,r))
to_plotMCG1<-to_plotMCG%>%group_by(Tenor, K, Price, sigmamean,sdgamma
                                   )%>%mutate(ImpliedVol=ImpliedVol(Price, K,Tenor))

to_plotMCG1$sdgamma=as.factor(to_plotMCG1$sdgamma)

plot<-ggplot(to_plotMCG1,aes(x=K, y =ImpliedVol, col=sdgamma))+ geom_point()+geom_line() +geom_line(aes(x=K, y=sigmamean, col=rep("0(BS)",dim(to_plotMCG1)[1])))+ xlab("Strike")+ylab("Implied Vol") +theme_classic()+theme(text=element_text(size=7))+ggtitle(paste("Smile Impliec Vol"))+theme_bw()+theme(legend.position = "bottom", legend.title = element_blank())  
plot

to_plot_input<-input%>%filter(Tenor==1)

plot2<-ggplot(to_plotMCG1,aes(x=K, y =ImpliedVol, col=sdgamma))+ geom_point()+
  geom_line() +
  geom_point(data=to_plot_input,aes(x=K, y=sigma, col=rep("Implied vol observed",dim(to_plot_input)[1])))+ 
  xlab("Strike")+ylab("Implied Vol") +theme_classic()+theme(text=element_text(size=7))+
  ggtitle(paste("Smile Impliec Vol"))+theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank())  
plot2

  dim(to_plot_input)[1]
  dim(to_plot_input)[1]

#               col=rep("0(BS)",dim(to_plotMCG1)[1])))
#x=K, y=sigmamean, 
#rep("0(BS)",dim(to_plotMCG1)[1])



"plot different gamma parametrizations