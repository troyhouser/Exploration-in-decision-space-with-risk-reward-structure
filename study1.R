rr = read.csv("~/Dropbox (University of Oregon)/risk_reward/study1_learning_choice_merged.csv")
#rr = read.csv("~/Dropbox (University of Oregon)/risk_reward/exp1_bdm_data_merged.csv")
kf = function(sub){
  df = rr[rr$subjid==sub,]
  df$qA=df$qB=df$sA=df$sB=0
  Q=S=rep(0,2)
  gamma=0.85
  at=bt=0
  for(i in 2:nrow(df)){
    muA = mean(df[1:i,"rewA"],na.rm=T)
    muB = mean(df[1:i,"rewB"],na.rm=T)
    tA = sd(df[1:i,"rewA"],na.rm=T)
    tB = sd(df[1:i,"rewB"],na.rm=T)
    if(df$choice[i] == "A"){
      at=at+1
      S[1] = abs(df$rewA[i] - muA) / at
      alpha = S[1]/(S[1] + tA)
      Q[1] = Q[1] + alpha*(df$rewA[i]-Q[1])
      df$qA[i] = Q[1]*df$probA[i]
      df$sA[i] = S[1]
      S[1] = S[1] - alpha*S[1]
      Q[2] = Q[2]*df$probB[i]
      df$qB[i] = Q[2];df$sB[i]=S[2]
    }else{
      bt=bt+1
      S[2] = abs(df$rewB[i] - muB) / bt
      alpha = S[2]/(S[2] + tB)
      Q[2] = Q[2] + alpha*(df$rewB[i]-Q[2])
      df$qB[i] = Q[2]*df$probB[i]
      df$sB[i] = S[2]
      S[2] = S[2] - alpha*S[2]
      Q[1] = Q[1]*df$probA[i]
      df$qA[i] = Q[1];df$sA[i]=S[1]
    }
  }
  df2 = data.frame(subj = df$subjid,
                   qA = df$qA,
                   qB = df$qB,
                   sA = df$sA,
                   sB = df$sB,
                   condition = df$condition,
                   choice = df$evHigherChosen)
  return(df2)
}

subs = unique(rr$subjid)

mfull = list()
m_woutEVD = list()
m_woutRU = list()
m_woutTS = list()
for(s in 1:length(subs)){
  res = kf(subs[s])
  res$evd = res$qA-res$qB
  res$ru = res$sA-res$sB
  res$tu = sqrt(res$sA+res$sB)
  res$ts = res$evd / res$tu
  mfull[[s]] = glm(choice~evd+ru+ts,family=binomial(link = "probit"),res)
  m_woutEVD[[s]] = glm(choice~ru+ts,family=binomial(link = "probit"),res)
  m_woutRU[[s]] = glm(choice~evd+ts,family=binomial(link = "probit"),res)
  m_woutTS[[s]] = glm(choice~evd+ru,family=binomial(link = "probit"),res)
}

bicmat = matrix(0,length(subs),4)
for(i in 1:length(subs)){
  bicmat[i,1] = BIC(mfull[[i]])
  bicmat[i,2] = BIC(m_woutEVD[[i]])
  bicmat[i,3] = BIC(m_woutRU[[i]])
  bicmat[i,4] = BIC(m_woutTS[[i]])
}
library(bmsR)
library(qpcR)
conds = aggregate(condition~subjid,rr,mean)
conds=conds$condition
bics_1 = bicmat[conds==1,]
bics_2 = bicmat[conds==2,]
m = akaike.weights(bics_1)$weights
m = matrix(m,nrow(bics_1),4,byrow=F)
pxp = VB_bms(log(m))
barplot(pxp$pxp)
m = akaike.weights(bics_2)$weights
m = matrix(m,nrow(bics_2),4,byrow=F)
pxp = VB_bms(log(m))
barplot(pxp$pxp)


mfull2 = list()
m_woutEVD2 = list()
m_woutRU2 = list()
m_woutTS2 = list()
rr$ru = (1-rr$probA) - (1-rr$probB)
rr$tu = sqrt((1-rr$probA)+(1-rr$probB))
rr$ts = rr$evDiff / rr$tu
for(s in 1:length(subs)){
  res = rr[rr$subjid==subs[s],]
  mfull2[[s]] = glm(evHigherChosen~evDiff+ru+ts,family=binomial(link = "probit"),res)
  m_woutEVD2[[s]] = glm(evHigherChosen~ru+ts,family=binomial(link = "probit"),res)
  m_woutRU2[[s]] = glm(evHigherChosen~evDiff+ts,family=binomial(link = "probit"),res)
  m_woutTS2[[s]] = glm(evHigherChosen~evDiff+ru,family=binomial(link = "probit"),res)
}
bicmat2 =aicmat = matrix(0,length(subs),4)
for(i in 1:length(subs)){
  bicmat2[i,1] = BIC(mfull2[[i]])
  bicmat2[i,2] = BIC(m_woutEVD2[[i]])
  bicmat2[i,3] = BIC(m_woutRU2[[i]])
  bicmat2[i,4] = BIC(m_woutTS2[[i]])
  aicmat[i,1] = AIC(mfull2[[i]])
  aicmat[i,2] = AIC(m_woutEVD2[[i]])
  aicmat[i,3] = AIC(m_woutRU2[[i]])
  aicmat[i,4] = AIC(m_woutTS2[[i]])
}

bic_tot = cbind(bicmat,bicmat2)
bics_1 = bic_tot[conds==1,]
bics_2 = bic_tot[conds==2,]
m = akaike.weights(bics_1)$weights
m = matrix(m,nrow(bics_1),8,byrow=F)
pxp = VB_bms(log(m))
barplot(pxp$pxp)
pxp1 = data.frame(PXP = pxp$pxp[5:8],
                 model = c("full","no EVD","no RU","no TS"))

m = akaike.weights(bics_2)$weights
m = matrix(m,nrow(bics_2),8,byrow=F)
pxp = VB_bms(log(m))
barplot(pxp$pxp)

pxp2 = data.frame(PXP = pxp$pxp[5:8],
                 model = c("full","no EVD","no RU","no TS"))
pxp = rbind(pxp1,pxp2)
pxp$condition = rep(1:2,each=4)
ggplot(pxp,aes(x=model,y=PXP,fill=model))+
  geom_bar(stat="identity")+
  scale_fill_viridis(discrete = T)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~condition)
aic = data.frame(AIC = c(aicmat[,1],aicmat[,2],aicmat[,3],aicmat[,4]),
                 subs = rep(1:62,4),
                 model = rep(c("full","no EVD","no RU","no TS"),each=62),
                 condition = rep(conds,4))

library(RColorBrewer)
library(ggbeeswarm)
library(viridis)
ggplot(aic,aes(x=model,y=AIC,fill=model))+
  geom_boxplot()+
  geom_beeswarm(aes(fill=model,colour=model),color="black",shape=21,dodge.width = 1.1,
                priority = "density")+
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~condition)

bic = data.frame(BIC = c(bicmat2[,1],bicmat2[,2],bicmat2[,3],bicmat2[,4]),
                 subs = rep(1:62,4),
                 model = rep(c("full","no EVD","no RU","no TS"),each=62),
                 condition = rep(conds,4))
ggplot(bic,aes(x=model,y=BIC,fill=model))+
  geom_boxplot()+
  geom_beeswarm(aes(fill=model,colour=model),color="black",shape=21,dodge.width = 1.1,
                priority = "density")+
  scale_fill_viridis(discrete = T)+
  scale_y_continuous(limits = c(),expand= c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~condition)
