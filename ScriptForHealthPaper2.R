library(rlist)
library(igraph)

library(boot)

#where functions and parameter csvs are located
sourcefolder <- "R"
paramsfolder <- "params"

#where the networks are located
networksfolder <-file.path(paramsfolder,"networks2")

#where you want to save the results (already created)
outputfolder<-file.path("output","results4")

source(file.path(sourcefolder,"FunctionsForHealthPaper.R"))

params1<-read.csv(file.path(paramsfolder,"model_params.csv"))
params2<-read.csv(file.path(paramsfolder,"model_params2.csv"))
params3<-read.csv(file.path(paramsfolder,"he_params.csv"))

net_params<-read.csv(file.path(paramsfolder,"network_params.csv"))

#h_che<-1

#start loop over networks
for(nt in 1:9){

par_id<-params1[nt,2]

pop_size=2000
ncomms=10
prop_belA=net_params$prop_belA[par_id]
prop_old=0.13
prop_young=0.63
prop_child=0.24

pop_info<-pop_gen(pop_size,ncomms,prop_belA,prop_old,prop_young,prop_child)

############################################

dis_input<-readRDS(file.path(networksfolder,paste0(params1[nt,3],"net_and_parents.RDS")))
info_input<-readRDS(file.path(networksfolder,paste0(params1[nt,2],"net_and_parents.RDS")))

parents<-info_input[[2]]
dis_mat<-dis_input[[1]]
info_mat<-info_input[[1]]

dis_mat<-par_ex(pop_info=pop_info,parents=parents,dis_mat=dis_mat)

############################################

for(md in 3:77){

# Here we define the prior beliefs of young adults (which will be used as probabilities in a bernoulli draw)
# e.g. currently there is a 50% chance a young adult of political belief A is concerned about the virus
A_concern_y<-params2[md,2]
B_concern_y<-params2[md,3]

#Here we define an additive effect of being old (to accomodate the fact they may be more likely to start concerned)
# e.g. There will be a 70% chance of an old adult of political belief A being concerned about the virus
A_concern_o<-params2[md,4]
B_concern_o<-params2[md,5]

#Here we define a daily extrinsic input into the belief of each political belief (I figured this would suffice to represent exposure to politicians/news/wider social media)
#N.B. These numbers are already defined on a logit scale. But our starting assumption is that concern of political belief A gets puched up and political belief B gets ushed down by these extrinsic factors
#(can obviously set to zero if preferred)
lA_ex<- params2[md,6]
lB_ex<- params2[md,7]

#Define a linear relationship between proportion of connections concerned (at the previous time-step) and concern levels in young adults
l_conc<-params2[md,8]
#and an additive effect used to calculate the same parameter for old adults
l_conc_o<-params2[md,9]

#Define a linear relationship between number of connections infected (at the previous time-step) and concern levels in young adults
l_inf<-params2[md,10]
#and an additive effect used to calculate the same parameter for old adults
l_inf_o<-params2[md,11]

for(r in 1:5){

#Define a threshold relationship whereby concern decreases while all immediate network conncections are fully healthy (at the previous time-step) and concern levels in young adults
l_hea<-params3[(nt-1)*77*5+(md-1)*5+r,3]
#and an additive effect used to calculate the same parameter for old adults
l_hea_o<-params2[md,13]

#h_che<-h_che+1

start<-concern_setup(A_concern_y,B_concern_y,A_concern_o,B_concern_o,pop_info,info_mat)

concern<-list()
belief<-list()

belief[[1]]<-start[[1]]
concern[[1]]<-start[[2]]

############################################

#Probability of becoming exposed having contacted an infectious individual (daily)
#Need to work out R0 based on other parameters
S_E<-0.3/mean(colSums(dis_mat))
#S_E<-0.2

#lambda for a Poisson draw for the length of this period
E_I1<-5.1

#probability of transitioning to serious case for young (daily)
yI1_I2<-0.01
#and same for old
oI1_I2<-0.05

#probability of transitioning to critical (HOSPITALISED) case for young (daily)
yI2_I3<-0.0125
#and same for old
oI2_I3<-0.025

#probability of death for young (daily)
yI3_D<-0.012
#and same for old
oI3_D<-0.092

#lambda for a Poisson draw for duration of a pre-symptomatic/mild infection - this is now misnamed as impossible to recover - simply transition to I2.
yI1_R<-6.7
#and same for old
oI1_R<-6.7

#lambda for a Poisson draw for duration of a serious infection
yI2_R<-10
#and same for old
oI2_R<-10

#lambda for a Poisson draw for duration of a critical/hospitalised infection
yI3_R<-4.2
#and same for old
oI3_R<-4.2

############################################

#start with 2 infected individuals
exp<-sample(1:pop_info$pop,5,replace=FALSE)

#create dataframe to store disease state
S<-rep(1,pop_info$pop)
E<-rep(0,pop_info$pop)
I1<-rep(0,pop_info$pop)
I2<-rep(0,pop_info$pop)
I3<-rep(0,pop_info$pop)
R<-rep(0,pop_info$pop)
D<-rep(0,pop_info$pop)
status<-data.frame(S,E,I1,I2,I3,R,D)

status$S[exp]<-0
status$E[exp]<-1

d_exp<-rep(NA,pop_info$pop)
d_inf1<-rep(NA,pop_info$pop)
d_inf2<-rep(NA,pop_info$pop)
d_inf3<-rep(NA,pop_info$pop)

d_exp[status$E==1]<-rpois(sum(status$E),E_I1)

time<-300

statuses<-list()
statuses[[1]]<-status

progression<-matrix(0,nr=time+1,nc=ncol(status))
progression[1,]<-colSums(status)


############################################

#be ready to change the network in this
for(t in 2:time){
  
  if(t==2){
    
    #time-step 2
    
    dis_mat<-network_rewire_concern(net=dis_mat,concern=start[[2]],concern.prev=rep(0,pop_info$pop),pop_info=pop_info,
                                    prop_cut=0.5,cut_to=0.001)
    dis_mat<-network_rewire_infectionS(net=dis_mat,status=statuses[[t-1]],
                                       pop_info=pop_info,cut_to=0.001)
    #dis_mat<-network_rewire_infectionM(net=dis_mat,status=statuses[[t-1]],
    #                                   pop_info=pop_info,cut_to=0.001)
    dis_mat<-network_rewire_infectionR(net=dis_mat,status=statuses[[t-1]],
                                       pop_info=pop_info,cut_to=0.001)
    
    inf<-cbind(sign(statuses[[t-1]]$I2+statuses[[t-1]]$I3),sign(statuses[[t-1]]$I2+statuses[[t-1]]$I3))
    
    current<-concern_timestep(pop_info=pop_info,net_b=info_mat,net_d=dis_mat,belief=start[[1]],concern=start[[2]],inf=inf,lA_ex,lB_ex,l_conc,l_conc_o,l_inf,l_inf_o,l_hea,l_hea_o)
    
    belief[[2]]<-current[[1]]
    concern[[2]]<-current[[2]]
    
    dis<-infection_timestep(pop_info=pop_info,status=statuses[[t-1]],net=dis_mat,d_exp=d_exp,d_inf1=d_inf1,d_inf2=d_inf2,d_inf3=d_inf3,S_E=S_E,E_I1=E_I1,yI1_I2=yI1_I2,oI1_I2=oI1_I2,yI2_I3=yI2_I3,oI2_I3=oI2_I3,yI3_D=yI3_D,oI3_D=oI3_D,yI1_R=yI1_R,oI1_R=oI1_R,yI2_R=yI2_R,oI2_R=oI2_R,yI3_R=yI3_R,oI3_R=oI3_R)
    
  }
  if(t>2){
    
    dis_mat<-network_rewire_concern(net=dis_mat,concern=current[[2]],concern.prev=concern[[t-2]],pop_info=pop_info,
                                    prop_cut=0.5,cut_to=0.001)
    dis_mat<-network_rewire_infectionS(net=dis_mat,status=statuses[[t-1]],
                                       pop_info=pop_info,cut_to=0.001)
    #dis_mat<-network_rewire_infectionM(net=dis_mat,status=statuses[[t-1]],
    #                                   pop_info=pop_info,cut_to=0.001)
    dis_mat<-network_rewire_infectionR(net=dis_mat,status=statuses[[t-1]],
                                       pop_info=pop_info,cut_to=0.001)
    
    inf<-cbind(sign(statuses[[t-1]]$I2+statuses[[t-1]]$I3),sign(statuses[[t-1]]$I2+statuses[[t-1]]$I3))
    current<-concern_timestep(pop_info=pop_info,net_b=info_mat,net_d=dis_mat,belief=current[[1]],
                              concern=current[[2]],inf=inf,lA_ex,lB_ex,l_conc,l_conc_o,l_inf,
                              l_inf_o,l_hea,l_hea_o)
    
    belief[[t]]<-current[[1]]
    concern[[t]]<-current[[2]]  
    
    dis<-infection_timestep(pop_info=pop_info,status=statuses[[t-1]],net=dis_mat,d_exp=dis[[2]],d_inf1=dis[[3]],d_inf2=dis[[4]],d_inf3=dis[[5]],S_E=S_E,E_I1=E_I1,yI1_I2=yI1_I2,oI1_I2=oI1_I2,yI2_I3=yI2_I3,oI2_I3=oI2_I3,yI3_D=yI3_D,oI3_D=oI3_D,yI1_R=yI1_R,oI1_R=oI1_R,yI2_R=yI2_R,oI2_R=oI2_R,yI3_R=yI3_R,oI3_R=oI3_R)
  }
  print(colSums(dis$status))
  progression[t,]<-colSums(dis$status)
  statuses[[t]]<-dis$status
  
  if(sum(colSums(statuses[[t]])[c(2,3,4,5)])==0){break()}
  
}


########################################
########################################

mod_concerns<-matrix(0,nr=10,nc=length(concern))

for(i in 1:length(concern)){
  mod_concerns[,i]<-aggregate(concern[[i]],by=list(pop_info$comms),mean)[,2]
}

mod_exps<-matrix(0,nr=10,nc=length(statuses))

for(i in 1:length(statuses)){
  mod_exps[,i]<-aggregate(statuses[[i]][,2],by=list(pop_info$comms),sum)[,2]
}

mod_infs<-matrix(0,nr=10,nc=length(statuses))

for(i in 1:length(statuses)){
  mod_infs[,i]<-aggregate(statuses[[i]][,4],by=list(pop_info$comms),sum)[,2]
}

mod_hosps<-matrix(0,nr=10,nc=length(statuses))

for(i in 1:length(statuses)){
  mod_hosps[,i]<-aggregate(statuses[[i]][,5],by=list(pop_info$comms),sum)[,2]
}

OUT<-list(mod_concerns,mod_exps,mod_infs,mod_hosps,l_hea)
names(OUT)<-c("concern","exps","infs","hosps","he")

saveRDS(OUT, (outputfolder,paste0("nets",params1$NetSelect[nt],"mods",params2$ModSelect[md],"rep",r,".RDS")))

###################################
###################################

cols=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")

plot(NULL,xlim=c(0,250),ylim=c(0,1))
for(i in 1:10){
  lines(x=seq(1,length(concern)),y=mod_concerns[i,],col=cols[i],lwd=3)
}

plot(NULL,xlim=c(0,250),ylim=c(0,50))
for(i in 1:10){
  lines(x=seq(1,length(statuses)),y=mod_infs[i,],col=cols[i],lwd=3)
}

} #end r loop

} #end md loop

} #end nt loop 






