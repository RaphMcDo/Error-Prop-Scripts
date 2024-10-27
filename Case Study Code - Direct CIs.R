#########################################################
############    MW-SH Modelling #########################
#########################################################

library(ggplot2)
library(dplyr)
library(INLA)
library(TMB)
library(sf)
library(raster)
library(mgcv)

load("all_mw.RData")
gb_2023_mw<-all_mw %>% st_as_sf(coords=c("lon","lat"))
st_crs(gb_2023_mw)<-4326
gb_2023_mw<-st_transform(gb_2023_mw,crs=32619)
gb_2023_mw$geometry<-gb_2023_mw$geometry/1000

gb_2023_surv<-st_as_sf(bigbank,coords=c("lon","lat"))
st_crs(gb_2023_surv)<-4326
gb_2023_surv<-st_transform(gb_2023_surv,crs=32619)
gb_2023_surv$geometry<-gb_2023_surv$geometry/1000
gb_2023_surv<-gb_2023_surv[which(gb_2023_surv$state=="live"),]

unique_surv_locs<-st_coordinates(gb_2023_surv)
unique_mw_locs<-unique(st_coordinates(gb_2023_mw))
mw_matches<-rep(NA,nrow(unique_mw_locs))
for (i in 1:nrow(unique_mw_locs)){
  mw_matches[i]<-which(unique_surv_locs[,1]==unique_mw_locs[i,1] & unique_surv_locs[,2]==unique_mw_locs[i,2])
}

coords_mw<-st_coordinates(gb_2023_mw)
gb_2023_mw$id_loc2<-rep(NA,nrow(gb_2023_mw))
for (i in 1:nrow(gb_2023_mw)){
  gb_2023_mw$id_loc2[i]<-mw_matches[which(unique_mw_locs[,1]==coords_mw[i,1] & unique_mw_locs[,2]==coords_mw[i,2])]
}

sub_mw<-gb_2023_mw

#####################################
######## LW MODELS ################
####################################


#Off method

# dyn.unload(dynlib("glmm_offshore_fixed_b"))
compile("glmm_offshore_fixed_b.cpp")
dyn.load(dynlib("glmm_offshore_fixed_b"))

tmb_data<-list(weight=sub_mw$wmw,
               heights=sub_mw$sh/100,
               tow_id=as.integer(as.factor(sub_mw$ID))-1)

par_list<-list(b=3,
               beta=20,
               log_phi=-1,
               log_epsilon=-1,
               tow_eff=rep(0,length(unique(sub_mw$ID))))

random<-c("tow_eff")

map<-list(b=factor(NA))

obj = MakeADFun(data=tmb_data, 
                parameters=par_list,
                map=map,
                random=random,
                DLL="glmm_offshore_fixed_b",
                silent = F)

Opt<-optimx::optimr(obj$par,obj$fn,obj$gr,
                    control=list(maxit=100000,maxeval=100000),
                    method="nlminb")

rep<-sdreport(obj,bias.correct=F)

Report<-obj$report()

se_comb<-rep(NA,length(Report$tow_eff))
for (i in 1:length(se_comb)){
  se_comb[i]<-sqrt(sum(rep$cov[c(which(names(rep$value)=="beta"),which(names(rep$value)=="tow_eff")[i]),
                              c(which(names(rep$value)=="beta"),which(names(rep$value)=="tow_eff")[i])]))
}

temp_fit<-data.frame(ID=unique(sub_mw$ID),CF=(Report$beta+Report$tow_eff),se=se_comb)

sub_mw$lon<-st_coordinates(sub_mw)[,1]
sub_mw$lat<-st_coordinates(sub_mw)[,2]
CF.data<-merge(sub_mw[!duplicated(sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit)

names(CF.data)<-c('ID','lon','lat','year','depth','tow','CF',"se")

CF.fit<-gam(CF~s(lon,lat)+s(depth),data=CF.data)

#for pred, get CF at each loc, so get mean lon and lat for survey
gb_2023_surv$lon<-st_coordinates(gb_2023_surv)[,1]
gb_2023_surv$lat<-st_coordinates(gb_2023_surv)[,2]

gb_2023_surv$ID<-paste(as.character(gb_2023_surv$cruise),as.character(gb_2023_surv$tow),sep=".")

pred.loc<-data.frame(depth=gb_2023_surv$depth,lon=gb_2023_surv$lon,lat=gb_2023_surv$lat)

Cf.pred<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
Cf.pred$CF<-predict(CF.fit,Cf.pred, se=T)$fit
Cf.pred$se_gam<-predict(CF.fit,Cf.pred, se=T)$se.fit
#Checked and the triangular elements are so small so 

surv.dat<-Cf.pred
surv.dat$ID<-gb_2023_surv$ID

# Pull out the ID and condition factor
tmp.dat<-subset(st_drop_geometry(CF.data),select=c("ID","CF","se"))
# Rename CF to CFh
names(tmp.dat)[2]<-"CFh"
# merge the two data sets, keeping all x values
surv.dat<-merge(surv.dat,tmp.dat,all.x=T)
# Replace any NA's in CFh with the original Condition Factor.
surv.dat$CFh[is.na(surv.dat$CFh)]<-surv.dat$CF[is.na(surv.dat$CFh)]
surv.dat$se[is.na(surv.dat$se)]<-surv.dat$se_gam[is.na(surv.dat$se)]

#Off LN method

# dyn.unload(dynlib("glmm_offshore_fixed_b_lognormal"))
compile("glmm_offshore_fixed_b_lognormal.cpp")
dyn.load(dynlib("glmm_offshore_fixed_b_lognormal"))

tmb_data<-list(weight=sub_mw$wmw,
               heights=sub_mw$sh/100,
               tow_id=as.integer(as.factor(sub_mw$ID))-1)

par_list<-list(b=3,
               beta=20,
               log_phi=-1,
               log_epsilon=-1,
               tow_eff=rep(0,length(unique(sub_mw$ID))))

random<-c("tow_eff")

map<-list(b=factor(NA))

obj2 = MakeADFun(data=tmb_data, 
                 parameters=par_list,
                 map=map,
                 random=random,
                 DLL="glmm_offshore_fixed_b_lognormal",
                 silent = F)

Opt2<-optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep2<-sdreport(obj2,bias.correct=F)

Report2<-obj2$report()

se_comb2<-rep(NA,length(Report2$tow_eff))
for (i in 1:length(se_comb)){
  se_comb2[i]<-sqrt(sum(rep2$cov[c(which(names(rep2$value)=="beta"),which(names(rep2$value)=="tow_eff")[i]),
                               c(which(names(rep2$value)=="beta"),which(names(rep2$value)=="tow_eff")[i])]))
}

temp_fit2<-data.frame(ID=unique(sub_mw$ID),CF=(Report2$beta+Report2$tow_eff),se=se_comb2)

CF.data2<-merge(sub_mw[!duplicated(sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit2)

names(CF.data2)<-c('ID','lon','lat','year','depth','tow','CF',"se")

CF.fit2<-gam(CF~s(lon,lat)+s(depth),data=CF.data2)

Cf.pred2<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
Cf.pred2$CF<-predict(CF.fit2,Cf.pred2, se=T)$fit
Cf.pred2$se_gam<-predict(CF.fit2,Cf.pred2, se=T)$se.fit

surv.dat2<-Cf.pred2
surv.dat2$ID<-gb_2023_surv$ID

# Pull out the ID and condition factor
tmp.dat2<-subset(st_drop_geometry(CF.data2),select=c("ID","CF","se"))
# Rename CF to CFh
names(tmp.dat2)[2]<-"CFh"
# merge the two data sets, keeping all x values
surv.dat2<-merge(surv.dat2,tmp.dat2,all.x=T)
# Replace any NA's in CFh with the original Condition Factor.
surv.dat2$CFh[is.na(surv.dat2$CFh)]<-surv.dat2$CF[is.na(surv.dat2$CFh)]
surv.dat2$se[is.na(surv.dat2$se)]<-surv.dat2$se_gam[is.na(surv.dat2$se)]

#Depth method

# dyn.unload(dynlib("depth_glmm_no_space"))
compile("depth_glmm_no_space.cpp")
dyn.load(dynlib("depth_glmm_no_space"))

tmb_data<-list(weight=sub_mw$wmw,
               varmat=as.matrix(cbind(1, log(sub_mw$sh/100))),
               depth=sub_mw$depth,
               ind_loc=sub_mw$id_loc-1)

nvar = ncol(tmb_data$varmat)
nobs = nrow(tmb_data$varmat)

parameters = list(
  beta = rep(0, nvar),
  beta_depth = rep(0, nvar),
  log_phi = log(0.1))

maps <- list(
  beta = factor(c(NA,2)),
  beta_depth=factor(c(1,NA))
)

obj3 = MakeADFun(data=tmb_data, 
                 parameters=parameters,
                 map=maps,
                 random=c(),
                 DLL="depth_glmm_no_space",
                 silent = F)

Opt3<-optimx::optimr(obj3$par,obj3$fn,obj3$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep3<-sdreport(obj3,bias.correct=F)

Report3<-obj3$report()

#Spatial

# dyn.unload(dynlib("spatial_depth_glmm_direct_cov"))
compile("spatial_depth_glmm_direct_cov.cpp")
dyn.load(dynlib("spatial_depth_glmm_direct_cov"))

tmb_data<-list(weight=sub_mw$wmw,
               varmat=as.matrix(cbind(1, log(sub_mw$sh/100))),
               depth=sub_mw$depth,
               ind_loc=sub_mw$id_loc2-1,
               locations=st_coordinates(gb_2023_surv))

parameters = list(
  beta = rep(0, nvar),
  beta_depth = rep(0, nvar),
  beta_s = rep(0,nrow(tmb_data$locations)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_phi = log(0.1))

maps <- list(
  beta = factor(c(NA,2)),
  beta_depth=factor(c(1,NA)),
  log_nu=factor(NA)
)

obj4 = MakeADFun(data=tmb_data, 
                 parameters=parameters,
                 map=maps,
                 random=c("beta_s"),
                 DLL="spatial_depth_glmm_direct_cov",
                 silent = F)

Opt4<-optimx::optimr(obj4$par,obj4$fn,obj4$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep4<-sdreport(obj4,bias.correct=F)

Report4<-obj4$report()

#Setup for Shell height model
sub_mw$size_bin<-as.integer(cut(sub_mw$sh,breaks=seq(60,160,by=5),
                                include.lowest=T))-1
sub_mw$size_bin_name<-cut(sub_mw$sh,breaks=seq(60,170,by=5),
                          include.lowest=T)

#Shell height portion of model (only in supplementary materials)

# dyn.unload(dynlib("spatial_depth_glmm_height"))
compile("spatial_depth_glmm_height.cpp")
dyn.load(dynlib("spatial_depth_glmm_height"))

non_corrected_heights<-read.csv("DR2023_13.csv")
non_corrected_heights<-subset(non_corrected_heights,BIN_ID>95)
non_corrected_heights$DEPTH_M<-non_corrected_heights$DEPTH_F*1.8288
tot_scallops<-sum(non_corrected_heights$LIVE_RAW)
long_heights<-data.frame(heights=rep(NA,tot_scallops),depth=rep(NA,tot_scallops),
                         tow_id=rep(NA,tot_scallops))
counter<-1
for (i in 1:nrow(non_corrected_heights)){
  if (non_corrected_heights$LIVE_RAW[i]>0) {
    for (g in 1:non_corrected_heights$LIVE_RAW[i]){
      long_heights[counter,]<-c(runif(1,non_corrected_heights$BIN_ID[i]-5,non_corrected_heights$BIN_ID[i]),
                                non_corrected_heights$DEPTH_M[i],
                                which(st_drop_geometry(gb_2023_surv$tow)==non_corrected_heights$TOW_NO[i]))
      counter<-counter+1
    }
  }
}

unique_means<-rep(NA,length(unique(long_heights$tow_id)))

for (i in 1:length(unique_means)){
  unique_means[i]<-which(long_heights$tow_id==unique(long_heights$tow_id)[i])[1]
}

tmb_data<-list(varmat_sh=as.matrix(cbind(rep(1,nrow(long_heights)))),
               depth_sh=long_heights$depth,
               ind_loc_sh=long_heights$tow_id-1,
               locations=st_coordinates(gb_2023_surv),
               s_heights=long_heights$heights/10,
               a_b_truncation=c(9.5,17),
               unique_means=unique_means,
               n_mod_loc=length(unique(long_heights$tow_id)))

parameters = list(
  beta_depth = c(0),
  beta_sh = c(14),
  beta_sh_s = rep(0,nrow(tmb_data$locations)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_upsilon = log(0.5))

maps <- list(
  log_nu=c(factor(NA))#
)

obj5 = MakeADFun(data=tmb_data,
                 parameters=parameters,
                 map=maps,
                 random=c("beta_sh_s"),
                 DLL="spatial_depth_glmm_height",
                 silent = F)

Opt5<-optimx::optimr(obj5$par,obj5$fn,obj5$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep5<-sdreport(obj5,bias.correct=F)

Report5<-obj5$report()

#Spatial b

# dyn.unload(dynlib("spatial_depth_glmm_offshore"))
compile("spatial_depth_glmm_offshore.cpp")
dyn.load(dynlib("spatial_depth_glmm_offshore"))

tmb_data<-list(weight=sub_mw$wmw,
               varmat=as.matrix(cbind(1, log(sub_mw$sh/100))),
               depth=sub_mw$depth,
               ind_loc=sub_mw$id_loc2-1,
               locations=st_coordinates(gb_2023_surv))

parameters = list(
  beta = rep(0, nvar),
  beta_depth = rep(0, nvar),
  beta_s = rep(0,nrow(tmb_data$locations)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_phi = log(0.1))

maps <- list(
  # beta = factor(rep(NA, nvar)),
  # beta = factor(c(2,NA)),
  beta = factor(c(NA,2)),
  # beta = factor(c(3,2)),
  beta_depth=factor(c(1,NA)),
  log_nu=factor(NA)
)

obj6 = MakeADFun(data=tmb_data, 
                 parameters=parameters,
                 map=maps,
                 random=c("beta_s"),
                 DLL="spatial_depth_glmm_offshore",
                 silent = F)

Opt6<-optimx::optimr(obj6$par,obj6$fn,obj6$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep6<-sdreport(obj6,bias.correct=F)

Report6<-obj6$report()

#Off, but estimating b

tmb_data<-list(weight=sub_mw$wmw,
               heights=sub_mw$sh/100,
               tow_id=as.integer(as.factor(sub_mw$ID))-1)

par_list<-list(b=3,
               beta=20,
               log_phi=-1,
               log_epsilon=-1,
               tow_eff=rep(0,length(unique(sub_mw$ID))))

random<-c("tow_eff")

map<-list()

obj7 = MakeADFun(data=tmb_data, 
                 parameters=par_list,
                 map=map,
                 random=random,
                 DLL="glmm_offshore_fixed_b",
                 silent = F)

Opt7<-optimx::optimr(obj7$par,obj7$fn,obj7$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep7<-sdreport(obj7,bias.correct=F)

Report7<-obj7$report()

se_comb7<-rep(NA,length(Report7$tow_eff))
for (i in 1:length(se_comb)){
  se_comb7[i]<-sqrt(sum(rep7$cov[c(which(names(rep7$value)=="beta"),which(names(rep7$value)=="tow_eff")[i],which(names(rep7$value)=="b")),
                                 c(which(names(rep7$value)=="beta"),which(names(rep7$value)=="tow_eff")[i],which(names(rep7$value)=="b"))]))
}

temp_fit7<-data.frame(ID=unique(sub_mw$ID),CF=(Report7$beta+Report7$tow_eff),se=se_comb7)

CF.data7<-merge(sub_mw[!duplicated(sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit7)

names(CF.data7)<-c('ID','lon','lat','year','depth','tow','CF',"se")

CF.fit7<-gam(CF~s(lon,lat)+s(depth),data=CF.data7)

Cf.pred7<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
Cf.pred7$CF<-predict(CF.fit7,Cf.pred7, se=T)$fit
Cf.pred7$se_gam<-predict(CF.fit7,Cf.pred7, se=T)$se.fit

surv.dat7<-Cf.pred7
surv.dat7$ID<-gb_2023_surv$ID

# Pull out the ID and condition factor
tmp.dat7<-subset(st_drop_geometry(CF.data7),select=c("ID","CF","se"))
# Rename CF to CFh
names(tmp.dat7)[2]<-"CFh"
# merge the two data sets, keeping all x values
surv.dat7<-merge(surv.dat7,tmp.dat7,all.x=T)
# Replace any NA's in CFh with the original Condition Factor.
surv.dat7$CFh[is.na(surv.dat7$CFh)]<-surv.dat7$CF[is.na(surv.dat7$CFh)]
surv.dat7$se[is.na(surv.dat7$se)]<-surv.dat7$se_gam[is.na(surv.dat7$se)]

#Off LN, estimating b

# dyn.unload(dynlib("glmm_offshore_fixed_b_lognormal"))
compile("glmm_offshore_fixed_b_lognormal.cpp")
dyn.load(dynlib("glmm_offshore_fixed_b_lognormal"))

tmb_data<-list(weight=sub_mw$wmw,
               heights=sub_mw$sh/100,
               tow_id=as.integer(as.factor(sub_mw$ID))-1)

par_list<-list(b=3,
               beta=20,
               log_phi=-1,
               log_epsilon=-1,
               tow_eff=rep(0,length(unique(sub_mw$ID))))

random<-c("tow_eff")

obj8 = MakeADFun(data=tmb_data, 
                 parameters=par_list,
                 map=map,
                 random=random,
                 DLL="glmm_offshore_fixed_b_lognormal",
                 silent = F)

Opt8<-optimx::optimr(obj8$par,obj8$fn,obj8$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep8<-sdreport(obj8,bias.correct=F)

Report8<-obj8$report()

se_comb8<-rep(NA,length(Report8$tow_eff))
for (i in 1:length(se_comb)){
  se_comb8[i]<-sqrt(sum(rep8$cov[c(which(names(rep8$value)=="beta"),which(names(rep8$value)=="tow_eff")[i],which(names(rep8$value)=="b")),
                                 c(which(names(rep8$value)=="beta"),which(names(rep8$value)=="tow_eff")[i],which(names(rep8$value)=="b"))]))
}

temp_fit8<-data.frame(ID=unique(sub_mw$ID),CF=(Report8$beta+Report8$tow_eff))

CF.data8<-merge(sub_mw[!duplicated(sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit8)

names(CF.data8)<-c('ID','lon','lat','year','depth','tow','CF',"se")

CF.fit8<-gam(CF~s(lon,lat)+s(depth),data=CF.data8)

Cf.pred8<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
Cf.pred8$CF<-predict(CF.fit8,Cf.pred8, se=T)$fit
Cf.pred8$se_gam<-predict(CF.fit8,Cf.pred8, se=T)$se.fit

surv.dat8<-Cf.pred8
surv.dat8$ID<-gb_2023_surv$ID

# Pull out the ID and condition factor
tmp.dat8<-subset(st_drop_geometry(CF.data8),select=c("ID","CF","se"))
# Rename CF to CFh
names(tmp.dat8)[2]<-"CFh"
# merge the two data sets, keeping all x values
surv.dat8<-merge(surv.dat,tmp.dat8,all.x=T)
# Replace any NA's in CFh with the original Condition Factor.
surv.dat8$CFh[is.na(surv.dat8$CFh)]<-surv.dat8$CF[is.na(surv.dat8$CFh)]
surv.dat8$se[is.na(surv.dat8$se)]<-surv.dat8$se_gam[is.na(surv.dat8$se)]

#Spatial Both

# dyn.unload(dynlib("spatial_both"))
compile("spatial_both.cpp")
dyn.load(dynlib("spatial_both"))

tmb_data<-list(weight=sub_mw$wmw,
               varmat=as.matrix(cbind(1, log(sub_mw$sh/100))),
               depth=sub_mw$depth,
               ind_loc=sub_mw$id_loc2-1,
               locations=st_coordinates(gb_2023_surv))

parameters = list(
  beta = rep(0, nvar),
  beta_depth = rep(0, nvar),
  beta_s = rep(0,nrow(tmb_data$locations)),
  beta_b_s = rep(0,nrow(tmb_data$locations)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_rho_b = log(10),
  log_sig_b = log(0.1),
  log_nu_b = log(1),
  log_phi = log(0.1))

maps <- list(
  beta = factor(c(NA,2)),
  beta_depth=factor(c(1,NA)),
  log_nu=factor(NA),
  log_nu_b=factor(NA)
)

non_r<-names(parameters)[-which(names(parameters) %in% c("beta_s","beta_b_s"))]

obj9 = MakeADFun(data=tmb_data, 
                 parameters=parameters,
                 map=maps,
                 random=c("beta_s","beta_b_s"),
                 DLL="spatial_both",
                 silent = F)

Opt9<-optimx::optimr(obj9$par,obj9$fn,obj9$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep9<-sdreport(obj9,bias.correct=F)

Report9<-obj9$report()

#Now gotta predict the 3 possible ways:
#Predict midpoint of each bin, scale up
#Predict all scallops
#Predict value at mean predicted shell height, scale that up

midpoints<-seq(0.975,1.675,by=0.05)
mw_models<-c("Off",
             "Off LN",
             "Depth",
             "Spatial",
             "Spatial b",
             "Offshore Estimating b",
             "Offshore LN Estimating",
             "Spatial Both")

n_loc<-nrow(gb_2023_surv)
match_gam_gb<-rep(NA,234)
for (i in 1:length(match_gam_gb)){
  match_gam_gb[i]<-which(gb_2023_surv$ID==surv.dat$ID[i])
}

#Will use lognorm package, cite Lo 2013
library(lognorm)

#Do classical stuff, Var(X+Y)=Var(X)+Var(Y)+2Cov(X,Y)
#Var(A+B+C) = Var(A)+Var(B)+Var(C)+2Cov(A,B)+2Cov(A,C)+2Cov(B,C)
#See supplementary materials for details

pred_midpoints<-data.frame(model=rep(mw_models,each=length(midpoints)*n_loc),
                          midpoints=rep(rep(midpoints,each=nrow(gb_2023_surv)),length(mw_models)),
                          se=rep(midpoints,each=length(mw_models)*nrow(gb_2023_surv)),
                          depth=rep(gb_2023_surv$depth,length(mw_models)*length(midpoints)),
                          loc=rep(1:234,length(mw_models)*length(midpoints)),
                          gam_loc=rep(match_gam_gb,length(mw_models)*length(midpoints)),
                          geometry=rep(gb_2023_surv$geometry,length(mw_models)*length(midpoints)),
                          pred_mid=rep(NA,length(rep(mw_models,each=length(midpoints)*n_loc))))  %>% st_as_sf()

summary_rep3_cov<-rep3$cov[c(which(names(rep3$value)=="beta_depth"),which(names(rep3$value)=="beta")),
                        c(which(names(rep3$value)=="beta_depth"),which(names(rep3$value)=="beta"))][c(1,4),c(1,4)]


for (i in 1:nrow(pred_midpoints)){
  if (pred_midpoints$model[i]==mw_models[1]) {
    pred_midpoints$pred_mid[i]<-(surv.dat$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
    pred_midpoints$se[i]<-sqrt((surv.dat$se[pred_midpoints$gam_loc][i])^2 * (pred_midpoints$midpoints[i]^3)^2)
  }
  if (pred_midpoints$model[i]==mw_models[2]) {
    pred_midpoints$pred_mid[i]<-(surv.dat2$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
    pred_midpoints$se[i]<-sqrt((surv.dat2$se[pred_midpoints$gam_loc][i])^2 * (pred_midpoints$midpoints[i]^3)^2)
  }
  if (pred_midpoints$model[i]==mw_models[3]) {
    pred_midpoints$pred_mid[i]<-exp((summary(rep3)[2,1])*log(pred_midpoints$depth[i])+(summary(rep3)[1,1])*log(pred_midpoints$midpoints[i]))
    pred_midpoints$se[i]<-sqrt(sum(c(summary_rep3_cov[1,1] * log(pred_midpoints$depth[i])^2,
                                     summary_rep3_cov[2,2] * log(pred_midpoints$midpoints[i])^2,
                                     2*(summary_rep3_cov[1,2]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i])))))
  }  
  if (pred_midpoints$model[i]==mw_models[4]) {
    temp_cov<-rep4$cov[c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[pred_midpoints$loc][i]),
                       c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[pred_midpoints$loc][i])]
    pred_midpoints$pred_mid[i]<-exp((summary(rep4)[4,1])*log(pred_midpoints$depth[i])+(Report4$beta_s[pred_midpoints$loc])[i]+(summary(rep4)[3,1])*log(pred_midpoints$midpoints[i]))
    pred_midpoints$se[i]<-sqrt(sum(c(temp_cov[1,1]*log(pred_midpoints$depth[i])^2,
                                     temp_cov[3,3],
                                     temp_cov[2,2]*log(pred_midpoints$midpoints[i])^2),
                                     2*temp_cov[1,3]*log(pred_midpoints$depth[i]),
                                     2*temp_cov[1,2]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[2,3]*log(pred_midpoints$midpoints[i])))
  }  
  if (pred_midpoints$model[i]==mw_models[5]) {
    temp_cov<-rep6$cov[c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[pred_midpoints$loc][i]),
                       c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[pred_midpoints$loc][i])]
    pred_midpoints$pred_mid[i]<-exp((summary(rep6)[4,1])*log(pred_midpoints$depth[i])+((Report6$beta_s[pred_midpoints$loc])[i]+(summary(rep6)[3,1]))*log(pred_midpoints$midpoints[i]))
    pred_midpoints$se[i]<-sqrt(sum(c(temp_cov[1,1]*log(pred_midpoints$depth[i])^2,
                                     temp_cov[3,3]*log(pred_midpoints$midpoints[i])^2,
                                     temp_cov[2,2]*log(pred_midpoints$midpoints[i])^2),
                                     2*temp_cov[1,3]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[1,2]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[2,3]*log(pred_midpoints$midpoints[i])^2))
  }  
  if (pred_midpoints$model[i]==mw_models[6]) {
    pred_midpoints$pred_mid[i]<-(surv.dat7$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^(summary(rep7)[1,1])
    pred_midpoints$se[i]<-sqrt((surv.dat7$se[pred_midpoints$gam_loc][i])^2 * (pred_midpoints$midpoints[i]^(summary(rep7)[1,1]))^2)
  }  
  if (pred_midpoints$model[i]==mw_models[7]) {
    pred_midpoints$pred_mid[i]<-(surv.dat8$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^(summary(rep8)[1,1])
    pred_midpoints$se[i]<-sqrt((surv.dat8$se[pred_midpoints$gam_loc][i])^2 * (pred_midpoints$midpoints[i]^(summary(rep8)[1,1]))^2)
  } 
  if (pred_midpoints$model[i]==mw_models[8]) {
    temp_cov<-rep9$cov[c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[pred_midpoints$loc][i],which(names(rep9$value)=="beta_b_s")[pred_midpoints$loc][i]),
                       c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[pred_midpoints$loc][i],which(names(rep9$value)=="beta_b_s")[pred_midpoints$loc][i])]
    pred_midpoints$pred_mid[i]<-exp((summary(rep9)[6,1])*log(pred_midpoints$depth[i])+Report9$beta_s[pred_midpoints$loc[i]]+((Report9$beta_b_s[pred_midpoints$loc])[i]+(summary(rep9)[5,1]))*log(pred_midpoints$midpoints[i]))
    pred_midpoints$se[i]<-sqrt(sum(c(temp_cov[1,1]*log(pred_midpoints$depth[i])^2,
                                     temp_cov[4,4]*log(pred_midpoints$midpoints[i])^2,
                                     temp_cov[2,2]*log(pred_midpoints$midpoints[i])^2,
                                     temp_cov[3,3],
                                     2*temp_cov[1,4]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[1,2]*log(pred_midpoints$depth[i])*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[2,4]*log(pred_midpoints$midpoints[i])^2,
                                     2*temp_cov[1,3]*log(pred_midpoints$depth[i]),
                                     2*temp_cov[2,3]*log(pred_midpoints$midpoints[i]),
                                     2*temp_cov[3,4]*log(pred_midpoints$midpoints[i]))))
  }
} 

long_heights$cut_heights<-cut(long_heights$heights/100,breaks=c(midpoints-0.025,1.7))
long_heights$pred_heights<-Report5$mu_sh/10#+exp(summary(rep5)[5,1])*(dnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1]))-dnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1])))/(pnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1]))-pnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1])))
agg_long_heights<-aggregate(heights~cut_heights+tow_id,data=long_heights,FUN=length)
agg_long_heights$mid<-midpoints[as.integer(agg_long_heights$cut_heights)]
join_frame<-data.frame(tow_id=rep(1:234,each=length(midpoints)),mid=rep(midpoints,234))
agg_long_heights<-left_join(join_frame,agg_long_heights,by=c("tow_id","mid"))
agg_long_heights$heights[which(is.na(agg_long_heights$heights))]<-0
colnames(agg_long_heights)[4]<-"total"

agg_long_heights_se<-data.frame(agg_long_heights,off=rep(NA,nrow(agg_long_heights)),
                                off_ln=rep(NA,nrow(agg_long_heights)),insh=rep(NA,nrow(agg_long_heights)),
                                spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                off_est=rep(NA,nrow(agg_long_heights)),off_ln_est=rep(NA,nrow(agg_long_heights)),
                                spat_both=rep(NA,nrow(agg_long_heights)))
agg_long_heights<-data.frame(agg_long_heights,off=rep(NA,nrow(agg_long_heights)),
                            off_ln=rep(NA,nrow(agg_long_heights)),insh=rep(NA,nrow(agg_long_heights)),
                            spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                            off_est=rep(NA,nrow(agg_long_heights)),off_ln_est=rep(NA,nrow(agg_long_heights)),
                            spat_both=rep(NA,nrow(agg_long_heights)))

for (model in 1:length(mw_models)){
  for (i in 1:nrow(gb_2023_surv)){
    if (model %in% c(1,2,6,7)){
      agg_long_heights[which(agg_long_heights$tow_id==i),model+4]<-agg_long_heights[which(agg_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$gam_loc==i),]$pred_mid
      agg_long_heights_se[which(agg_long_heights_se$tow_id==i),model+4]<-agg_long_heights_se[which(agg_long_heights_se$tow_id==i),]$total * pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$gam_loc==i),]$se
    } else{
      for (j in 1:length(agg_long_heights[which(agg_long_heights$tow_id==i),model+4])){
        estims<-estimateSumLognormal(mu=rep(log(pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$pred_mid[j]),agg_long_heights[which(agg_long_heights$tow_id==i),]$total[j]),
                                     sigma=rep(pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$se[j],agg_long_heights[which(agg_long_heights$tow_id==i),]$total[j])^2)
        agg_long_heights[which(agg_long_heights$tow_id==i),model+4][j]<-exp(estims[1])
        agg_long_heights_se[which(agg_long_heights_se$tow_id==i),model+4][j]<-sqrt(estims[2])
      }
    }
  }
}

for (i in 5:12){
  agg_long_heights[,i][which(is.na(agg_long_heights[,i]))]<-0
}

agg_midpoints_scallops<-data.frame(tow_id=1:234,off=NA,off_ln=NA,insh=NA,spat=NA,spat_off=NA,
                                   off_est=NA,off_ln_est=NA,spat_both=NA)
agg_midpoints_scallops_se<-data.frame(tow_id=1:234,off=NA,off_ln=NA,insh=NA,spat=NA,spat_off=NA,
                                   off_est=NA,off_ln_est=NA,spat_both=NA)

for (tow in 1:234){
  for (model in 1:8){
    if (model %in% c(3,4,5,8)){
      estims<-estimateSumLognormal(mu=log(agg_long_heights[which(agg_long_heights$tow_id==tow),model+4][which(agg_long_heights[which(agg_long_heights$tow_id==tow),model+4]>0)]),
                                   sigma=agg_long_heights_se[which(agg_long_heights_se$tow_id==tow),model+4][which(!is.na(agg_long_heights_se[which(agg_long_heights_se$tow_id==tow),model+4]))]^2)
      agg_midpoints_scallops[tow,model+1]<-exp(estims[1])
      agg_midpoints_scallops_se[tow,model+1]<-sqrt(estims[2])
    }
  }
}
agg_midpoints_scallops$off<-aggregate(off~tow_id,data=agg_long_heights,FUN=sum)$off
agg_midpoints_scallops$off_ln<-aggregate(off_ln~tow_id,data=agg_long_heights,FUN=sum)$off_ln
agg_midpoints_scallops$off_est<-aggregate(off_est~tow_id,data=agg_long_heights,FUN=sum)$off_est
agg_midpoints_scallops$off_ln_est<-aggregate(off_ln_est~tow_id,data=agg_long_heights,FUN=sum)$off_ln_est

agg_midpoints_scallops_se$off<-aggregate(off~tow_id,data=agg_long_heights_se,FUN=sum)$off
agg_midpoints_scallops_se$off_ln<-aggregate(off_ln~tow_id,data=agg_long_heights_se,FUN=sum)$off_ln
agg_midpoints_scallops_se$off_est<-aggregate(off_est~tow_id,data=agg_long_heights_se,FUN=sum)$off_est
agg_midpoints_scallops_se$off_ln_est<-aggregate(off_ln_est~tow_id,data=agg_long_heights_se,FUN=sum)$off_ln_est

#All scallops

rematch_gam<-rep(NA,nrow(long_heights))
for (i in 1:nrow(long_heights)){
  rematch_gam[i]<-which(match_gam_gb==long_heights$tow_id[i])
}

new_long_heights<-long_heights
new_long_heights_se<-long_heights

for (model in 1:length(mw_models)){
  if (model == 1) {
    new_long_heights$off<-(surv.dat$CFh[rematch_gam])*(new_long_heights$heights/100)^3
    new_long_heights_se$off<-sqrt((surv.dat$se[rematch_gam]^2)*((new_long_heights$heights/100)^3)^2)
  }  
  if (model == 2) {
    new_long_heights$off_ln<-(surv.dat2$CFh[rematch_gam])*(new_long_heights$heights/100)^3
    new_long_heights_se$off_ln<-sqrt((surv.dat2$se[rematch_gam]^2)*((new_long_heights$heights/100)^3)^2)
  }  
  if (model == 3) {
    new_long_heights$insh<-exp((summary(rep3)[2,1])*log(new_long_heights$depth)+(summary(rep3)[1,1])*log(new_long_heights$heights/100))
    new_long_heights_se$insh<-sqrt(rowSums(data.frame(c(summary_rep3_cov[1,1] * log(new_long_heights$depth)^2),
                                                      c(summary_rep3_cov[2,2] *log(new_long_heights$heights/100)^2),
                                                      c(2*(summary_rep3_cov[1,2]* log(new_long_heights$depth)*log(new_long_heights$heights/100))))))
  }  
  if (model == 4) {
    for (i in 1:nrow(new_long_heights)){
      temp_cov<-rep4$cov[c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[new_long_heights$tow_id][i]),
                         c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[new_long_heights$tow_id][i])]
      new_long_heights$spat[i]<-exp((summary(rep4)[4,1])*log(new_long_heights$depth[i])+(Report4$beta_s[new_long_heights$tow_id][i])+(summary(rep4)[3,1])*log(new_long_heights$heights[i]/100))
      new_long_heights_se$spat[i]<-sqrt(sum(c(temp_cov[1,1]*log(new_long_heights$depth[i])^2,
                                           temp_cov[3,3],
                                           temp_cov[2,2]*log(new_long_heights$heights[i]/100)^2,
                                           2*temp_cov[1,3]*log(new_long_heights$depth[i]),
                                           2*temp_cov[1,2]*log(new_long_heights$depth[i])*log(new_long_heights$heights[i]/100),
                                           2*temp_cov[2,3]*log(new_long_heights$heights[i]/100))))
    }
  }  
  if (model == 5) {
    for (i in 1:nrow(new_long_heights)){
      temp_cov<-rep6$cov[c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[new_long_heights$tow_id][i]),
                         c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[new_long_heights$tow_id][i])]
      new_long_heights$spat_off[i]<-exp((summary(rep6)[4,1])*log(new_long_heights$depth[i])+((Report6$beta_s[new_long_heights$tow_id][i])+(summary(rep6)[3,1]))*log(new_long_heights$heights[i]/100))
      new_long_heights_se$spat_off[i]<-sqrt(sum(c(temp_cov[1,1]*log(new_long_heights$depth[i])^2),
                                                c(temp_cov[3,3]*log(new_long_heights$heights[i]/100)^2),
                                                c(temp_cov[2,2]*log(new_long_heights$heights[i]/100)^2),
                                                c(2*temp_cov[1,3]*log(new_long_heights$depth[i])*log(new_long_heights$heights[i]/100)),
                                                c(2*temp_cov[1,2]*log(new_long_heights$depth[i])*log(new_long_heights$heights[i]/100)),
                                                c(2*temp_cov[2,3]*log(new_long_heights$heights[i]/100)^2))) 
    }
  }  
  if (model == 6) {
    new_long_heights$off_est<-(surv.dat7$CFh[rematch_gam])*(new_long_heights$heights/100)^(summary(rep7)[1,1])
    new_long_heights_se$off_est<-sqrt(surv.dat7$se[rematch_gam]^2 * (new_long_heights$heights/100)^summary(rep7)[1,1])
  }  
  if (model == 7) {
    new_long_heights$off_ln_est<-(surv.dat8$CFh[rematch_gam])*(new_long_heights$heights/100)^(summary(rep8)[1,1])
    new_long_heights_se$off_ln_est<-sqrt(surv.dat8$se[rematch_gam]^2 * (new_long_heights$heights/100)^summary(rep8)[1,1])
  }
  if (model == 8){
    for (i in 1:nrow(new_long_heights)){
      temp_cov<-rep9$cov[c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[new_long_heights$tow_id][i],which(names(rep9$value)=="beta_b_s")[new_long_heights$tow_id][i]),
                         c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[new_long_heights$tow_id][i],which(names(rep9$value)=="beta_b_s")[new_long_heights$tow_id][i])]
      new_long_heights$spat_both[i]<-exp((summary(rep9)[6,1])*log(new_long_heights$depth[i])+Report9$beta_s[new_long_heights$tow_id][i]+((Report9$beta_b_s[new_long_heights$tow_id])[i]+(summary(rep9)[5,1]))*log(new_long_heights$heights[i]/100))
      new_long_heights_se$spat_both[i]<-sqrt(sum(c(temp_cov[1,1]*log(new_long_heights$depth[i])^2,
                                                   temp_cov[4,4]*log(new_long_heights$heights[i]/100)^2,
                                                   temp_cov[2,2]*log(new_long_heights$heights[i]/100)^2,
                                                   temp_cov[3,3],
                                                   2*temp_cov[1,4]*log(new_long_heights$depth[i])*log(new_long_heights$heights[i]/100),
                                                   2*temp_cov[1,2]*log(new_long_heights$depth[i])*log(new_long_heights$heights[i]/100),
                                                   2*temp_cov[2,4]*log(new_long_heights$heights[i]/100)^2,
                                                   2*temp_cov[1,3]*log(new_long_heights$depth[i]),
                                                   2*temp_cov[2,3]*log(new_long_heights$heights[i]/100),
                                                   2*temp_cov[3,4]*log(new_long_heights$heights[i]/100)))) 
    }
  }
}

blep<-data.frame(heights=mean(new_long_heights$heights),depth=mean(new_long_heights$depth),
                 tow_id=c(1:234)[-which(1:234 %in% unique(new_long_heights$tow_id))],
                 pred_heights=NA,
                 cut_heights=NA,off=0,off_ln=0,insh=0,spat=0,
                 spat_off=0,off_est=0,off_ln_est=0,spat_both=0)
new_long_heights<-rbind(new_long_heights,blep)
new_long_heights_se<-rbind(new_long_heights_se,blep)

agg_indiv_scallops<-data.frame(tow_id=1:234,off=NA,off_ln=NA,insh=NA,spat=NA,spat_off=NA,
                               off_est=NA,off_ln_est=NA,spat_both=NA)
agg_indiv_scallops_se<-data.frame(tow_id=1:234,off=NA,off_ln=NA,insh=NA,spat=NA,spat_off=NA,
                               off_est=NA,off_ln_est=NA,spat_both=NA)
for (tow in 1:234){
  for (model in 1:8){
    if (model %in% c(3,4,5,8)){
      estims=estimateSumLognormal(mu=log(new_long_heights[which(new_long_heights$tow_id==tow),model+5][which(new_long_heights[which(new_long_heights$tow_id==tow),model+5]>0)]),
                                  sigma=new_long_heights_se[which(new_long_heights_se$tow_id==tow),model+5][which(!is.na(new_long_heights_se[which(new_long_heights_se$tow_id==tow),model+5]))]^2)
      agg_indiv_scallops[tow,model+1]<-exp(estims[1])
      agg_indiv_scallops_se[tow,model+1]<-sqrt(estims[2])
    }
  }
}

agg_indiv_scallops$off<-aggregate(off~tow_id,data=new_long_heights,FUN=sum)$off
agg_indiv_scallops$off_ln<-aggregate(off_ln~tow_id,data=new_long_heights,FUN=sum)$off_ln
agg_indiv_scallops$off_est<-aggregate(off_est~tow_id,data=new_long_heights,FUN=sum)$off_est
agg_indiv_scallops$off_ln_est<-aggregate(off_ln_est~tow_id,data=new_long_heights,FUN=sum)$off_ln_est

agg_indiv_scallops_se$off<-aggregate(off~tow_id,data=new_long_heights_se,FUN=sum)$off
agg_indiv_scallops_se$off_ln<-aggregate(off_ln~tow_id,data=new_long_heights_se,FUN=sum,na.action=na.omit)$off_ln
agg_indiv_scallops_se$off_est<-aggregate(off_est~tow_id,data=new_long_heights_se,FUN=sum,na.action=na.omit)$off_est
agg_indiv_scallops_se$off_ln_est<-aggregate(off_ln_est~tow_id,data=new_long_heights_se,FUN=sum,na.action=na.omit)$off_ln_est


#Only predicting on predicted mean shell heights

small_heights<-long_heights[!duplicated(long_heights$tow_id),c(2,3,5)]
tow_ns<-aggregate(heights~tow_id,data=long_heights,FUN=length)
small_heights<-left_join(tow_ns,small_heights)
colnames(small_heights)[2]<-"n"

heights_var<-data.frame(tow_id=long_heights$tow_id[unique_means],
                       pred_height=Report5$unique_mu/10,
                       var=(rep5$sd[which(names(rep5$value)=="unique_mu")]^2)/100)
heights_var<-heights_var[order(heights_var$tow_id),]



all_small_heights<-small_heights
all_small_heights_se<-small_heights
all_small_heights_comb_se<-small_heights

all_small_heights_comb_se<-left_join(all_small_heights_comb_se,heights_var[,c(1,3)],by="tow_id")
all_small_heights_comb_se$cubed_var<-(9*all_small_heights_comb_se$pred_heights^4*all_small_heights_comb_se$var)+(36*all_small_heights_comb_se$pred_heights^2*all_small_heights_comb_se$var^2)+(15*sqrt(all_small_heights_comb_se$var)^6)
all_small_heights_comb_se$cubed_pred_heights<-all_small_heights_comb_se$pred_heights^3 + 3*all_small_heights_comb_se$pred_heights*all_small_heights_comb_se$var

match_gam_gb2<-match_gam_gb[which((1:234 %in% small_heights$tow_id))]

for (model in 1:length(mw_models)){
  if (model == 1) {
    # all_small_heights$off<-((surv.dat$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)
    all_small_heights$off<-((surv.dat$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)+(0.5*surv.dat$CFh[match_gam_gb2]*2*3*all_small_heights$pred_heights*heights_var$var)
    all_small_heights_se$off<-sqrt((surv.dat$se[match_gam_gb2]^2)*((all_small_heights$pred_heights)^3)^2)
    all_small_heights_comb_se$off<-sqrt(surv.dat$CFh[match_gam_gb2]^2 *all_small_heights_comb_se$cubed_var + all_small_heights_comb_se$cubed_var*surv.dat$se[match_gam_gb2]^2 +surv.dat$se[match_gam_gb2]^2 *all_small_heights_comb_se$cubed_pred_heights^2)
  }  
  if (model == 2) {
    # all_small_heights$off_ln<-((surv.dat2$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)
    all_small_heights$off_ln<-((surv.dat2$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)+(0.5*surv.dat2$CFh[match_gam_gb2]*2*3*all_small_heights$pred_heights*heights_var$var)
    all_small_heights_se$off_ln<-sqrt((surv.dat2$se[match_gam_gb2]^2)*((all_small_heights$pred_heights)^3)^2)
    all_small_heights_comb_se$off_ln<-sqrt(surv.dat2$CFh[match_gam_gb2]^2 *all_small_heights_comb_se$cubed_var + all_small_heights_comb_se$cubed_var*surv.dat2$se[match_gam_gb2]^2 +surv.dat2$se[match_gam_gb2]^2 *all_small_heights_comb_se$cubed_pred_heights^2)
  }  
  if (model == 3) {
    # all_small_heights$insh<-(exp((summary(rep3)[2,1])*log(all_small_heights$depth)+(summary(rep3)[1,1])*log(all_small_heights$pred_heights)))
    all_small_heights$insh<-(exp((summary(rep3)[2,1])*log(all_small_heights$depth)+(summary(rep3)[1,1])*log(all_small_heights$pred_heights)))+(0.5*exp((summary(rep3)[2,1])*log(all_small_heights$depth))*(summary(rep3)[1,1]-1)*(summary(rep3)[1,1])*exp((summary(rep3)[1,1]-2)*log(all_small_heights$pred_heights))*heights_var$var)
    all_small_heights_se$insh<-sqrt(rowSums(data.frame(c(summary_rep3_cov[1,1] * log(all_small_heights$depth)^2),
                                                       c(summary_rep3_cov[2,2] *log(all_small_heights$pred_heights)^2),
                                                       c(2*(summary_rep3_cov[1,2]* log(all_small_heights$depth)*log(all_small_heights$pred_heights))))))
    var_bheight<-summary_rep3_cov[2,2] * log(all_small_heights$pred_heights)^2 + summary_rep3_cov[2,2]*(all_small_heights$pred_heights^(-2)*all_small_heights_comb_se$var)+ (all_small_heights$pred_heights^(-2)*all_small_heights_comb_se$var)*summary(rep3)[1,1]^2
    all_small_heights_comb_se$insh<-sqrt(rowSums(data.frame(c(summary_rep3_cov[1,1] * log(all_small_heights$depth)^2),
                                                       c(var_bheight),
                                                       c(2*(summary_rep3_cov[1,2]* log(all_small_heights$depth)*log(all_small_heights$pred_heights))))))
  #Correlation between 2 parameters is -0.15, so I'm basically gonna ignore its impact and leave it as is because it does very very little
  }
  if (model == 4) {
    for (i in 1:nrow(all_small_heights)){
      temp_cov<-rep4$cov[c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[all_small_heights$tow_id][i]),
                         c(which(names(rep4$value)=="beta_depth")[1],which(names(rep4$value)=="beta")[2],which(names(rep4$value)=="beta_s")[all_small_heights$tow_id][i])]
      # all_small_heights$spat[i]<-(exp((summary(rep4)[4,1])*log(all_small_heights$depth[i])+(Report4$beta_s[all_small_heights$tow_id][i])+(summary(rep4)[3,1])*log(all_small_heights$pred_heights[i])))
      all_small_heights$spat[i]<-(exp((summary(rep4)[4,1])*log(all_small_heights$depth[i])+(Report4$beta_s[all_small_heights$tow_id][i])+(summary(rep4)[3,1])*log(all_small_heights$pred_heights[i])))+(0.5*exp((summary(rep4)[4,1])*log(all_small_heights$depth[i])+(Report4$beta_s[all_small_heights$tow_id][i]))*(summary(rep4)[3,1]-1)*summary(rep4)[3,1]*exp((summary(rep4)[3,1]-2)*log(all_small_heights$pred_heights[i]))*heights_var$var[i])
      var_bheight<-temp_cov[2,2] * log(all_small_heights$pred_heights[i])^2 + temp_cov[2,2]*(all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])+ (all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])*summary(rep4)[3,1]^2
      all_small_heights_se$spat[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2,
                                               temp_cov[3,3],
                                               temp_cov[2,2]*log(all_small_heights$pred_heights[i])^2,
                                               2*temp_cov[1,3]*log(all_small_heights$depth[i]),
                                               2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                               2*temp_cov[2,3]*log(all_small_heights$pred_heights[i]))))
      all_small_heights_comb_se$spat[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2,
                                            temp_cov[3,3],
                                            var_bheight,
                                            2*temp_cov[1,3]*log(all_small_heights$depth[i]),
                                            2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                            2*temp_cov[2,3]*log(all_small_heights$pred_heights[i]))))
    #Correlation between beta (b) and other parameters is negligible, so ima leave those ones as is
    }
  }
  if (model == 5) {
    for (i in 1:nrow(all_small_heights)){
      temp_cov<-rep6$cov[c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[all_small_heights$tow_id][i]),
                         c(which(names(rep6$value)=="beta_depth")[1],which(names(rep6$value)=="beta")[2],which(names(rep6$value)=="beta_s")[all_small_heights$tow_id][i])]
      # all_small_heights$spat_off[i]<-(exp((summary(rep6)[4,1])*log(all_small_heights$depth[i])+((Report6$beta_s[all_small_heights$tow_id][i])+(summary(rep6)[3,1]))*log(all_small_heights$pred_heights[i]))) 
      all_small_heights$spat_off[i]<-(exp((summary(rep6)[4,1])*log(all_small_heights$depth[i])+((Report6$beta_s[all_small_heights$tow_id][i])+(summary(rep6)[3,1]))*log(all_small_heights$pred_heights[i]))) + (0.5*exp(summary(rep6)[4,1]*log(all_small_heights$depth[i]))*((Report6$beta_s[all_small_heights$tow_id][i])+(summary(rep6)[3,1])-1)*((Report6$beta_s[all_small_heights$tow_id][i])+(summary(rep6)[3,1]))*exp(((Report6$beta_s[all_small_heights$tow_id][i])+(summary(rep6)[3,1])-2)*log(all_small_heights$pred_heights[i]))*heights_var$var[i])
      var_field_height<-temp_cov[3,3] * log(all_small_heights$pred_heights[i])^2 + temp_cov[3,3]*(all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])+ (all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])*rep6$value[which(names(rep6$value)=="beta_s")][all_small_heights$tow_id][i]^2
      var_bheight<-temp_cov[2,2] * log(all_small_heights$pred_heights[i])^2 + temp_cov[2,2]*(all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])+ (all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])*summary(rep6)[3,1]^2
      all_small_heights_se$spat_off[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2),
                                                 c(temp_cov[2,2]*log(all_small_heights$pred_heights[i])^2),
                                                 c(temp_cov[3,3]*log(all_small_heights$pred_heights[i])^2),
                                                 c(2*temp_cov[1,3]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i])),
                                                 c(2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i])),
                                                 c(2*temp_cov[2,3]*log(all_small_heights$pred_heights[i])^2))) 
      all_small_heights_comb_se$spat_off[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2),
                                                c(var_field_height),
                                                c(var_bheight),
                                                c(2*temp_cov[1,3]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i])),
                                                c(2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i])),
                                                c(2*temp_cov[2,3]*log(all_small_heights$pred_heights[i])^2))) 
    }
  }
  if (model == 6) {
    # all_small_heights$off_est<-((surv.dat7$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep7)[1,1]))
    all_small_heights$off_est<-((surv.dat7$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep7)[1,1]))+(0.5*(surv.dat7$CFh[match_gam_gb2])*(summary(rep7)[1,1]-1)*summary(rep7)[1,1]*(all_small_heights$pred_heights)^(summary(rep7)[1,1]-2)*heights_var$var)
    all_small_heights_se$off_est<-sqrt((surv.dat7$se[match_gam_gb2]^2 * (all_small_heights$pred_heights)^summary(rep7)[1,1]))
    all_small_heights_comb_se$off_est<-NA
  }
  if (model == 7) {
    # all_small_heights$off_ln_est<-((surv.dat8$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep8)[1,1]))
    all_small_heights$off_ln_est<-((surv.dat8$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep8)[1,1]))+(0.5*(surv.dat8$CFh[match_gam_gb2])*(summary(rep8)[1,1]-1)*summary(rep8)[1,1]*(all_small_heights$pred_heights)^(summary(rep8)[1,1]-2)*heights_var$var)
    all_small_heights_se$off_ln_est<-sqrt((surv.dat8$se[match_gam_gb2]^2 * (all_small_heights$pred_heights)^summary(rep8)[1,1]))
    all_small_heights_comb_se$off_ln_est<-NA
  }    
  if (model == 8){
    for (i in 1:nrow(all_small_heights)){
      temp_cov<-rep9$cov[c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[all_small_heights$tow_id][i],which(names(rep9$value)=="beta_b_s")[all_small_heights$tow_id][i]),
                         c(which(names(rep9$value)=="beta_depth")[1],which(names(rep9$value)=="beta")[2],which(names(rep9$value)=="beta_s")[all_small_heights$tow_id][i],which(names(rep9$value)=="beta_b_s")[all_small_heights$tow_id][i])]
      var_field_height<-temp_cov[4,4] * log(all_small_heights$pred_heights[i])^2 + temp_cov[4,4]*(all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])+ (all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])*rep9$value[which(names(rep9$value)=="beta_b_s")][all_small_heights$tow_id][i]^2
      var_bheight<-temp_cov[2,2] * log(all_small_heights$pred_heights[i])^2 + temp_cov[2,2]*(all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])+ (all_small_heights$pred_heights[i]^(-2)*all_small_heights_comb_se$var[i])*summary(rep9)[5,1]^2
      # all_small_heights$spat_both[i]<-exp((summary(rep9)[6,1])*log(all_small_heights$depth[i])+Report9$beta_s[all_small_heights$tow_id][i]+((Report9$beta_b_s[all_small_heights$tow_id])[i]+(summary(rep9)[5,1]))*log(all_small_heights$pred_heights[i]))
      all_small_heights$spat_both[i]<-exp((summary(rep9)[6,1])*log(all_small_heights$depth[i])+Report9$beta_s[all_small_heights$tow_id][i]+((Report9$beta_b_s[all_small_heights$tow_id])[i]+(summary(rep9)[5,1]))*log(all_small_heights$pred_heights[i]))+(0.5*exp((summary(rep9)[6,1])*log(all_small_heights$depth[i])+Report9$beta_s[all_small_heights$tow_id][i])*((Report9$beta_b_s[all_small_heights$tow_id])[i]+(summary(rep9)[5,1])-1)*((Report9$beta_b_s[all_small_heights$tow_id])[i]+(summary(rep9)[5,1]))*exp(((Report9$beta_b_s[all_small_heights$tow_id])[i]+(summary(rep9)[5,1])-2)*log(all_small_heights$pred_heights[i]))*heights_var$var[i])
      all_small_heights_se$spat_both[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2,
                                                   temp_cov[4,4]*log(all_small_heights$pred_heights[i])^2,
                                                   temp_cov[2,2]*log(all_small_heights$pred_heights[i])^2,
                                                   temp_cov[3,3],
                                                   2*temp_cov[1,4]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                                   2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                                   2*temp_cov[2,4]*log(all_small_heights$pred_heights[i])^2,
                                                   2*temp_cov[1,3]*log(all_small_heights$depth[i]),
                                                   2*temp_cov[2,3]*log(all_small_heights$pred_heights[i]),
                                                   2*temp_cov[3,4]*log(all_small_heights$pred_heights[i]))))
      all_small_heights_comb_se$spat_both[i]<-sqrt(sum(c(temp_cov[1,1]*log(all_small_heights$depth[i])^2,
                                                         c(var_field_height),
                                                         c(var_bheight),
                                                         temp_cov[3,3],
                                                         2*temp_cov[1,4]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                                         2*temp_cov[1,2]*log(all_small_heights$depth[i])*log(all_small_heights$pred_heights[i]),
                                                         2*temp_cov[2,4]*log(all_small_heights$pred_heights[i])^2,
                                                         2*temp_cov[1,3]*log(all_small_heights$depth[i]),
                                                         2*temp_cov[2,3]*log(all_small_heights$pred_heights[i]),
                                                         2*temp_cov[3,4]*log(all_small_heights$pred_heights[i]))))
    }
  }
} 

final_small_heights<-all_small_heights
final_small_heights_estim_se<-all_small_heights_se
final_small_heights_se<-all_small_heights_se
final_small_heights_comb<-all_small_heights
final_small_heights_comb_se<-all_small_heights_comb_se
for (col in 5:12){
  if (col %in% c(7:9,12)){
    for (i in 1:nrow(final_small_heights)){
      estims<-estimateSumLognormal(mu=rep(log(final_small_heights[i,col]),final_small_heights$n[i]),
                                   sigma=rep(final_small_heights_se[i,col],final_small_heights$n[i])^2)
      final_small_heights[i,col]<-exp(estims[1])
      final_small_heights_se[i,col]<-sqrt(estims[2])
      estims2<-estimateSumLognormal(mu=rep(log(final_small_heights[i,col]),final_small_heights$n[i]),
                                   sigma=rep(final_small_heights_comb_se[i,col+3],final_small_heights$n[i])^2)
      final_small_heights_comb[i,col]<-exp(estims2[1])
      final_small_heights_comb_se[i,col+3]<-sqrt(estims2[2])
    }
  } else {
   final_small_heights[,col]<-final_small_heights[,col]*final_small_heights$n
   final_small_heights_se[,col]<-sqrt((final_small_heights_se[,col])^2*final_small_heights$n^2)
   if (col %in% c(5:6)) final_small_heights_comb_se[,col+3]<-sqrt((final_small_heights_comb_se[,col+3])^2*final_small_heights$n^2)
  }
}

final_small_heights<-left_join(data.frame(tow_id=1:234),final_small_heights)
final_small_heights[is.na(final_small_heights)]<-0
final_small_heights_se<-left_join(data.frame(tow_id=1:234),final_small_heights_se)
final_small_heights_comb_se<-left_join(data.frame(tow_id=1:234),final_small_heights_comb_se)


#Got them all together, now just need to plot them
all_sums<-data.frame(tow=rep(agg_midpoints_scallops$tow_id,length(mw_models)),
                     val=c(agg_midpoints_scallops$off/1000,agg_midpoints_scallops$off_ln/1000,
                           agg_midpoints_scallops$insh/1000,agg_midpoints_scallops$spat/1000,
                           agg_midpoints_scallops$spat_off/1000,agg_midpoints_scallops$off_est/1000,
                           agg_midpoints_scallops$off_ln_est/1000,agg_midpoints_scallops$spat_both/1000),
                     sigma=c(sqrt(agg_midpoints_scallops_se$off^2/1000^2),sqrt(agg_midpoints_scallops_se$off_ln^2/1000^2),
                             agg_midpoints_scallops_se$insh,agg_midpoints_scallops_se$spat,
                             agg_midpoints_scallops_se$spat_off,sqrt(agg_midpoints_scallops_se$off_est^2/1000^2),
                             sqrt(agg_midpoints_scallops_se$off_ln_est^2/1000^2),agg_midpoints_scallops_se$spat_both),
                     high=c((agg_midpoints_scallops$off/1000)+1.96*sqrt(agg_midpoints_scallops_se$off^2/1000^2),
                            (agg_midpoints_scallops$off_ln/1000)+1.96*sqrt(agg_midpoints_scallops_se$off_ln^2/1000^2),
                            exp(log(agg_midpoints_scallops$insh)-log(1000)+1.96*agg_midpoints_scallops_se$insh),
                            exp(log(agg_midpoints_scallops$spat)-log(1000)+1.96*agg_midpoints_scallops_se$spat),
                            exp(log(agg_midpoints_scallops$spat_off)-log(1000)+1.96*agg_midpoints_scallops_se$spat_off),
                            (agg_midpoints_scallops$off_est/1000)+1.96*sqrt(agg_midpoints_scallops_se$off_est^2/1000^2),
                            (agg_midpoints_scallops$off_ln_est/1000)+1.96*sqrt(agg_midpoints_scallops_se$off_ln_est^2/1000^2),
                            exp(log(agg_midpoints_scallops$spat_both)-log(1000)+1.96*agg_midpoints_scallops_se$spat_both)),
                     low=c((agg_midpoints_scallops$off/1000)-1.96*sqrt(agg_midpoints_scallops_se$off^2/1000^2),
                           (agg_midpoints_scallops$off_ln/1000)-1.96*sqrt(agg_midpoints_scallops_se$off_ln^2/1000^2),
                           exp(log(agg_midpoints_scallops$insh)-log(1000)-1.96*agg_midpoints_scallops_se$insh),
                           exp(log(agg_midpoints_scallops$spat)-log(1000)-1.96*agg_midpoints_scallops_se$spat),
                           exp(log(agg_midpoints_scallops$spat_off)-log(1000)-1.96*agg_midpoints_scallops_se$spat_off),
                           (agg_midpoints_scallops$off_est/1000)-1.96*sqrt(agg_midpoints_scallops_se$off_est^2/1000^2),
                           (agg_midpoints_scallops$off_ln_est/1000)-1.96*sqrt(agg_midpoints_scallops_se$off_ln_est^2/1000^2),
                           exp(log(agg_midpoints_scallops$spat_both)-log(1000)-1.96*agg_midpoints_scallops_se$spat_both)),
                     model=rep(mw_models,each=nrow(agg_midpoints_scallops)),
                     step2=rep("Midpoint",nrow(agg_midpoints_scallops)*length(mw_models)),
                     geometry=rep(gb_2023_surv$geometry,length(mw_models)))

all_sums<-rbind(all_sums,
                data.frame(tow=rep(agg_indiv_scallops$tow_id,length(mw_models)),
                           val=c(agg_indiv_scallops$off/1000,agg_indiv_scallops$off_ln/1000,
                                 agg_indiv_scallops$insh/1000,agg_indiv_scallops$spat/1000,
                                 agg_indiv_scallops$spat_off/1000,agg_indiv_scallops$off_est/1000,
                                 agg_indiv_scallops$off_ln_est/1000,agg_indiv_scallops$spat_both/1000),
                           sigma=c(sqrt(agg_indiv_scallops_se$off^2/1000^2),sqrt(agg_indiv_scallops_se$off_ln^2/1000^2),
                                   agg_indiv_scallops_se$insh,agg_indiv_scallops_se$spat,
                                   agg_indiv_scallops_se$spat_off,sqrt(agg_indiv_scallops_se$off_est^2/1000^2),
                                   sqrt(agg_indiv_scallops_se$off_ln_est^2/1000^2),agg_indiv_scallops_se$spat_both),
                           high=c((agg_indiv_scallops$off/1000)+1.96*sqrt(agg_indiv_scallops_se$off^2/1000^2),
                                  (agg_indiv_scallops$off_ln/1000)+1.96*sqrt(agg_indiv_scallops_se$off_ln^2/1000^2),
                                  exp(log(agg_indiv_scallops$insh)-log(1000)+1.96*agg_indiv_scallops_se$insh),
                                  exp(log(agg_indiv_scallops$spat)-log(1000)+1.96*agg_indiv_scallops_se$spat),
                                  exp(log(agg_indiv_scallops$spat_off)-log(1000)+1.96*agg_indiv_scallops_se$spat_off),
                                  (agg_indiv_scallops$off_est/1000)+1.96*sqrt(agg_indiv_scallops_se$off_est^2/1000^2),
                                  (agg_indiv_scallops$off_ln_est/1000)+1.96*sqrt(agg_indiv_scallops_se$off_ln_est^2/1000^2),
                                  exp(log(agg_indiv_scallops$spat_both)-log(1000)+1.96*agg_indiv_scallops_se$spat_both)),
                           low=c((agg_indiv_scallops$off/1000)-1.96*sqrt(agg_indiv_scallops_se$off^2/1000^2),
                                 (agg_indiv_scallops$off_ln/1000)-1.96*sqrt(agg_indiv_scallops_se$off_ln^2/1000^2),
                                 exp(log(agg_indiv_scallops$insh)-log(1000)-1.96*agg_indiv_scallops_se$insh),
                                 exp(log(agg_indiv_scallops$spat)-log(1000)-1.96*agg_indiv_scallops_se$spat),
                                 exp(log(agg_indiv_scallops$spat_off)-log(1000)-1.96*agg_indiv_scallops_se$spat_off),
                                 (agg_indiv_scallops$off_est/1000)-1.96*sqrt(agg_indiv_scallops_se$off_est^2/1000^2),
                                 (agg_indiv_scallops$off_ln_est/1000)-1.96*sqrt(agg_indiv_scallops_se$off_ln_est^2/1000^2),
                                 exp(log(agg_indiv_scallops$spat_both)-log(1000)-1.96*agg_indiv_scallops_se$spat_both)),
                           model=rep(mw_models,each=nrow(agg_indiv_scallops)),
                           step2=rep("Individual",nrow(agg_indiv_scallops)*length(mw_models)),
                           geometry=rep(gb_2023_surv$geometry,length(mw_models))))

all_sums<-rbind(all_sums,
                data.frame(tow=rep(final_small_heights$tow_id,length(mw_models)),
                           val=c(final_small_heights$off/1000,final_small_heights$off_ln/1000,
                                 final_small_heights$insh/1000,final_small_heights$spat/1000,
                                 final_small_heights$spat_off/1000,final_small_heights$off_est/1000,
                                 final_small_heights$off_ln_est/1000,final_small_heights$spat_both/1000),
                           sigma=c(sqrt(final_small_heights_se$off^2/1000^2),sqrt(final_small_heights_se$off_ln^2/1000^2),
                                   final_small_heights_comb_se$insh,final_small_heights_comb_se$spat,
                                   final_small_heights_comb_se$spat_off,sqrt(final_small_heights_se$off_est^2/1000^2),
                                   sqrt(final_small_heights_se$off_ln_est^2/1000^2),final_small_heights_comb_se$spat_both),
                           high=c((final_small_heights$off/1000)+1.96*sqrt(final_small_heights_se$off^2/1000^2),
                                  (final_small_heights$off_ln/1000)+1.96*sqrt(final_small_heights_se$off_ln^2/1000^2),
                                  exp((log(final_small_heights$insh)-log(1000)+1.96*final_small_heights_comb_se$insh)),
                                  exp((log(final_small_heights$spat)-log(1000)+1.96*final_small_heights_comb_se$spat)),
                                  exp((log(final_small_heights$spat_off)-log(1000)+1.96*final_small_heights_comb_se$spat_off)),
                                  (final_small_heights$off_est/1000)+1.96*sqrt(final_small_heights_se$off_est^2/1000^2),
                                  (final_small_heights$off_ln_est/1000)+1.96*sqrt(final_small_heights_se$off_ln_est^2/1000^2),
                                  exp((log(final_small_heights$spat_both)-log(1000)+1.96*final_small_heights_comb_se$spat_both))),
                           low=c((final_small_heights$off/1000)-1.96*sqrt(final_small_heights_se$off^2/1000^2),
                                 (final_small_heights$off_ln/1000)-1.96*sqrt(final_small_heights_se$off_ln^2/1000^2),
                                 exp((log(final_small_heights$insh)-log(1000)-1.96*final_small_heights_comb_se$insh)),
                                 exp((log(final_small_heights$spat)-log(1000)-1.96*final_small_heights_comb_se$spat)),
                                 exp((log(final_small_heights$spat_off)-log(1000)-1.96*final_small_heights_comb_se$spat_off)),
                                 (final_small_heights$off_est/1000)-1.96*sqrt(final_small_heights_se$off_est^2/1000^2),
                                 (final_small_heights$off_ln_est/1000)-1.96*sqrt(final_small_heights_se$off_ln_est^2/1000^2),
                                 exp((log(final_small_heights$spat_both)-log(1000)-1.96*final_small_heights_comb_se$spat_both))),
                           model=rep(mw_models,each=nrow(final_small_heights)),
                           step2=rep("Mean Height",nrow(final_small_heights)*length(mw_models)),
                           geometry=rep(gb_2023_surv$geometry,length(mw_models)))) %>% 
  st_as_sf()

all_sums_sub<-subset(all_sums,model %in% mw_models[c(1,2,3,4,5,8)])

tows_dist<-st_drop_geometry(gb_2023_surv[,c(2,6,7,8,9)])

tow_start<-st_as_sf(tows_dist,coords=c("slon","slat"))
st_crs(tow_start)<-4326
tow_start<-tow_start %>% st_transform(crs=32619)
tow_end<-st_as_sf(tows_dist,coords=c("elon","elat"))
st_crs(tow_end)<-4326
tow_end<-tow_end %>% st_transform(crs=32619)

actual_dists<-data.frame(tow_id=gb_2023_surv$tow,dist=diag(st_distance(tow_start,tow_end)))

tot_area<-4252.073

#For lognormal, multiplying on natural scale is same as adding on log-scale
#Means that the variance doesn't change
transf_val<-((1/(as.numeric(actual_dists$dist)*2.4384/10^6))*tot_area)/1000000

transf_sums<-data.frame(tow=rep(agg_midpoints_scallops$tow_id,length(mw_models)),
                     val=c(agg_midpoints_scallops$off*transf_val,agg_midpoints_scallops$off_ln*transf_val,
                           agg_midpoints_scallops$insh*transf_val,agg_midpoints_scallops$spat*transf_val,
                           agg_midpoints_scallops$spat_off*transf_val,agg_midpoints_scallops$off_est*transf_val,
                           agg_midpoints_scallops$off_ln_est*transf_val,agg_midpoints_scallops$spat_both*transf_val),
                     sigma=c(sqrt(agg_midpoints_scallops_se$off^2*transf_val^2),sqrt(agg_midpoints_scallops_se$off_ln^2*transf_val^2),
                             agg_midpoints_scallops_se$insh,agg_midpoints_scallops_se$spat,agg_midpoints_scallops_se$spat_off,
                             sqrt(agg_midpoints_scallops_se$off_est^2*transf_val^2),
                             sqrt(agg_midpoints_scallops_se$off_ln_est^2*transf_val^2),
                             agg_midpoints_scallops_se$spat_both),
                     high=c((agg_midpoints_scallops$off*transf_val)+1.96*sqrt(agg_midpoints_scallops_se$off^2*transf_val^2),
                            (agg_midpoints_scallops$off_ln*transf_val)+1.96*sqrt(agg_midpoints_scallops_se$off_ln^2*transf_val^2),
                            exp(log(agg_midpoints_scallops$insh)+log(transf_val)+1.96*agg_midpoints_scallops_se$insh),
                            exp(log(agg_midpoints_scallops$spat)+log(transf_val)+1.96*agg_midpoints_scallops_se$spat),
                            exp(log(agg_midpoints_scallops$spat_off)+log(transf_val)+1.96*agg_midpoints_scallops_se$spat_off),
                            (agg_midpoints_scallops$off_est*transf_val)+1.96*sqrt(agg_midpoints_scallops_se$off_est^2*transf_val^2),
                            (agg_midpoints_scallops$off_ln_est*transf_val)+1.96*sqrt(agg_midpoints_scallops_se$off_ln_est^2*transf_val^2),
                            exp(log(agg_midpoints_scallops$spat_both)+log(transf_val)+1.96*agg_midpoints_scallops_se$spat_both)),
                     low=c((agg_midpoints_scallops$off*transf_val)-1.96*sqrt(agg_midpoints_scallops_se$off^2*transf_val^2),
                           (agg_midpoints_scallops$off_ln*transf_val)-1.96*(agg_midpoints_scallops_se$off_ln^2*transf_val^2),
                           exp(log(agg_midpoints_scallops$insh)+log(transf_val)-1.96*agg_midpoints_scallops_se$insh),
                           exp(log(agg_midpoints_scallops$spat)+log(transf_val)-1.96*agg_midpoints_scallops_se$spat),
                           exp(log(agg_midpoints_scallops$spat_off)+log(transf_val)-1.96*agg_midpoints_scallops_se$spat_off),
                           (agg_midpoints_scallops$off_est*transf_val)-1.96*(agg_midpoints_scallops_se$off_est^2*transf_val^2),
                           (agg_midpoints_scallops$off_ln_est*transf_val)-1.96*(agg_midpoints_scallops_se$off_ln_est^2*transf_val^2),
                           exp(log(agg_midpoints_scallops$spat_both)+log(transf_val)-1.96*agg_midpoints_scallops_se$spat_both)),
                     model=rep(mw_models,each=nrow(agg_midpoints_scallops)),
                     step2=rep("Midpoint",nrow(agg_midpoints_scallops)*length(mw_models)),
                     geometry=rep(gb_2023_surv$geometry,length(mw_models)))

transf_sums<-rbind(transf_sums,
                data.frame(tow=rep(agg_indiv_scallops$tow_id,length(mw_models)),
                           val=c(agg_indiv_scallops$off*transf_val,agg_indiv_scallops$off_ln*transf_val,
                                 agg_indiv_scallops$insh*transf_val,agg_indiv_scallops$spat*transf_val,
                                 agg_indiv_scallops$spat_off*transf_val,agg_indiv_scallops$off_est*transf_val,
                                 agg_indiv_scallops$off_ln_est*transf_val,agg_indiv_scallops$spat_both*transf_val),
                           sigma=c(sqrt(agg_indiv_scallops_se$off^2*transf_val^2),sqrt(agg_indiv_scallops_se$off_ln^2*transf_val^2),
                                   agg_indiv_scallops_se$insh,agg_indiv_scallops_se$spat,agg_indiv_scallops_se$spat_off,
                                   sqrt(agg_indiv_scallops_se$off_est^2*transf_val^2),sqrt(agg_indiv_scallops_se$off_ln_est^2*transf_val^2),
                                   agg_indiv_scallops_se$spat_both),
                           high=c((agg_indiv_scallops$off*transf_val)+1.96*sqrt(agg_indiv_scallops_se$off^2*transf_val^2),
                                  (agg_indiv_scallops$off_ln*transf_val)+1.96*sqrt(agg_indiv_scallops_se$off_ln^2*transf_val^2),
                                  exp(log(agg_indiv_scallops$insh)+log(transf_val)+1.96*agg_indiv_scallops_se$insh),
                                  exp(log(agg_indiv_scallops$spat)+log(transf_val)+1.96*agg_indiv_scallops_se$spat),
                                  exp(log(agg_indiv_scallops$spat_off)+log(transf_val)+1.96*agg_indiv_scallops_se$spat_off),
                                  (agg_indiv_scallops$off_est*transf_val)+1.96*sqrt(agg_indiv_scallops_se$off_est^2*transf_val^2),
                                  (agg_indiv_scallops$off_ln_est*transf_val)+1.96*sqrt(agg_indiv_scallops_se$off_ln_est^2*transf_val^2),
                                  exp(log(agg_indiv_scallops$spat_both)-log(transf_val)+1.96*agg_indiv_scallops_se$spat_both)),
                           low=c((agg_indiv_scallops$off*transf_val)-1.96*sqrt(agg_indiv_scallops_se$off^2*transf_val^2),
                                 (agg_indiv_scallops$off_ln*transf_val)-1.96*(agg_indiv_scallops_se$off_ln^2*transf_val^2),
                                 exp(log(agg_indiv_scallops$insh)+log(transf_val)-1.96*agg_indiv_scallops_se$insh),
                                 exp(log(agg_indiv_scallops$spat)+log(transf_val)-1.96*agg_indiv_scallops_se$spat),
                                 exp(log(agg_indiv_scallops$spat_off)+log(transf_val)-1.96*agg_indiv_scallops_se$spat_off),
                                 (agg_indiv_scallops$off_est*transf_val)-1.96*(agg_indiv_scallops_se$off_est^2*transf_val^2),
                                 (agg_indiv_scallops$off_ln_est*transf_val)-1.96*(agg_indiv_scallops_se$off_ln_est^2*transf_val^2),
                                 exp(log(agg_indiv_scallops$spat_both)+log(transf_val)-1.96*agg_indiv_scallops_se$spat_both)),
                           model=rep(mw_models,each=nrow(agg_indiv_scallops)),
                           step2=rep("Individual",nrow(agg_indiv_scallops)*length(mw_models)),
                           geometry=rep(gb_2023_surv$geometry,length(mw_models))))

transf_sums<-rbind(transf_sums,
                data.frame(tow=rep(final_small_heights$tow_id,length(mw_models)),
                           val=c(final_small_heights$off*transf_val,final_small_heights$off_ln*transf_val,
                                 final_small_heights$insh*transf_val,final_small_heights$spat*transf_val,
                                 final_small_heights$spat_off*transf_val,final_small_heights$off_est*transf_val,
                                 final_small_heights$off_ln_est*transf_val,final_small_heights$spat_both*transf_val),
                           sigma=c(sqrt(final_small_heights_se$off^2*transf_val^2),sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                   final_small_heights_comb_se$insh,final_small_heights_comb_se$spat,
                                   final_small_heights_comb_se$spat_off,sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                   sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),final_small_heights_comb_se$spat_both),
                           high=c((final_small_heights$off*transf_val)+1.96*sqrt(final_small_heights_se$off^2*transf_val^2),
                                  (final_small_heights$off*transf_val)+1.96*sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                  exp(log(final_small_heights$insh)+log(transf_val)+1.96*final_small_heights_comb_se$insh),
                                  exp(log(final_small_heights$spat)+log(transf_val)+1.96*final_small_heights_comb_se$spat),
                                  exp(log(final_small_heights$spat_off)+log(transf_val)+1.96*final_small_heights_comb_se$spat_off),
                                  (final_small_heights$off_est*transf_val)+1.96*sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                  (final_small_heights$off_ln_est*transf_val)+1.96*sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),
                                  exp(log(final_small_heights$spat_both)+log(transf_val)+1.96*final_small_heights_comb_se$spat_both)),
                           low=c((final_small_heights$off*transf_val)-1.96*sqrt(final_small_heights_se$off^2*transf_val^2),
                                 (final_small_heights$off*transf_val)-1.96*sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                 exp(log(final_small_heights$insh)+log(transf_val)-1.96*final_small_heights_comb_se$insh),
                                 exp(log(final_small_heights$spat)+log(transf_val)-1.96*final_small_heights_comb_se$spat),
                                 exp(log(final_small_heights$spat_off)+log(transf_val)-1.96*final_small_heights_comb_se$spat_off),
                                 (final_small_heights$off_est*transf_val)-1.96*sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                 (final_small_heights$off_ln_est*transf_val)-1.96*sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),
                                 exp(log(final_small_heights$spat_both)+log(transf_val)-1.96*final_small_heights_comb_se$spat_both)),
                           model=rep(mw_models,each=nrow(final_small_heights)),
                           step2=rep("Mean Height",nrow(final_small_heights)*length(mw_models)),
                           geometry=rep(gb_2023_surv$geometry,length(mw_models))))

transf_sums<-rbind(transf_sums,
                   data.frame(tow=rep(final_small_heights$tow_id,length(mw_models)),
                              val=c(final_small_heights$off*transf_val,final_small_heights$off_ln*transf_val,
                                    final_small_heights$insh*transf_val,final_small_heights$spat*transf_val,
                                    final_small_heights$spat_off*transf_val,final_small_heights$off_est*transf_val,
                                    final_small_heights$off_ln_est*transf_val,final_small_heights$spat_both*transf_val),
                              sigma=c(sqrt(final_small_heights_se$off^2*transf_val^2),sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                      final_small_heights_se$insh,final_small_heights_se$spat,
                                      final_small_heights_se$spat_off,sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                      sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),final_small_heights_se$spat_both),
                              high=c((final_small_heights$off*transf_val)+1.96*sqrt(final_small_heights_se$off^2*transf_val^2),
                                     (final_small_heights$off*transf_val)+1.96*sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                     exp(log(final_small_heights$insh)+log(transf_val)+1.96*final_small_heights_se$insh),
                                     exp(log(final_small_heights$spat)+log(transf_val)+1.96*final_small_heights_se$spat),
                                     exp(log(final_small_heights$spat_off)+log(transf_val)+1.96*final_small_heights_se$spat_off),
                                     (final_small_heights$off_est*transf_val)+1.96*sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                     (final_small_heights$off_ln_est*transf_val)+1.96*sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),
                                     exp(log(final_small_heights$spat_both)+log(transf_val)+1.96*final_small_heights_se$spat_both)),
                              low=c((final_small_heights$off*transf_val)-1.96*sqrt(final_small_heights_se$off^2*transf_val^2),
                                    (final_small_heights$off*transf_val)-1.96*sqrt(final_small_heights_se$off_ln^2*transf_val^2),
                                    exp(log(final_small_heights$insh)+log(transf_val)-1.96*final_small_heights_se$insh),
                                    exp(log(final_small_heights$spat)+log(transf_val)-1.96*final_small_heights_se$spat),
                                    exp(log(final_small_heights$spat_off)+log(transf_val)-1.96*final_small_heights_se$spat_off),
                                    (final_small_heights$off_est*transf_val)-1.96*sqrt(final_small_heights_se$off_est^2*transf_val^2),
                                    (final_small_heights$off_ln_est*transf_val)-1.96*sqrt(final_small_heights_se$off_ln_est^2*transf_val^2),
                                    exp(log(final_small_heights$spat_both)+log(transf_val)-1.96*final_small_heights_se$spat_both)),
                              model=rep(mw_models,each=nrow(final_small_heights)),
                              step2=rep("Mean Height (w/o Error Propagation)",nrow(final_small_heights)*length(mw_models)),
                              geometry=rep(gb_2023_surv$geometry,length(mw_models)))) %>% st_as_sf()

all_methods<-c("Midpoint","Individual","Mean Height","Mean Height (w/o Error Propagation)")

calc_tot<-data.frame(off=rep(NA,4),off_ln=rep(NA,4),
                     insh=rep(NA,4),spat=rep(NA,4),
                     spat_off=rep(NA,4),off_estim=rep(NA,4),
                     off_ln_estim=rep(NA,4),spat_both=rep(NA,4),
                     method=c("Midpoint","Individual","Mean Height","Mean Height (w/o Error Propagation)"))
calc_tot_se<-data.frame(off=rep(NA,4),off_ln=rep(NA,4),
                     insh=rep(NA,4),spat=rep(NA,4),
                     spat_off=rep(NA,4),off_estim=rep(NA,4),
                     off_ln_estim=rep(NA,4),spat_both=rep(NA,4),
                     method=c("Midpoint","Individual","Mean Height","Mean Height (w/o Error Propagation)"))

for (col in 1:8){
  for (method in 1:4){
    if (col %in% c(3:5,8)){
      temp_frame<-subset(transf_sums,model==mw_models[col] & step2==all_methods[method] & val>0)
      estims<-estimateSumLognormal(mu=log(temp_frame$val),sigma=temp_frame$sigma^2)
      calc_tot[method,col]<-(exp(estims[1])/234)
      calc_tot_se[method,col]<-sqrt(estims[2])
    } else {
      temp_frame<-subset(transf_sums,model==mw_models[col] & step2==all_methods[method])
      calc_tot[method,col]<-mean(temp_frame$val)
      calc_tot_se[method,col]<-sum(temp_frame$sigma^2,na.rm=T)/(234^2)
    }
  }
}

calc_tot<-as.matrix(calc_tot[,-9])
calc_tot_se<-as.matrix(calc_tot_se[,-9])

calc_tot_low<-calc_tot
calc_tot_high<-calc_tot
for (i in 1:8){
  if (i %in% c(3:5,8)){
    calc_tot_low[,i]<-exp(log(calc_tot[,i])-1.96*calc_tot_se[,i])
    calc_tot_high[,i]<-exp(log(calc_tot[,i])+1.96*calc_tot_se[,i])
  } else {
    calc_tot_low[,i]<-calc_tot[,i]-1.96*calc_tot_se[,i]
    calc_tot_high[,i]<-calc_tot[,i]+1.96*calc_tot_se[,i]
  }
}


tot_frame<-data.frame(model=rep(mw_models,4),
                      method=rep(all_methods,each=length(mw_models)),
                      low=c(calc_tot_low[1,],calc_tot_low[2,],calc_tot_low[3,],calc_tot_low[4,]),
                      est=c(calc_tot[1,],calc_tot[2,],calc_tot[3,],calc_tot[4,]),
                      high=c(calc_tot_high[1,],calc_tot_high[2,],calc_tot_high[3,],calc_tot_high[4,]))

library(forcats)
tot_frame$model<-as.factor(tot_frame$model)
tot_frame<-tot_frame %>% mutate(model=fct_relevel(model,mw_models))

tot_frame$size<-tot_frame$high-tot_frame$low

straight_mean_calc<-aggregate(val~model,data=subset(transf_sums,step2==all_methods[1]),FUN=sum)
straight_mean_calc$mean<-straight_mean_calc$val/234
temp_sums<-subset(transf_sums,step2==all_methods[1])
temp_sums$val[which(is.na(temp_sums$val))]<-0
straight_mean_calc$sd<-aggregate(val~model,data=temp_sums,FUN=sd)$val
straight_mean_calc$se<-straight_mean_calc$sd/sqrt(234)
straight_mean_calc$high<-straight_mean_calc$mean+1.96*straight_mean_calc$se
straight_mean_calc$low<-straight_mean_calc$mean-1.96*straight_mean_calc$se

