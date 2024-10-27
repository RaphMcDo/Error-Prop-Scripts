
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


load("./LW_Work/all_mw.RData")
gb_2023_mw<-all_mw %>% st_as_sf(coords=c("lon","lat"))
st_crs(gb_2023_mw)<-4326
gb_2023_mw<-st_transform(gb_2023_mw,crs=32619)
gb_2023_mw$geometry<-gb_2023_mw$geometry/1000

gb_2023_surv<-st_as_sf(bigbank[which(bigbank$year == 2023),],coords=c("lon","lat"))
gb_2023_surv<-subset(gb_2023_surv, bank="GBb")
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

#Residuals for all dem bois
resids_off<-sub_mw$wmw-Report$mu
resids_off_gam<-CF.fit$residuals

resids_off_ln<-log(sub_mw$wmw)-(Report2$mu-(exp(rep2$par.fixed[2])^2)/2)
resids_off_ln_gam<-CF.fit2$residuals

resids_in_ln<-log(sub_mw$wmw)-(Report3$mu-(exp(rep3$par.fixed[3])^2)/2)

resids_spat_in_ln<-log(sub_mw$wmw)-(Report4$mu-(exp(rep4$par.fixed[5])^2)/2)

resids_spat_height<-Report5$residuals

resids_spat_off_ln<-log(sub_mw$wmw)-(Report6$mu-(exp(rep6$par.fixed[5])^2)/2)

resids_off_estim<-sub_mw$wmw-Report7$mu
resids_off_estim_gam<-CF.fit7$residuals

resids_off_estim_ln<-log(sub_mw$wmw)-(Report8$mu-(exp(rep8$par.fixed[3])^2)/2)
resids_off_estim_ln_gam<-CF.fit8$residuals

resids_spat_both<-log(sub_mw$wmw)-(Report9$mu-exp(rep9$par.fixed[7])^2/2)

all_resids<-data.frame(resids=c(resids_off,resids_off_gam,
                                resids_off_ln,resids_off_ln_gam,
                                resids_in_ln,
                                resids_spat_in_ln,
                                resids_spat_height,
                                resids_spat_off_ln,
                                resids_spat_both,
                                resids_off_estim,resids_off_estim_gam,
                                resids_off_estim_ln,resids_off_estim_ln_gam),
                       model=as.factor(c(rep(1,length(resids_off)),rep(2,length(resids_off_gam)),
                                         rep(3,length(resids_off_ln)),rep(4,length(resids_off_ln_gam)),
                                         rep(5,length(resids_in_ln)),
                                         rep(6,length(resids_spat_in_ln)),
                                         rep(7,length(resids_spat_height)),
                                         rep(8,length(resids_spat_off_ln)),
                                         rep(9,length(resids_spat_both)),
                                         rep(10,length(resids_off_estim)),rep(11,length(resids_off_estim_gam)),
                                         rep(12,length(resids_off_estim_ln)),rep(13,length(resids_off_estim_ln_gam)))),
                       fitted=c(Report$mu,CF.fit$fitted.values,
                                Report2$mu,CF.fit2$fitted.values,
                                Report3$mu,
                                Report4$mu,
                                Report5$mu_sh,
                                Report6$mu,
                                Report9$mu,
                                Report7$mu,CF.fit7$fitted.values,
                                Report8$mu,CF.fit8$fitted.values),
                       geometry=c(sub_mw$geometry,sub_mw[!duplicated(sub_mw$ID),]$geometry,
                                  sub_mw$geometry,sub_mw[!duplicated(sub_mw$ID),]$geometry,
                                  sub_mw$geometry,
                                  sub_mw$geometry,
                                  gb_2023_surv$geometry[long_heights$tow_id],
                                  sub_mw$geometry,
                                  sub_mw$geometry,
                                  sub_mw$geometry,sub_mw[!duplicated(sub_mw$ID),]$geometry,
                                  sub_mw$geometry,sub_mw[!duplicated(sub_mw$ID),]$geometry),
                       ids=c(sub_mw$ID,sub_mw[!duplicated(sub_mw$ID),]$ID,
                             sub_mw$ID,sub_mw[!duplicated(sub_mw$ID),]$ID,
                             sub_mw$ID,
                             sub_mw$ID,
                             st_drop_geometry(gb_2023_surv)$ID[long_heights$tow_id],
                             sub_mw$ID,
                             sub_mw$ID,
                             sub_mw$ID,sub_mw[!duplicated(sub_mw$ID),]$ID,
                             sub_mw$ID,sub_mw[!duplicated(sub_mw$ID),]$ID)) %>% st_as_sf()

resid_names<-c("Off LMM","Off GAM",
               "Off LN LMM","Off LN GAM",
               "Depth",
               "Spatial",
               "Shell Height",
               "Spatial b",
               "Spatial Both",
               "Offshore Estimating b (GLMM)", "Offshore Estimating b (GAM)",
               "Offshore LN Estimating b (GLMM)", "Offshore LN Estimating b (GAM)")
resid_names2<-resid_names[c(1,3,5:9)]
resid_names3<-resid_names[c(1,3,5,6,8,9)]
mw_names<-resid_names[c(1:9)]

facet_labeller <- function(variable,value){
  return(resid_names[value])
}
facet_labeller2 <- function(variable,value){
  return(resid_names2[value])
}
facet_labeller3 <- function(variable,value){
  return(resid_names3[value])
}
mw_labeller<-function(variable,value){
  return(mw_names[value])
}

all_resids2<-subset(all_resids, model %in% c(1,3,5:9))
all_resids2$model<-as.factor(as.integer(all_resids2$model))
mw_resids<-subset(all_resids,model %in% c(1:6,8,9))

# big_resid_plot<-ggplot(data=all_resids,aes(sample=resids))+
# big_resid_plot<-ggplot(data=all_resids2,aes(sample=resids))+
# big_resid_plot<-ggplot(data=mw_resids,aes(sample=resids))+
big_resid_plot<-ggplot(data=subset(all_resids,model==7),aes(sample=resids))+
  stat_qq()+stat_qq_line()+
  ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+
  theme_bw()+
  facet_wrap(~model,scales="free",labeller=facet_labeller)
  # facet_wrap(~model,scales="free",labeller=facet_labeller2)
  # facet_wrap(~model,scales="free",labeller=mw_labeller)

# big_scedast_plot<-ggplot(data=all_resids,aes(y=resids,x=fitted))+
# big_scedast_plot<-ggplot(data=all_resids2,aes(y=resids,x=fitted))+
big_scedast_plot<-ggplot(data=mw_resids,aes(y=resids,x=fitted))+
# big_scedast_plot<-ggplot(data=subset(all_resids,model==7),aes(y=resids,x=fitted))+
  geom_point()+
  ylab("Residuals")+xlab("Fitted Values")+
  theme_bw()+
# facet_wrap(~model,scales="free",labeller=facet_labeller)
# facet_wrap(~model,scales="free",labeller=facet_labeller2)
  facet_wrap(~model,scales="free",labeller=mw_labeller)

mean_resids_tow<-aggregate(resids~model+ids,data=all_resids,FUN=mean)
mean_resids_tow2<-aggregate(resids~model+ids,data=all_resids,FUN=length) %>% rename(n=resids)
mean_resids_tow<-left_join(cbind(mean_resids_tow,mean_resids_tow2$n),all_resids[!duplicated(all_resids$ids),c("ids","geometry")],by="ids") %>% rename(n="mean_resids_tow2$n") %>% st_as_sf()
median_resids_tow<-aggregate(resids~model+ids,data=all_resids,FUN=median)
median_resids_tow2<-aggregate(resids~model+ids,data=all_resids,FUN=length) %>% rename(n=resids)
median_resids_tow<-left_join(cbind(median_resids_tow,median_resids_tow2$n),all_resids[!duplicated(all_resids$ids),c("ids","geometry")],by="ids") %>% rename(n="median_resids_tow2$n") %>% st_as_sf()
sub_mean_resids_tow<-aggregate(resids~model+ids,data=all_resids2,FUN=mean)
sub_mean_resids_tow2<-aggregate(resids~model+ids,data=all_resids2,FUN=length) %>% rename(n=resids)
sub_mean_resids_tow<-left_join(cbind(sub_mean_resids_tow,sub_mean_resids_tow2$n),all_resids[!duplicated(all_resids$ids),c("ids","geometry")],by="ids") %>% rename(n="sub_mean_resids_tow2$n")  %>% st_as_sf()
sub_median_resids_tow<-aggregate(resids~model+ids,data=all_resids2,FUN=median)
sub_median_resids_tow2<-aggregate(resids~model+ids,data=all_resids2,FUN=length) %>% rename(n=resids)
sub_median_resids_tow<-left_join(cbind(sub_median_resids_tow,sub_median_resids_tow2$n),all_resids[!duplicated(all_resids$ids),c("ids","geometry")],by="ids")%>% rename(n="sub_median_resids_tow2$n")  %>% st_as_sf()

spat_mean_resid<-ggplot()+
  # geom_sf(data=mean_resids_tow,aes(col=resids,size=n))+
  # geom_sf(data=sub_mean_resids_tow,aes(col=resids,size=n))+
  # geom_sf(data=subset(mean_resids_tow,model %in% c(1:6,8,9)),aes(col=resids,size=n))+
  geom_sf(data=subset(mean_resids_tow,model %in% c(7)),aes(col=resids,size=n))+
  scale_color_gradient2(name="Mean Residual",midpoint=mean(mean_resids_tow$resids,na.rm=T),
                        low="blue",mid="white",high="red")+
  scale_size(name="Number of observations")+
  facet_wrap(~model,labeller=facet_labeller)
  # facet_wrap(~model,labeller=facet_labeller2)
  # facet_wrap(~model,labeller=mw_labeller)

spat_median_resid<-ggplot()+
  # geom_sf(data=median_resids_tow,aes(col=resids,size=n))+
  geom_sf(data=sub_median_resids_tow,aes(col=resids,size=n))+
  scale_color_gradient2(name="median Residual",midpoint=mean(median_resids_tow$resids,na.rm=T),
                        low="blue",mid="white",high="red")+
  # facet_wrap(~model,labeller=facet_labeller)
  facet_wrap(~model,labeller=facet_labeller2)

#Will do CV for all of them as well, but will only show subset
library(caret)

#Height CV
test_heights<-list()
pred_heights<-list()
height_obj<-list()
height_Opt<-list()
height_rep<-list()
height_Report<-list()
height_MSPE<-list()

parameters = list(
  # beta = matrix(rep(0, nvar*2),ncol=2),
  beta_depth = c(0),
  beta_sh = c(14),
  beta_sh_s = rep(0,nrow(tmb_data$locations)),
  # beta_size_bin = matrix(c(seq(62.5,157.5,by=5)/10,
  #                          rep(0,length(seq(62.5,157.5,by=5)))),ncol=2),
  # log_rho = rep(log(10),2),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_upsilon = log(0.5))

#Does stratified on factors, so if I give it factors of tow, then should work
set.seed(921)
folds<-createFolds(as.factor(long_heights$tow_id),k=10,returnTrain = T)

for (fold in 1:10){

  test_heights[[fold]]<-long_heights[folds[[fold]],]
  pred_heights[[fold]]<-long_heights[-folds[[fold]],]

  tmb_data<-list(varmat_sh=as.matrix(cbind(rep(1,nrow(test_heights[[fold]])))),
                 depth_sh=test_heights[[fold]]$depth,
                 ind_loc_sh=test_heights[[fold]]$tow_id-1,
                 locations=st_coordinates(gb_2023_surv),
                 s_heights=test_heights[[fold]]$heights/10,
                 a_b_truncation=c(9.5,17))

  maps <- list(
    log_nu=c(factor(NA))#
  )

  height_obj[[fold]] = MakeADFun(data=tmb_data,
                   parameters=parameters,
                   map=maps,
                   random=c("beta_sh_s"),
                   DLL="spatial_depth_glmm_height",
                   silent = F)

  height_Opt[[fold]]<-optimx::optimr(height_obj[[fold]]$par,height_obj[[fold]]$fn,height_obj[[fold]]$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  height_rep[[fold]]<-sdreport(height_obj[[fold]],bias.correct=F)

  height_Report[[fold]]<-height_obj[[fold]]$report()

  upsilon<-exp(height_rep[[fold]]$value[which(names(height_rep[[fold]]$value)=="log_upsilon")])
  mus<-(height_Report[[fold]]$beta_sh + height_Report[[fold]]$beta_sh_s[pred_heights[[fold]]$tow_id]+height_Report[[fold]]$beta_depth*pred_heights[[fold]]$depth)

  pred_heights[[fold]]$pred_height<-mus +(((dnorm(9.5,mus,upsilon))-(dnorm(17,mus,upsilon)))/(pnorm(17,mus,upsilon)-pnorm(9.5,mus,upsilon)))*upsilon
  pred_heights[[fold]]$squared_error<-(pred_heights[[fold]]$heights/10-pred_heights[[fold]]$pred_height)^2

}

mw_models<-c(mw_names[c(1,3,5:6,8,9)],
             "Offshore Estimating b",
             "Offshore LN Estimating")#,
             #"Inshore Depth Tow LN GLMM")
mw_models[3]<-"Inshore Depth LN LMM"

tmb_models<-c("glmm_offshore_fixed_b",
              "glmm_offshore_fixed_b_lognormal",
              "depth_glmm_no_space",
              "spatial_depth_glmm_direct_cov",
              "spatial_depth_glmm_offshore",
              "spatial_both",
              "glmm_offshore_fixed_b",
              "glmm_offshore_fixed_b_lognormal")#,
              #"depth_glmm_tow_eff")

set.seed(0493)
folds<-createFolds(as.factor(sub_mw$ID),k=10,returnTrain = T)

# dyn.unload(dynlib("./LW_Work/depth_glmm_tow_eff"))
compile("./LW_Work/depth_glmm_tow_eff.cpp")
dyn.load(dynlib("./LW_Work/depth_glmm_tow_eff"))

test_weights<-list()
pred_weights<-list()
weight_obj<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
weight_Opt<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
weight_rep<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
weight_Report<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
weight_MSPE<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())

for (fold in 1:10){
  for (model in 1:length(mw_models)){

test_weights[[fold]]<-sub_mw[folds[[fold]],]
pred_weights[[fold]]<-sub_mw[-folds[[fold]],]

tmb_data<-list(weight=test_weights[[fold]]$wmw,
               heights=test_weights[[fold]]$sh/100,
               tow_id=as.integer(as.factor(test_weights[[fold]]$ID))-1,
               weight=test_weights[[fold]]$wmw,
               varmat=as.matrix(cbind(1, log(test_weights[[fold]]$sh/100))),
               depth=test_weights[[fold]]$depth,
               ind_loc=test_weights[[fold]]$id_loc2-1,
               locations=st_coordinates(gb_2023_surv))

    if (startsWith(mw_models[model],"Off")){
      parameters<-list(b=3,
                     beta=20,
                     log_phi=-1,
                     log_epsilon=-1,
                     tow_eff=rep(0,length(unique(sub_mw$ID))))

      random<-c("tow_eff")

      if (grepl("Estimating",mw_models[model],fixed=T)) maps<-list() else maps<-list(b=factor(NA))
    } else if (startsWith(mw_models[model],"Insh")){
      parameters = list(
        beta = rep(0, nvar),
        beta_depth = rep(0, nvar),
        log_phi = log(0.1))

      maps <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA))
      )

      random<-c()

    } else if (grepl("Both",mw_models[model],fixed=T)) {
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
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA),
        log_nu_b=factor(NA)
      )

      random=c("beta_s","beta_b_s")
    } else {
    parameters = list(
    beta = rep(0, nvar),
    beta_depth = rep(0, nvar),
    beta_s = rep(0,nrow(tmb_data$locations)),
    log_rho = log(10),
    log_sig = log(0.1),
    log_nu = log(1),
    log_phi = log(0.1),
    log_epsilon = log(0.1))#,
    # tow_eff = rep(0,length(unique(sub_mw$ID))))

    maps <- list(
      # beta = factor(rep(NA, nvar)),
      # beta = factor(c(2,NA)),
      beta = factor(c(NA,2)),
      # beta = factor(c(3,2)),
      beta_depth=factor(c(1,NA)),
      log_nu=factor(NA)
    )

    random<-c("beta_s")
    }

  weight_obj[[model]][[fold]] = MakeADFun(data=tmb_data,
                   parameters=parameters,
                   map=maps,
                   random=random,
                   DLL=tmb_models[model],
                   silent = F)

  weight_Opt[[model]][[fold]]<-optimx::optimr(weight_obj[[model]][[fold]]$par,weight_obj[[model]][[fold]]$fn,weight_obj[[model]][[fold]]$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")
  while (weight_Opt[[model]][[fold]]$message=="iteration limit reached without convergence (10)"){
    weight_obj[[model]][[fold]]$par<-weight_obj[[model]][[fold]]$env$last.par.best[-which(names(weight_obj[[model]][[fold]]$env$last.par.best)==random)]
    weight_Opt[[model]][[fold]]<-optimx::optimr(weight_obj[[model]][[fold]]$par,weight_obj[[model]][[fold]]$fn,weight_obj[[model]][[fold]]$gr,
                                                control=list(maxit=100000,maxeval=100000),
                                                method="nlminb")
  }

  weight_rep[[model]][[fold]]<-sdreport(weight_obj[[model]][[fold]],bias.correct=F)

  weight_Report[[model]][[fold]]<-weight_obj[[model]][[fold]]$report()

  }
}

save(test_weights,pred_weights,weight_obj,weight_Opt,weight_rep,weight_Report,file="./LW_Work/weight_cv_output.RData")

weight_predictions<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list(),in_tow_ln=list())

for (fold in 1:10){
  for (model in 1:length(mw_models)){
    weight_predictions[[model]][[fold]]<-pred_weights[[fold]]
    weight_predictions[[model]][[fold]]$sh<-weight_predictions[[model]][[fold]]$sh/100
    # weight_predictions[[model]][[fold]]$pred_weight<-rep(NA,nrow(weight_predictions[[model]][[fold]]))
    if (model==1) weight_predictions[[model]][[fold]]$pred_weight<-(summary(weight_rep[[model]][[fold]])[1,1]+weight_Report[[model]][[fold]]$tow_eff[as.integer(as.factor(pred_weights[[fold]]$ID))])*(weight_predictions[[model]][[fold]]$sh)^3
    if (model==2) weight_predictions[[model]][[fold]]$pred_weight<-exp(log((summary(weight_rep[[model]][[fold]])[1,1]+weight_Report[[model]][[fold]]$tow_eff[as.integer(as.factor(pred_weights[[fold]]$ID))])*(weight_predictions[[model]][[fold]]$sh)^3)-(exp(summary(weight_rep[[model]][[fold]])[2,1])^2)/2)
    if (model==3) weight_predictions[[model]][[fold]]$pred_weight<-exp((summary(weight_rep[[model]][[fold]])[2,1]*log(weight_predictions[[model]][[fold]]$depth) + summary(weight_rep[[model]][[fold]])[1,1]*log(weight_predictions[[model]][[fold]]$sh))-(exp(summary(weight_rep[[model]][[fold]])[3,1])^2)/2)
    if (model==4) weight_predictions[[model]][[fold]]$pred_weight<-exp((summary(weight_rep[[model]][[fold]])[4,1]*log(weight_predictions[[model]][[fold]]$depth) + weight_Report[[model]][[fold]]$beta_s[pred_weights[[fold]]$id_loc2] + summary(weight_rep[[model]][[fold]])[3,1]*log(weight_predictions[[model]][[fold]]$sh))-(exp(summary(weight_rep[[model]][[fold]])[5,1])^2)/2)
    if (model==5) weight_predictions[[model]][[fold]]$pred_weight<-exp((summary(weight_rep[[model]][[fold]])[4,1]*log(weight_predictions[[model]][[fold]]$depth) + (summary(weight_rep[[model]][[fold]])[3,1] + weight_Report[[model]][[fold]]$beta_s[pred_weights[[fold]]$id_loc2])*log(weight_predictions[[model]][[fold]]$sh))-(exp(summary(weight_rep[[model]][[fold]])[5,1])^2)/2)
    if (model==6) weight_predictions[[model]][[fold]]$pred_weight<-exp((summary(weight_rep[[model]][[fold]])[6,1]*log(weight_predictions[[model]][[fold]]$depth) + weight_Report[[model]][[fold]]$beta_s[pred_weights[[fold]]$id_loc2] +(summary(weight_rep[[model]][[fold]])[5,1]+weight_Report[[model]][[fold]]$beta_b_s[pred_weights[[fold]]$id_loc2])*log(weight_predictions[[model]][[fold]]$sh)))-(exp(summary(weight_rep[[model]][[fold]])[7,1])^2)/2
    if (model==7) weight_predictions[[model]][[fold]]$pred_weight<-(summary(weight_rep[[model]][[fold]])[2,1]+weight_Report[[model]][[fold]]$tow_eff[as.integer(as.factor(pred_weights[[fold]]$ID))])*(weight_predictions[[model]][[fold]]$sh)^summary(weight_rep[[model]][[fold]])[1,1]
    if (model==8) weight_predictions[[model]][[fold]]$pred_weight<-exp(log((summary(weight_rep[[model]][[fold]])[2,1]+weight_Report[[model]][[fold]]$tow_eff[as.integer(as.factor(pred_weights[[fold]]$ID))])*(weight_predictions[[model]][[fold]]$sh)^summary(weight_rep[[model]][[fold]])[1,1])-(exp(summary(weight_rep[[model]][[fold]])[3,1])^2)/2)
    if (model==9) weight_predictions[[model]][[fold]]$pred_weight<-exp((summary(weight_rep[[model]][[fold]])[2,1]*log(weight_predictions[[model]][[fold]]$depth) + summary(weight_rep[[model]][[fold]])[5:119,1][as.integer(as.factor(weight_predictions[[model]][[fold]]$ID))] + summary(weight_rep[[model]][[fold]])[1,1]*log(weight_predictions[[model]][[fold]]$sh))-(exp(summary(weight_rep[[model]][[fold]])[3,1])^2)/2)
    weight_predictions[[model]][[fold]]$squared_error<-(weight_predictions[[model]][[fold]]$wmw-weight_predictions[[model]][[fold]]$pred_weight)^2

  }
}

mw_models[3]<-"Depth"

all_pred_heights<-cbind(pred_heights[[1]],fold=rep(1,nrow(pred_heights[[1]])))
for (fold in 2:10){
  all_pred_heights<-rbind(all_pred_heights,cbind(pred_heights[[fold]],fold=rep(fold,nrow(pred_heights[[fold]]))))
}

all_pred_heights$geometry<-gb_2023_surv[all_pred_heights$tow_id,]$geometry
all_pred_heights<-st_as_sf(all_pred_heights)

for (model in 1:length(mw_models)){
  for (fold in 1:10){
    if (model == 1 && fold ==1) {
      all_pred_weights<-cbind(weight_predictions[[model]][[fold]],
                              model=rep(mw_models[model],nrow(weight_predictions[[model]][[fold]])),
                              fold=rep(fold,nrow(weight_predictions[[model]][[fold]])))
    } else {
      all_pred_weights<-rbind(all_pred_weights,
                              cbind(weight_predictions[[model]][[fold]],
                                    model=rep(mw_models[model],nrow(weight_predictions[[model]][[fold]])),
                                    fold=rep(fold,nrow(weight_predictions[[model]][[fold]]))))
    }
  } 
}

library(forcats)
all_pred_weights$model<-as.factor(all_pred_weights$model)
all_pred_weights<-all_pred_weights %>% mutate(model=fct_relevel(model,mw_models))

all_pred_heights$rsq_e<-sqrt(all_pred_heights$squared_error)/10
all_pred_weights$rsq_e<-sqrt(all_pred_weights$squared_error)
heights_overall_rmspe<-mean(all_pred_heights$rsq_e)
weights_overall_rmspe<-aggregate(rsq_e~model,data=all_pred_weights,FUN=mean)
heights_fold_rmspe<-aggregate(rsq_e~fold,data=all_pred_heights,FUN=mean)
weights_fold_rmspe<-aggregate(rsq_e~fold+model,data=all_pred_weights,FUN=mean)
weights_tow_rmspe<-aggregate(rsq_e~ID+model,data=all_pred_weights,FUN=mean)
weights_tow_rmspe2<-aggregate(rsq_e~ID+model,data=all_pred_weights,FUN=length) %>% rename(n=rsq_e)
weights_tow_rmspe<-cbind(weights_tow_rmspe,n=weights_tow_rmspe2$n)

weights_tow_rmspe<-left_join(as.data.frame(weights_tow_rmspe),gb_2023_surv[,c("geometry","ID")],by="ID") %>% st_as_sf()

weight_model_rmspe_plot<-ggplot()+
  geom_point(data=weights_overall_rmspe,aes(x=model,y=rsq_e),cex=2)+
  # geom_point(data=subset(weights_overall_rmspe,model %in% mw_models[1:5]),aes(x=model,y=rsq_e),cex=2)+
  ylab("RMSPE")+xlab("Fold")+
  theme_bw()


weights_overall_rmspe

weight_fold_rmspe_plot<-ggplot()+
  geom_point(data=weights_fold_rmspe,aes(x=fold,y=rsq_e,col=model),cex=2)+
  # geom_point(data=subset(weights_fold_rmspe,model %in% mw_models[1:6]),
  # aes(x=fold,y=rsq_e,shape=model,col=model),cex=2)+
  ylab("RMSPE")+xlab("Fold")+
  scale_color_viridis_d(name="Model")+
  scale_shape_discrete(name="Model")+
  theme_bw()

weight_tow_spat_rmspe_plot<-ggplot()+
  # geom_sf(data=weights_tow_rmspe,aes(col=rsq_e,size=n))+
  geom_sf(data=subset(weights_tow_rmspe,model %in% mw_models[1:6]),aes(col=rsq_e,size=n))+
  scale_color_viridis_c(name="RMSPE")+
  scale_size(name="Number of sampled scallops")+
  facet_wrap(~model)+
  theme_bw()

heights_tow_rmspe<-aggregate(rsq_e~tow_id,data=all_pred_heights,FUN=mean)
heights_tow_rmspe2<-aggregate(rsq_e~tow_id,data=all_pred_heights,FUN=length) %>% rename(n=rsq_e)
heights_tow_rmspe<-cbind(heights_tow_rmspe,n=heights_tow_rmspe2$n)

heights_tow_rmspe$geometry<-rep(NA,nrow(heights_tow_rmspe))
for (i in 1:nrow(heights_tow_rmspe)){
  heights_tow_rmspe$geometry[i]<-all_pred_heights$geometry[which(heights_tow_rmspe$tow_id[i]==all_pred_heights$tow_id)[1]]
}
heights_tow_rmspe<-st_as_sf(heights_tow_rmspe)

height_tow_spat_rmspe_plot<-ggplot()+
  geom_sf(data=heights_tow_rmspe,aes(col=rsq_e,size=n))+
  scale_color_viridis_c(name="RMSPE")+
  scale_size(name="Number of fully-recruited \nscallops caught")+
  theme_bw()


#Predicted mean heights per tow
a_b_truncation=c(9.5,17)
non_cor_sh_sf<-long_heights
non_cor_sh_sf$geometry<-gb_2023_surv$geometry[long_heights$tow_id]
non_cor_sh_sf<-st_as_sf(non_cor_sh_sf)
# st_crs(non_cor_sh_sf)<-4326
non_cor_sh_sf$mean_pred_height<-Report5$mu_sh+exp(summary(rep5)[5,1])*(dnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1]))-dnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1])))/(pnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1]))-pnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1])))
non_cor_sh_sf<-non_cor_sh_sf[order(non_cor_sh_sf$tow_id),]

sub_sh_sf<-non_cor_sh_sf[!duplicated(non_cor_sh_sf$tow_id),]
sub_sh_sf$n<-aggregate(mean_pred_height~tow_id,data=st_drop_geometry(non_cor_sh_sf),FUN=length)$mean_pred_height

spat_mean_pred_height<-ggplot()+
  geom_sf(data=sub_sh_sf,aes(col=mean_pred_height,size=n))+
  scale_color_viridis_c(name="Predicted Mean \nShell Height",
                        limits=c(9.5,15.2))+
  scale_size(name="Number of fully-recruited \nscallops caught")+
  theme_bw()

# #Raw commercial size mean shell heights
sub_sh_sf$mean_obs_height<-aggregate(heights~tow_id,data=st_drop_geometry(non_cor_sh_sf),FUN=mean)$heights

spat_mean_obs_height<-ggplot()+
  geom_sf(data=sub_sh_sf,aes(col=mean_obs_height/10,size=n))+
  scale_color_viridis_c(name="Observed Mean \nShell Height",
                        limits=c(9.5,15.2))+
  scale_size(name="Number of fully-recruited \nscallops caught")+
  theme_bw()

#Doing leave-1-loc-out CV, since we have to predict in unsampled areas
test_tcv_heights<-list()
pred_tcv_heights<-list()
tcv_height_obj<-list()
tcv_height_Opt<-list()
tcv_height_rep<-list()
tcv_height_Report<-list()
tcv_height_MSPE<-list()

parameters = list(
  # beta = matrix(rep(0, nvar*2),ncol=2),
  beta_depth = c(0),
  beta_sh = c(14),
  beta_sh_s = rep(0,nrow(gb_2023_surv)),
  # beta_size_bin = matrix(c(seq(62.5,157.5,by=5)/10,
  #                          rep(0,length(seq(62.5,157.5,by=5)))),ncol=2),
  # log_rho = rep(log(10),2),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_upsilon = log(0.5))


for (tow_fold in 1:length(unique(long_heights$tow_id))){

  test_tcv_heights[[tow_fold]]<-long_heights[-which(long_heights$tow_id==unique(long_heights$tow_id)[tow_fold]),]
  pred_tcv_heights[[tow_fold]]<-long_heights[which(long_heights$tow_id==unique(long_heights$tow_id)[tow_fold]),]

  tmb_data<-list(varmat_sh=as.matrix(cbind(rep(1,nrow(test_tcv_heights[[tow_fold]])))),
                 depth_sh=test_tcv_heights[[tow_fold]]$depth,
                 ind_loc_sh=test_tcv_heights[[tow_fold]]$tow_id-1,
                 locations=st_coordinates(gb_2023_surv),
                 s_heights=test_tcv_heights[[tow_fold]]$heights/10,
                 a_b_truncation=c(9.5,17))

  maps <- list(
    log_nu=c(factor(NA))#
  )

  tcv_height_obj[[tow_fold]] = MakeADFun(data=tmb_data,
                   parameters=parameters,
                   map=maps,
                   random=c("beta_sh_s"),
                   DLL="spatial_depth_glmm_height",
                   silent = F)

  tcv_height_Opt[[tow_fold]]<-optimx::optimr(tcv_height_obj[[tow_fold]]$par,tcv_height_obj[[tow_fold]]$fn,tcv_height_obj[[tow_fold]]$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  tcv_height_rep[[tow_fold]]<-sdreport(tcv_height_obj[[tow_fold]],bias.correct=F)

  tcv_height_Report[[tow_fold]]<-tcv_height_obj[[tow_fold]]$report()

  upsilon<-exp(tcv_height_rep[[tow_fold]]$value[which(names(tcv_height_rep[[tow_fold]]$value)=="log_upsilon")])
  mus<-(tcv_height_Report[[tow_fold]]$beta_sh + tcv_height_Report[[tow_fold]]$beta_sh_s[pred_tcv_heights[[tow_fold]]$tow_id]+tcv_height_Report[[tow_fold]]$beta_depth*pred_tcv_heights[[tow_fold]]$depth)

  pred_tcv_heights[[tow_fold]]$pred_height<-mus +(((dnorm(9.5,mus,upsilon))-(dnorm(17,mus,upsilon)))/(pnorm(17,mus,upsilon)-pnorm(9.5,mus,upsilon)))*upsilon
  pred_tcv_heights[[tow_fold]]$squared_error<-(pred_tcv_heights[[tow_fold]]$heights/10-pred_tcv_heights[[tow_fold]]$pred_height)^2

}

tot_height_tow_folds<-length(unique(long_heights$tow_id))

all_tow_pred_heights<-pred_tcv_heights[[1]]
for (i in 2:tot_height_tow_folds){
  all_tow_pred_heights<-rbind(all_tow_pred_heights,pred_tcv_heights[[i]])
}

tcv_height_rmspe<-sqrt(mean(all_tow_pred_heights$squared_error))/10

tcv_tow_height_rmspe<-aggregate(squared_error~tow_id,data=all_tow_pred_heights,FUN=mean)

tcv_tow_height_rmspe$geometry<-gb_2023_surv$geometry[tcv_tow_height_rmspe$tow_id]
tcv_tow_height_rmspe<-st_as_sf(tcv_tow_height_rmspe)
tcv_tow_height_rmspe$n<-aggregate(squared_error~tow_id,data=all_tow_pred_heights,FUN=length)$squared_error

tcv_height_tow_rmspe_plot<-ggplot()+
  geom_sf(data=tcv_tow_height_rmspe,aes(col=sqrt(squared_error),size=n))+
  scale_color_viridis_c(name="RMSPE")+
  scale_size(name="Number of scallops")+
  theme_bw()

tcv_height_tow_rel_rmspe_plot<-ggplot()+
  geom_sf(data=tcv_tow_height_rmspe,aes(col=sqrt(squared_error)/n),cex=2)+
  scale_color_viridis_c(name="RMSPE divided by n")+
  scale_size(name="Number of scallops")+
  theme_bw()

tmb_models2<-tmb_models
tmb_models2[6]<-"spat_both_less_adreport"

set.seed(1823)
test_tcv_weights<-list()
pred_tcv_weights<-list()
tcv_weight_obj<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_weight_Opt<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_weight_rep<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_weight_Report<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_weight_MSPE<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_gam_fit<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())
tcv_gam_pred<-list(off=list(),off_ln=list(),in_ln=list(),spat=list(),spat_off=list(),spat_both=list(),off_estim=list(),off_ln_estim=list())

load("./LW_Work/tcv_weight_output.RData")

for (tow_fold in 63:length(unique(sub_mw$ID))){
for (tow_fold in 56){
  for (model in 6:length(mw_models)){

    test_tcv_weights[[tow_fold]]<-sub_mw[-which(sub_mw$ID==unique(sub_mw$ID)[tow_fold]),]
    pred_tcv_weights[[tow_fold]]<-sub_mw[which(sub_mw$ID==unique(sub_mw$ID)[tow_fold]),]

    tmb_data<-list(weight=test_tcv_weights[[tow_fold]]$wmw,
                   heights=test_tcv_weights[[tow_fold]]$sh/100,
                   tow_id=as.integer(as.factor(test_tcv_weights[[tow_fold]]$ID))-1,
                   weight=test_tcv_weights[[tow_fold]]$wmw,
                   varmat=as.matrix(cbind(1, log(test_tcv_weights[[tow_fold]]$sh/100))),
                   depth=test_tcv_weights[[tow_fold]]$depth,
                   ind_loc=test_tcv_weights[[tow_fold]]$id_loc2-1,
                   locations=st_coordinates(gb_2023_surv))

    if (startsWith(mw_models[model],"Off")){
      parameters<-list(b=3,
                       beta=20,
                       log_phi=-1,
                       log_epsilon=-1,
                       tow_eff=rep(0,length(unique(sub_mw$ID))-1))

      random<-c("tow_eff")

      if (grepl("Estimating",mw_models[model],fixed=T)) maps<-list() else maps<-list(b=factor(NA))
    } else if (startsWith(mw_models[model],"Insh")){
      parameters = list(
        beta = rep(0, nvar),
        beta_depth = rep(0, nvar),
        log_phi = log(0.1))

      maps <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA))
      )

      random<-c()

    } else if (grepl("Both",mw_models[model],fixed=T)) {
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
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA),
        log_nu_b=factor(NA)
      )

      random=c("beta_s","beta_b_s")
    } else {
      parameters = list(
        beta = rep(0, nvar),
        beta_depth = rep(0, nvar),
        beta_s = rep(0,nrow(tmb_data$locations)),
        log_rho = log(10),
        log_sig = log(0.1),
        log_nu = log(1),
        log_phi = log(0.1),
        log_epsilon = log(0.1))#,
      # tow_eff = rep(0,length(unique(sub_mw$ID))))

      maps <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA)
      )

      random<-c("beta_s")
    }

    tcv_weight_obj[[model]][[tow_fold]] = MakeADFun(data=tmb_data,
                     parameters=parameters,
                     map=maps,
                     random=random,
                     DLL=tmb_models2[model],
                     silent = F)

    tcv_weight_Opt[[model]][[tow_fold]]<-optimx::optimr(tcv_weight_obj[[model]][[tow_fold]]$par,tcv_weight_obj[[model]][[tow_fold]]$fn,tcv_weight_obj[[model]][[tow_fold]]$gr,
                         control=list(maxit=100000,maxeval=100000),
                         method="nlminb")
    while (tcv_weight_Opt[[model]][[tow_fold]]$message=="iteration limit reached without convergence (10)"){
      tcv_weight_obj[[model]][[tow_fold]]$par<-tcv_weight_obj[[model]][[tow_fold]]$env$last.par.best[-which(names(tcv_weight_obj[[model]][[tow_fold]]$env$last.par.best)==random)]
      tcv_weight_Opt[[model]][[tow_fold]]<-optimx::optimr(tcv_weight_obj[[model]][[tow_fold]]$par,tcv_weight_obj[[model]][[tow_fold]]$fn,tcv_weight_obj[[model]][[tow_fold]]$gr,
                                                  control=list(maxit=100000,maxeval=100000),
                                                  method="nlminb")
    }

    tcv_weight_rep[[model]][[tow_fold]]<-sdreport(tcv_weight_obj[[model]][[tow_fold]],bias.correct=F)

    tcv_weight_Report[[model]][[tow_fold]]<-tcv_weight_obj[[model]][[tow_fold]]$report()

    if (startsWith(mw_models[model],"Off")) {
      temp_fit<-data.frame(ID=unique(sub_mw$ID)[-tow_fold],CF=(tcv_weight_Report[[model]][[tow_fold]]$beta+tcv_weight_Report[[model]][[tow_fold]]$tow_eff))

      CF.data<-merge(sub_mw[-which(sub_mw$ID==unique(sub_mw$ID)[tow_fold]),][!duplicated(sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit)

      names(CF.data)<-c('ID','lon','lat','year','depth','tow','CF')

      tcv_gam_fit[[model]][[tow_fold]]<-gam(CF~s(lon,lat)+s(depth),data=CF.data)

      Cf.pred<-data.frame(depth=unique(pred_tcv_weights[[tow_fold]]$depth),lon=unique(pred_tcv_weights[[tow_fold]]$lon),lat=unique(pred_tcv_weights[[tow_fold]]$lat))
      tcv_gam_pred[[model]][[tow_fold]]<-predict(tcv_gam_fit[[model]][[tow_fold]],Cf.pred, se=T)$fit
    }

  }
}

county<-1
for (i in 1:5){
  temp_test_tcv_weights<-list()
  temp_pred_tcv_weights<-list()
  temp_tcv_weight_obj<-list(list(),list(),list(),list(),list(),list(),list(),list())
  temp_tcv_weight_Opt<-list(list(),list(),list(),list(),list(),list(),list(),list())
  temp_tcv_weight_rep<-list(list(),list(),list(),list(),list(),list(),list(),list())
  temp_tcv_weight_Report<-list(list(),list(),list(),list(),list(),list(),list(),list())
  temp_tcv_gam_fit<-list(list(),list(),list(),list(),list(),list(),list(),list())
  temp_tcv_gam_pred<-list(list(),list(),list(),list(),list(),list(),list(),list())
  blep<-7
  if (i==5) blep<-13
for (j in 1:blep){
  temp_test_tcv_weights[[j]]<-test_tcv_weights[[j+county]]
  temp_pred_tcv_weights[[j]]<-pred_tcv_weights[[j+county]]
  for (model in 1:8){
    temp_tcv_weight_obj[[model]][[j]]<-tcv_weight_obj[[model]][[j+county]]
    temp_tcv_weight_Opt[[model]][[j]]<-tcv_weight_Opt[[model]][[j+county]]
    temp_tcv_weight_rep[[model]][[j]]<-tcv_weight_rep[[model]][[j+county]]
    temp_tcv_weight_Report[[model]][[j]]<-tcv_weight_Report[[model]][[j+county]]
    if (model %in% c(1,2,7,8)){
      temp_tcv_gam_fit[[model]][[j]]<-tcv_gam_fit[[model]][[j+county]]
      temp_tcv_gam_pred[[model]][[j]]<-tcv_gam_pred[[model]][[j+county]]
    }
  }
}
save(temp_test_tcv_weights,temp_pred_tcv_weights,temp_tcv_weight_obj,
     temp_tcv_weight_Opt,temp_tcv_weight_rep,temp_tcv_weight_Report,
     temp_tcv_gam_fit,temp_tcv_gam_pred,
     # file=paste0("./LW_Work/tcv_weight_output",i+2,".RData"))
     file="./LW_Work/tcv_weight_output2.RData")
  county<-county+10
}

# for (tow_fold in 1:length(unique(sub_mw$ID))){
for (tow_fold in 103:115){
  for (model in 1:length(mw_models)){
    tcv_weight_predictions[[model]][[tow_fold]]<-pred_tcv_weights[[tow_fold]]
    tcv_weight_predictions[[model]][[tow_fold]]$sh<-tcv_weight_predictions[[model]][[tow_fold]]$sh/100
    if (model==1) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-rep(tcv_gam_pred[[model]][[tow_fold]],nrow(tcv_weight_predictions[[model]][[tow_fold]]))*(tcv_weight_predictions[[model]][[tow_fold]]$sh)^3
    if (model==2) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp(log((rep(tcv_gam_pred[[model]][[tow_fold]],nrow(tcv_weight_predictions[[model]][[tow_fold]])))*(tcv_weight_predictions[[model]][[tow_fold]]$sh)^3)-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[2,1])^2)/2)
    if (model==3) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp((summary(tcv_weight_rep[[model]][[tow_fold]])[2,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$depth) + summary(tcv_weight_rep[[model]][[tow_fold]])[1,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$sh))-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[3,1])^2)/2)
    if (model==4) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp((summary(tcv_weight_rep[[model]][[tow_fold]])[4,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$depth) + tcv_weight_Report[[model]][[tow_fold]]$beta_s[pred_tcv_weights[[tow_fold]]$id_loc2] + summary(tcv_weight_rep[[model]][[tow_fold]])[3,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$sh))-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[5,1])^2)/2)
    if (model==5) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp((summary(tcv_weight_rep[[model]][[tow_fold]])[4,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$depth) + (summary(tcv_weight_rep[[model]][[tow_fold]])[3,1]+ tcv_weight_Report[[model]][[tow_fold]]$beta_s[pred_tcv_weights[[tow_fold]]$id_loc2])*log(tcv_weight_predictions[[model]][[tow_fold]]$sh))-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[5,1])^2)/2)
    if (model==6) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp((summary(tcv_weight_rep[[model]][[tow_fold]])[6,1]*log(tcv_weight_predictions[[model]][[tow_fold]]$depth) + tcv_weight_Report[[model]][[tow_fold]]$beta_s[pred_tcv_weights[[tow_fold]]$id_loc2] +(summary(tcv_weight_rep[[model]][[tow_fold]])[5,1]+tcv_weight_Report[[model]][[tow_fold]]$beta_b_s[pred_tcv_weights[[tow_fold]]$id_loc2])*log(tcv_weight_predictions[[model]][[tow_fold]]$sh)))-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[7,1])^2)/2
    if (model==7) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-(rep(tcv_gam_pred[[model]][[tow_fold]],nrow(tcv_weight_predictions[[model]][[tow_fold]])))*(tcv_weight_predictions[[model]][[tow_fold]]$sh)^summary(tcv_weight_rep[[model]][[tow_fold]])[1,1]
    if (model==8) tcv_weight_predictions[[model]][[tow_fold]]$pred_weight<-exp(log((rep(tcv_gam_pred[[model]][[tow_fold]],nrow(tcv_weight_predictions[[model]][[tow_fold]])))*(tcv_weight_predictions[[model]][[tow_fold]]$sh)^summary(tcv_weight_rep[[model]][[tow_fold]])[1,1])-(exp(summary(tcv_weight_rep[[model]][[tow_fold]])[3,1])^2)/2)
    tcv_weight_predictions[[model]][[tow_fold]]$squared_error<-(tcv_weight_predictions[[model]][[tow_fold]]$wmw-tcv_weight_predictions[[model]][[tow_fold]]$pred_weight)^2

  }
}


load("./LW_Work/tcv_weight_preds_MSPE.RData")


for (model in 1:length(tcv_weight_predictions)){
  for (fold in 1:length(unique(sub_mw$ID))){
    if (model==1 & fold ==1){
      all_tcv_weight_pred<-data.frame(tcv_weight_predictions[[model]][[fold]],model=rep(model,nrow(tcv_weight_predictions[[model]][[fold]])),fold=rep(fold,nrow(tcv_weight_predictions[[model]][[fold]])))
    } else {
      all_tcv_weight_pred<-rbind(all_tcv_weight_pred,data.frame(tcv_weight_predictions[[model]][[fold]],model=rep(model,nrow(tcv_weight_predictions[[model]][[fold]])),fold=rep(fold,nrow(tcv_weight_predictions[[model]][[fold]]))))
    }
  }
}

all_tcv_weight_pred$model<-as.factor(mw_models[all_tcv_weight_pred$model])
all_tcv_weight_pred<-all_tcv_weight_pred %>% mutate(model=fct_relevel(model,mw_models))


tcv_weight_rmspe<-aggregate(squared_error~model,data=all_tcv_weight_pred,FUN=mean)
tcv_weight_rmspe$rmspe<-sqrt(tcv_weight_rmspe$squared_error)

tcv_tow_weight_rmspe<-aggregate(squared_error~model+ID,data=all_tcv_weight_pred,FUN=mean)
tcv_tow_weight_rmspe$n<-aggregate(squared_error~model+ID,data=all_tcv_weight_pred,FUN=length)$squared_error

tcv_tow_weight_rmspe$geometry<-rep(NA,nrow(tcv_tow_weight_rmspe))
for (i in 1:nrow(tcv_tow_weight_rmspe)){
  tcv_tow_weight_rmspe$geometry[i]<-gb_2023_surv$geometry[which(gb_2023_surv$ID==tcv_tow_weight_rmspe$ID[i])]
}
tcv_tow_weight_rmspe<-st_as_sf(tcv_tow_weight_rmspe)


tcv_weight_tow_rmspe_plot<-ggplot()+
  # geom_sf(data=tcv_tow_weight_rmspe,aes(col=sqrt(squared_error),size=n))+
  geom_sf(data=subset(tcv_tow_weight_rmspe,model %in% mw_models[1:6]),aes(col=sqrt(squared_error),size=n))+
  scale_color_viridis_c(name="RMSPE")+
  scale_size(name="Number of sampled scallops")+
  facet_wrap(~model)+
  theme_bw()

tcv_weight_tow_rel_rmspe_plot<-ggplot()+
  # geom_sf(data=tcv_tow_weight_rmspe,aes(col=sqrt(squared_error)/n),cex=2)+
  geom_sf(data=subset(tcv_tow_weight_rmspe,model %in% mw_models[1:6]),aes(col=sqrt(squared_error)/n),cex=2)+
  scale_color_viridis_c(name="RMSPE divided by n")+
  # scale_size(name="Number of scallops")+
  facet_wrap(~model)+
  theme_bw()

final_weight_rmspe<-data.frame(model=tcv_weight_rmspe$model,kmean=weights_overall_rmspe$rsq_e,logo=tcv_weight_rmspe$rmspe)

final_weight_rmspe

