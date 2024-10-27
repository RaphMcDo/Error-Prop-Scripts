
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
gb_2023_mw<-all_mw[which(all_mw$year==2023),] %>% st_as_sf(coords=c("lon","lat"))
st_crs(gb_2023_mw)<-4326
gb_2023_mw<-st_transform(gb_2023_mw,crs=32619)
gb_2023_mw$geometry<-gb_2023_mw$geometry/10000

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

non_corrected_heights<-read.csv("./LW_Work/DR2023_13.csv")
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

tows_dist<-st_drop_geometry(gb_2023_surv[,c(2,6,7,8,9)])

tow_start<-st_as_sf(tows_dist,coords=c("slon","slat"))
st_crs(tow_start)<-4326
tow_start<-tow_start %>% st_transform(crs=32619)
tow_end<-st_as_sf(tows_dist,coords=c("elon","elat"))
st_crs(tow_end)<-4326
tow_end<-tow_end %>% st_transform(crs=32619)

actual_dists<-data.frame(tow_id=gb_2023_surv$tow,dist=diag(st_distance(tow_start,tow_end)))

tot_area<-sum(st_area(GB_surv_sf))

midpoints<-seq(0.975,1.675,by=0.05)
mw_models<-c("Offshore Current",
             "Offshore Current LN",
             "Inshore Depth LN GLM",
             "Spatial Depth GLMM",
             "Spatial Offshore LN GLMM",
             "Offshore Estimating b",
             "Offshore LN Estimating",
             "Spatial Both LN")

methods<-c("Midpoint","Individual","Mean Height")

#Models
dyn.load(dynlib("./LW_Work/glmm_offshore_fixed_b"))
dyn.load(dynlib("./LW_Work/glmm_offshore_fixed_b_lognormal"))
dyn.load(dynlib("./LW_Work/depth_glmm_no_space"))
dyn.load(dynlib("./LW_Work/spatial_depth_glmm_direct_cov"))
dyn.load(dynlib("./LW_Work/spatial_depth_glmm_height"))
dyn.load(dynlib("./LW_Work/spatial_depth_glmm_offshore"))
dyn.load(dynlib("./LW_Work/spatial_both"))

#Shell height portion of model!

unique_means<-rep(NA,length(unique(long_heights$tow_id)))

for (i in 1:length(unique_means)){
  unique_means[i]<-which(long_heights$tow_id==unique(long_heights$tow_id)[i])[1]
}

tmb_data5<-list(varmat_sh=as.matrix(cbind(rep(1,nrow(long_heights)))),
                depth_sh=long_heights$depth,
                ind_loc_sh=long_heights$tow_id-1,
                locations=st_coordinates(gb_2023_surv),
                s_heights=long_heights$heights/10,
                a_b_truncation=c(9.5,17),
                unique_means=unique_means,
                n_mod_loc=length(unique(long_heights$tow_id)))

parameters5 = list(
  beta_depth = c(0),
  beta_sh = c(14),
  beta_sh_s = rep(0,nrow(tmb_data5$locations)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_upsilon = log(0.5))

maps5 <- list(
  log_nu=c(factor(NA))#
)

obj5 = MakeADFun(data=tmb_data5,
                 parameters=parameters5,
                 map=maps5,
                 random=c("beta_sh_s"),
                 DLL="spatial_depth_glmm_height",
                 silent = F)

Opt5<-optimx::optimr(obj5$par,obj5$fn,obj5$gr,
                     control=list(maxit=100000,maxeval=100000),
                     method="nlminb")

rep5<-sdreport(obj5,bias.correct=F)

Report5<-obj5$report()

#Prolly need a large number, since sometimes there are 3 models in there
nboot<-1000
n_loc<-nrow(gb_2023_surv)

tot_frame<-list()

for(boot in 1:nboot){

  #Resampling meat weight
  boot_sub_mw<-sub_mw[sample(1:nrow(sub_mw),nrow(sub_mw),replace=T),]

  #OfF

  tmb_data<-list(weight=boot_sub_mw$wmw,
                 heights=boot_sub_mw$sh/100,
                 tow_id=as.integer(as.factor(boot_sub_mw$ID))-1)

  par_list<-list(b=3,
                 beta=20,
                 log_phi=-1,
                 log_epsilon=-1,
                 tow_eff=rep(0,length(unique(boot_sub_mw$ID))))

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

  temp_fit<-data.frame(ID=unique(boot_sub_mw$ID),CF=(Report$beta+Report$tow_eff))

  boot_sub_mw$lon<-st_coordinates(boot_sub_mw)[,1]
  boot_sub_mw$lat<-st_coordinates(boot_sub_mw)[,2]
  CF.data<-merge(boot_sub_mw[!duplicated(boot_sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit)

  names(CF.data)<-c('ID','lon','lat','year','depth','tow','CF')

  CF.fit<-gam(CF~s(lon,lat)+s(depth),data=CF.data)

  #for pred, get CF at each loc, so get mean lon and lat for survey
  gb_2023_surv$lon<-st_coordinates(gb_2023_surv)[,1]
  gb_2023_surv$lat<-st_coordinates(gb_2023_surv)[,2]

  gb_2023_surv$ID<-paste(as.character(gb_2023_surv$cruise),as.character(gb_2023_surv$tow),sep=".")

  pred.loc<-data.frame(depth=gb_2023_surv$depth,lon=gb_2023_surv$lon,lat=gb_2023_surv$lat)

  Cf.pred<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
  Cf.pred$CF<-predict(CF.fit,Cf.pred, se=T)$fit
  Cf.pred.se<-predict(CF.fit,Cf.pred, se=T)$se.fit

  surv.dat<-Cf.pred
  surv.dat$ID<-gb_2023_surv$ID

  # Pull out the ID and condition factor
  tmp.dat<-subset(st_drop_geometry(CF.data),select=c("ID","CF"))
  # Rename CF to CFh
  names(tmp.dat)[2]<-"CFh"
  # merge the two data sets, keeping all x values
  surv.dat<-merge(surv.dat,tmp.dat,all.x=T)
  # Replace any NA's in CFh with the original Condition Factor.
  surv.dat$CFh[is.na(surv.dat$CFh)]<-surv.dat$CF[is.na(surv.dat$CFh)]

  #Off LN

  tmb_data2<-list(weight=boot_sub_mw$wmw,
                  heights=boot_sub_mw$sh/100,
                  tow_id=as.integer(as.factor(boot_sub_mw$ID))-1)

  par_list2<-list(b=3,
                  beta=20,
                  log_phi=-1,
                  log_epsilon=-1,
                  tow_eff=rep(0,length(unique(boot_sub_mw$ID))))

  random2<-c("tow_eff")

  map2<-list(b=factor(NA))

  obj2 = MakeADFun(data=tmb_data2,
                   parameters=par_list2,
                   map=map2,
                   random=random2,
                   DLL="glmm_offshore_fixed_b_lognormal",
                   silent = F)

  Opt2<-optimx::optimr(obj2$par,obj2$fn,obj2$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep2<-sdreport(obj2,bias.correct=F)

  Report2<-obj2$report()

  temp_fit2<-data.frame(ID=unique(boot_sub_mw$ID),CF=(Report2$beta+Report2$tow_eff))

  CF.data2<-merge(boot_sub_mw[!duplicated(boot_sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit2)

  names(CF.data2)<-c('ID','lon','lat','year','depth','tow','CF')

  CF.fit2<-gam(CF~s(lon,lat)+s(depth),data=CF.data2)

  Cf.pred2<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
  Cf.pred2$CF<-predict(CF.fit2,Cf.pred2, se=T)$fit
  Cf.pred2.se<-predict(CF.fit2,Cf.pred2, se=T)$se.fit

  surv.dat2<-Cf.pred2
  surv.dat2$ID<-gb_2023_surv$ID

  # Pull out the ID and condition factor
  tmp.dat2<-subset(st_drop_geometry(CF.data2),select=c("ID","CF"))
  # Rename CF to CFh
  names(tmp.dat2)[2]<-"CFh"
  # merge the two data sets, keeping all x values
  surv.dat2<-merge(surv.dat2,tmp.dat2,all.x=T)
  # Replace any NA's in CFh with the original Condition Factor.
  surv.dat2$CFh[is.na(surv.dat2$CFh)]<-surv.dat2$CF[is.na(surv.dat2$CFh)]


  ##Depth

  tmb_data3<-list(weight=boot_sub_mw$wmw,
                  varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                  depth=boot_sub_mw$depth,
                  ind_loc=boot_sub_mw$id_loc-1)

  nvar = ncol(tmb_data3$varmat)
  nobs = nrow(tmb_data3$varmat)

  parameters3 = list(
    beta = rep(0, nvar),
    beta_depth = rep(0, nvar),
    log_phi = log(0.1))

  maps3 <- list(
    beta = factor(c(NA,2)),
    beta_depth=factor(c(1,NA))
  )

  obj3 = MakeADFun(data=tmb_data3,
                   parameters=parameters3,
                   map=maps3,
                   random=c(),
                   DLL="depth_glmm_no_space",
                   silent = F)

  Opt3<-optimx::optimr(obj3$par,obj3$fn,obj3$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep3<-sdreport(obj3,bias.correct=F)

  Report3<-obj3$report()

  #Spatial

  tmb_data4<-list(weight=boot_sub_mw$wmw,
                  varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                  depth=boot_sub_mw$depth,
                  ind_loc=boot_sub_mw$id_loc2-1,
                  locations=st_coordinates(gb_2023_surv))

  parameters4 = list(
    beta = rep(0, nvar),
    beta_depth = rep(0, nvar),
    beta_s = rep(0,nrow(tmb_data4$locations)),
    log_rho = log(10),
    log_sig = log(0.1),
    log_nu = log(1),
    log_phi = log(0.1))

  maps4 <- list(
    beta = factor(c(NA,2)),
    beta_depth=factor(c(1,NA)),
    log_nu=factor(NA)
  )

  obj4 = MakeADFun(data=tmb_data4,
                   parameters=parameters4,
                   map=maps4,
                   random=c("beta_s"),
                   DLL="spatial_depth_glmm_direct_cov",
                   silent = F)

  Opt4<-optimx::optimr(obj4$par,obj4$fn,obj4$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep4<-sdreport(obj4,bias.correct=F)

  Report4<-obj4$report()

  boot_sub_mw$size_bin<-as.integer(cut(boot_sub_mw$sh,breaks=seq(60,160,by=5),
                                  include.lowest=T))-1
  boot_sub_mw$size_bin_name<-cut(boot_sub_mw$sh,breaks=seq(60,170,by=5),
                            include.lowest=T)

  #Spatial b

  tmb_data6<-list(weight=boot_sub_mw$wmw,
                  varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                  depth=boot_sub_mw$depth,
                  ind_loc=boot_sub_mw$id_loc2-1,
                  locations=st_coordinates(gb_2023_surv))

  parameters6 = list(
    beta = rep(0, nvar),
    beta_depth = rep(0, nvar),
    beta_s = rep(0,nrow(tmb_data6$locations)),
    log_rho = log(10),
    log_sig = log(0.1),
    log_nu = log(1),
    log_phi = log(0.1))

  maps6 <- list(
    beta = factor(c(NA,2)),
    beta_depth=factor(c(1,NA)),
    log_nu=factor(NA)
  )

  obj6 = MakeADFun(data=tmb_data6,
                   parameters=parameters6,
                   map=maps6,
                   random=c("beta_s"),
                   DLL="spatial_depth_glmm_offshore",
                   silent = F)

  Opt6<-optimx::optimr(obj6$par,obj6$fn,obj6$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep6<-sdreport(obj6,bias.correct=F)

  Report6<-obj6$report()

  #Off estimating b

  tmb_data7<-list(weight=boot_sub_mw$wmw,
                  heights=boot_sub_mw$sh/100,
                  tow_id=as.integer(as.factor(boot_sub_mw$ID))-1)

  par_list7<-list(b=3,
                  beta=20,
                  log_phi=-1,
                  log_epsilon=-1,
                  tow_eff=rep(0,length(unique(boot_sub_mw$ID))))

  random7<-c("tow_eff")

  map7<-list()

  obj7 = MakeADFun(data=tmb_data7,
                   parameters=par_list7,
                   map=map7,
                   random=random7,
                   DLL="glmm_offshore_fixed_b",
                   silent = F)

  Opt7<-optimx::optimr(obj7$par,obj7$fn,obj7$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep7<-sdreport(obj7,bias.correct=F)

  Report7<-obj7$report()

  temp_fit7<-data.frame(ID=unique(boot_sub_mw$ID),CF=(Report7$beta+Report7$tow_eff))

  CF.data7<-merge(boot_sub_mw[!duplicated(boot_sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit7)

  names(CF.data7)<-c('ID','lon','lat','year','depth','tow','CF')

  CF.fit7<-gam(CF~s(lon,lat)+s(depth),data=CF.data7)

  Cf.pred7<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
  Cf.pred7$CF<-predict(CF.fit7,Cf.pred7, se=T)$fit
  Cf.pred7.se<-predict(CF.fit7,Cf.pred7, se=T)$se.fit

  surv.dat7<-Cf.pred7
  surv.dat7$ID<-gb_2023_surv$ID

  # Pull out the ID and condition factor
  tmp.dat7<-subset(st_drop_geometry(CF.data7),select=c("ID","CF"))
  # Rename CF to CFh
  names(tmp.dat7)[2]<-"CFh"
  # merge the two data sets, keeping all x values
  surv.dat7<-merge(surv.dat7,tmp.dat7,all.x=T)
  # Replace any NA's in CFh with the original Condition Factor.
  surv.dat7$CFh[is.na(surv.dat7$CFh)]<-surv.dat7$CF[is.na(surv.dat7$CFh)]

  #Off LN estimating b

  tmb_data8<-list(weight=boot_sub_mw$wmw,
                  heights=boot_sub_mw$sh/100,
                  tow_id=as.integer(as.factor(boot_sub_mw$ID))-1)

  par_list8<-list(b=3,
                  beta=20,
                  log_phi=-1,
                  log_epsilon=-1,
                  tow_eff=rep(0,length(unique(boot_sub_mw$ID))))

  random8<-c("tow_eff")

  obj8 = MakeADFun(data=tmb_data8,
                   parameters=par_list8,
                   map=map7,
                   random=random8,
                   DLL="glmm_offshore_fixed_b_lognormal",
                   silent = F)

  Opt8<-optimx::optimr(obj8$par,obj8$fn,obj8$gr,
                       control=list(maxit=100000,maxeval=100000),
                       method="nlminb")

  rep8<-sdreport(obj8,bias.correct=F)

  Report8<-obj8$report()

  temp_fit8<-data.frame(ID=unique(boot_sub_mw$ID),CF=(Report8$beta+Report8$tow_eff))

  CF.data8<-merge(boot_sub_mw[!duplicated(boot_sub_mw$ID),c('ID','lon','lat','year','depth','tow')],temp_fit8)

  names(CF.data8)<-c('ID','lon','lat','year','depth','tow','CF')

  CF.fit8<-gam(CF~s(lon,lat)+s(depth),data=CF.data8)

  Cf.pred8<-data.frame(depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
  Cf.pred8$CF<-predict(CF.fit8,Cf.pred8, se=T)$fit
  Cf.pred8.se<-predict(CF.fit8,Cf.pred8, se=T)$se.fit

  surv.dat8<-Cf.pred8
  surv.dat8$ID<-gb_2023_surv$ID

  # Pull out the ID and condition factor
  tmp.dat8<-subset(st_drop_geometry(CF.data8),select=c("ID","CF"))
  # Rename CF to CFh
  names(tmp.dat8)[2]<-"CFh"
  # merge the two data sets, keeping all x values
  surv.dat8<-merge(surv.dat,tmp.dat8,all.x=T)
  # Replace any NA's in CFh with the original Condition Factor.
  surv.dat8$CFh[is.na(surv.dat8$CFh)]<-surv.dat8$CF[is.na(surv.dat8$CFh)]

  #Spatial Both
  
  tmb_data<-list(weight=boot_sub_mw$wmw,
                 varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                 depth=boot_sub_mw$depth,
                 ind_loc=boot_sub_mw$id_loc2-1,
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

  match_gam_gb<-rep(NA,234)
  for (i in 1:length(match_gam_gb)){
    match_gam_gb[i]<-which(gb_2023_surv$ID==surv.dat$ID[i])
  }

  rematch_gam<-rep(NA,nrow(long_heights))
  for (i in 1:nrow(long_heights)){
    rematch_gam[i]<-which(match_gam_gb==long_heights$tow_id[i])
  }

  pred_midpoints<-data.frame(model=rep(mw_models,each=length(midpoints)*n_loc),
                             midpoints=rep(rep(midpoints,each=nrow(gb_2023_surv)),length(mw_models)),
                             se=rep(midpoints,each=length(mw_models)*nrow(gb_2023_surv)),
                             depth=rep(gb_2023_surv$depth,length(mw_models)*length(midpoints)),
                             loc=rep(1:234,length(mw_models)*length(midpoints)),
                             gam_loc=rep(match_gam_gb,length(mw_models)*length(midpoints)),
                             geometry=rep(gb_2023_surv$geometry,length(mw_models)*length(midpoints)),
                             pred_mid=rep(NA,length(rep(mw_models,each=length(midpoints)*n_loc))))  %>% st_as_sf()

  for (i in 1:nrow(pred_midpoints)){
    if (pred_midpoints$model[i]==mw_models[1]) pred_midpoints$pred_mid[i]<-(surv.dat$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
    if (pred_midpoints$model[i]==mw_models[2]) pred_midpoints$pred_mid[i]<-(surv.dat2$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
    if (pred_midpoints$model[i]==mw_models[3]) pred_midpoints$pred_mid[i]<-exp((summary(rep3)[2,1])*log(pred_midpoints$depth[i])+(summary(rep3)[1,1])*log(pred_midpoints$midpoints[i]))
    if (pred_midpoints$model[i]==mw_models[4]) pred_midpoints$pred_mid[i]<-exp((summary(rep4)[4,1])*log(pred_midpoints$depth[i])+(Report4$beta_s[pred_midpoints$loc])[i]+(summary(rep4)[3,1])*log(pred_midpoints$midpoints[i]))
    if (pred_midpoints$model[i]==mw_models[5]) pred_midpoints$pred_mid[i]<-exp((summary(rep6)[4,1])*log(pred_midpoints$depth[i])+((Report6$beta_s[pred_midpoints$loc])[i]+(summary(rep6)[3,1]))*log(pred_midpoints$midpoints[i]))
    if (pred_midpoints$model[i]==mw_models[6]) pred_midpoints$pred_mid[i]<-(surv.dat7$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^(summary(rep7)[1,1])
    if (pred_midpoints$model[i]==mw_models[7]) pred_midpoints$pred_mid[i]<-(surv.dat8$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^(summary(rep8)[1,1])
    if (pred_midpoints$model[i]==mw_models[8]) pred_midpoints$pred_mid[i]<-exp((summary(rep9)[6,1])*log(pred_midpoints$depth[i])+Report9$beta_s[pred_midpoints$loc[i]]+((Report9$beta_b_s[pred_midpoints$loc])[i]+(summary(rep9)[5,1]))*log(pred_midpoints$midpoints[i]))
  }

  long_heights$cut_heights<-cut(long_heights$heights/100,breaks=c(midpoints-0.025,1.7))
  agg_long_heights<-aggregate(heights~cut_heights+tow_id,data=long_heights,FUN=length)
  agg_long_heights$mid<-midpoints[as.integer(agg_long_heights$cut_heights)]
  join_frame<-data.frame(tow_id=rep(1:234,each=length(midpoints)),mid=rep(midpoints,234))
  agg_long_heights<-left_join(join_frame,agg_long_heights,by=c("tow_id","mid"))
  agg_long_heights$heights[which(is.na(agg_long_heights$heights))]<-0
  colnames(agg_long_heights)[4]<-"total"

  temp_agg_long_heights<-data.frame(agg_long_heights,off=rep(NA,nrow(agg_long_heights)),
                                    off_ln=rep(NA,nrow(agg_long_heights)),insh=rep(NA,nrow(agg_long_heights)),
                                    spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                    off_est=rep(NA,nrow(agg_long_heights)),off_ln_est=rep(NA,nrow(agg_long_heights)),
                                    spat_both=rep(NA,nrow(agg_long_heights)))

  for (model in 1:length(mw_models)){
    for (i in 1:nrow(gb_2023_surv)){
      if (model %in% c(1,2,6,7)){
        temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),model+4]<-temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$gam_loc==i),]$pred_mid
      } else{
        temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),model+4]<-temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$pred_mid
      }
    }
  }

  agg_midpoints_scallops<-data.frame(tow_id=aggregate(off~tow_id,data=temp_agg_long_heights,FUN=sum)$tow_id,
                                     off=aggregate(off~tow_id,data=temp_agg_long_heights,FUN=sum)$off,
                                     off_ln=aggregate(off_ln~tow_id,data=temp_agg_long_heights,FUN=sum)$off_ln,
                                     insh=aggregate(insh~tow_id,data=temp_agg_long_heights,FUN=sum)$insh,
                                     spat=aggregate(spat~tow_id,data=temp_agg_long_heights,FUN=sum)$spat,
                                     spat_off=aggregate(spat_off~tow_id,data=temp_agg_long_heights,FUN=sum)$spat_off,
                                     off_est=aggregate(off_est~tow_id,data=temp_agg_long_heights,FUN=sum)$off_est,
                                     off_ln_est=aggregate(off_ln_est~tow_id,data=temp_agg_long_heights,FUN=sum)$off_ln_est,
                                     spat_both=aggregate(spat_both~tow_id,data=temp_agg_long_heights,FUN=sum)$spat_both)

  long_heights$pred_heights<-Report5$mu_sh/10#+exp(summary(rep5)[5,1])*(dnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1]))-dnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1])))/(pnorm(17,Report5$mu_sh,exp(summary(rep5)[5,1]))-pnorm(9.5,Report5$mu_sh,exp(summary(rep5)[5,1])))
  new_long_heights<-long_heights
  new_long_heights<-new_long_heights[order(new_long_heights$tow_id),]

  for (model in 1:length(mw_models)){
    if (model == 1) new_long_heights$off<-(surv.dat$CFh[rematch_gam])*(new_long_heights$heights/100)^3
    if (model == 2) new_long_heights$off_ln<-(surv.dat2$CFh[rematch_gam])*(new_long_heights$heights/100)^3
    if (model == 3) new_long_heights$insh<-exp((summary(rep3)[2,1])*log(new_long_heights$depth)+(summary(rep3)[1,1])*log(new_long_heights$heights/100))
    if (model == 4) new_long_heights$spat<-exp((summary(rep4)[4,1])*log(new_long_heights$depth)+(Report4$beta_s[new_long_heights$tow_id])+(summary(rep4)[3,1])*log(new_long_heights$heights/100))
    if (model == 5) new_long_heights$spat_off<-exp((summary(rep6)[4,1])*log(new_long_heights$depth)+((Report6$beta_s[new_long_heights$tow_id])+(summary(rep6)[3,1]))*log(new_long_heights$heights/100))
    if (model == 6) new_long_heights$off_est<-(surv.dat7$CFh[rematch_gam])*(new_long_heights$heights/100)^(summary(rep7)[1,1])
    if (model == 7) new_long_heights$off_ln_est<-(surv.dat8$CFh[rematch_gam])*(new_long_heights$heights/100)^(summary(rep8)[1,1])
    if (model == 8) new_long_heights$spat_both<-exp((summary(rep9)[6,1])*log(new_long_heights$depth)+Report9$beta_s[new_long_heights$tow_id]+((Report9$beta_b_s[new_long_heights$tow_id])+(summary(rep9)[5,1]))*log(new_long_heights$heights/100))
  }

  blep<-data.frame(heights=mean(new_long_heights$heights),depth=mean(new_long_heights$depth),
                   tow_id=c(1:234)[-which(1:234 %in% unique(new_long_heights$tow_id))],
                   pred_heights=NA,
                   cut_heights=NA,off=0,off_ln=0,insh=0,spat=0,
                   spat_off=0,off_est=0,off_ln_est=0,spat_both=0)
  new_long_heights<-rbind(new_long_heights,blep)

  agg_indiv_scallops<-data.frame(tow_id=aggregate(off~tow_id,data=new_long_heights,FUN=sum)$tow_id,
                                 off=aggregate(off~tow_id,data=new_long_heights,FUN=sum)$off,
                                 off_ln=aggregate(off_ln~tow_id,data=new_long_heights,FUN=sum)$off_ln,
                                 insh=aggregate(insh~tow_id,data=new_long_heights,FUN=sum)$insh,
                                 spat=aggregate(spat~tow_id,data=new_long_heights,FUN=sum)$spat,
                                 spat_off=aggregate(spat_off~tow_id,data=new_long_heights,FUN=sum)$spat_off,
                                 off_est=aggregate(off_est~tow_id,data=new_long_heights,FUN=sum)$off_est,
                                 off_ln_est=aggregate(off_ln_est~tow_id,data=new_long_heights,FUN=sum)$off_ln_est,
                                 spat_both=aggregate(spat_both~tow_id,data=new_long_heights,FUN=sum)$spat_both)

  small_heights<-long_heights[!duplicated(long_heights$tow_id),c(2,3,5)]
  tow_ns<-aggregate(heights~tow_id,data=long_heights,FUN=length)
  # tow_ns<-rbind(tow_ns,data.frame(tow_id=c(1:234)[-which(1:234 %in% unique(tow_ns$tow_id))],
  #                                 heights=0))
  small_heights<-left_join(tow_ns,small_heights)
  colnames(small_heights)[2]<-"n"
  # small_heights[is.na(small_heights)]<-0

  all_small_heights<-small_heights
  match_gam_gb2<-match_gam_gb[which((1:234 %in% small_heights$tow_id))]

  for (model in 1:length(mw_models)){
    if (model == 1) all_small_heights$off<-((surv.dat$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)*all_small_heights$n
    if (model == 2) all_small_heights$off_ln<-((surv.dat2$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^3)*all_small_heights$n
    if (model == 3) all_small_heights$insh<-(exp((summary(rep3)[2,1])*log(all_small_heights$depth)+(summary(rep3)[1,1])*log(all_small_heights$pred_heights)))*all_small_heights$n
    if (model == 4) all_small_heights$spat<-(exp((summary(rep4)[4,1])*log(all_small_heights$depth)+(Report4$beta_s[all_small_heights$tow_id])+(summary(rep4)[3,1])*log(all_small_heights$pred_heights)))*all_small_heights$n
    if (model == 5) all_small_heights$spat_off<-(exp((summary(rep6)[4,1])*log(all_small_heights$depth)+((Report6$beta_s[all_small_heights$tow_id])+(summary(rep6)[3,1]))*log(all_small_heights$pred_heights)))*all_small_heights$n
    if (model == 6) all_small_heights$off_est<-((surv.dat7$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep7)[1,1]))*all_small_heights$n
    if (model == 7) all_small_heights$off_ln_est<-((surv.dat8$CFh[match_gam_gb2])*(all_small_heights$pred_heights)^(summary(rep8)[1,1]))*all_small_heights$n
    if (model == 8) all_small_heights$spat_both<-exp((summary(rep9)[6,1])*log(all_small_heights$depth)+Report9$beta_s[all_small_heights$tow_id]+((Report9$beta_b_s[all_small_heights$tow_id])+(summary(rep9)[5,1]))*log(all_small_heights$pred_heights))*all_small_heights$n
  }

  all_sums<-data.frame(tow=rep(agg_midpoints_scallops$tow_id,length(mw_models)),
                       val=c(agg_midpoints_scallops$off,agg_midpoints_scallops$off_ln,
                             agg_midpoints_scallops$insh,agg_midpoints_scallops$spat,
                             agg_midpoints_scallops$spat_off,agg_midpoints_scallops$off_est,
                             agg_midpoints_scallops$off_ln_est,agg_midpoints_scallops$spat_both),
                       model=rep(mw_models,each=nrow(agg_midpoints_scallops)),
                       step2=rep("Midpoint",nrow(agg_midpoints_scallops)*length(mw_models)),
                       geometry=rep(gb_2023_surv$geometry,length(mw_models)))
  temp_frame<-left_join(data.frame(tow_id=1:234),agg_indiv_scallops,by="tow_id")
  temp_frame[is.na(temp_frame)]<-0
  all_sums<-rbind(all_sums,
                  data.frame(tow=rep(temp_frame$tow_id,length(mw_models)),
                             val=c(temp_frame$off,temp_frame$off_ln,
                                   temp_frame$insh,temp_frame$spat,
                                   temp_frame$spat_off,temp_frame$off_est,
                                   temp_frame$off_ln_est,temp_frame$spat_both),
                             model=rep(mw_models,each=nrow(temp_frame)),
                             step2=rep("Individual",nrow(temp_frame)*length(mw_models)),
                             geometry=rep(gb_2023_surv$geometry,length(mw_models))))
  temp_frame<-left_join(data.frame(tow_id=1:234),all_small_heights[,c(1,5:12)],by="tow_id")
  temp_frame[is.na(temp_frame)]<-0
  all_sums<-rbind(all_sums,
                  data.frame(tow=rep(temp_frame$tow_id,length(mw_models)),
                             val=c(temp_frame$off,temp_frame$off_ln,
                                   temp_frame$insh,temp_frame$spat,
                                   temp_frame$spat_off,temp_frame$off_est,
                                   temp_frame$off_ln_est,temp_frame$spat_both),
                             model=rep(mw_models,each=nrow(temp_frame)),
                             step2=rep("Mean Height",nrow(temp_frame)*length(mw_models)),
                             geometry=rep(gb_2023_surv$geometry,length(mw_models)))) %>%
    st_as_sf()
  all_sums$val<-all_sums$val/1000

  calc_tot<-list(c(rep(NA,8),"Midpoint"),c(rep(NA,8),"Individual"),c(rep(NA,8),"Mean Height"))

  transf_val<-((1/(as.numeric(actual_dists$dist)*2.4384/10^6))*tot_area)/1000

  for (method in 1:3){
    for (model in 1:length(mw_models)){
      calc_tot[[method]][model]<-mean((st_drop_geometry(all_sums)[which(all_sums$model==mw_models[model] & all_sums$step2==methods[method]),]$val*transf_val))
    }
  }

  tot_frame[[boot]]<-data.frame(model=rep(mw_models,3),
                        method=rep(methods,each=length(mw_models)),
                        est=c(as.numeric(calc_tot[[1]][-9]),as.numeric(calc_tot[[2]][-9]),as.numeric(calc_tot[[3]][-9])))

}

final_tot_frame<-tot_frame[[1]]
for (i in 2:nboot){
  final_tot_frame<-rbind(final_tot_frame,tot_frame[[i]])
}

mw_models2<-c("Off","Off LN", "Depth","Spatial",
             "Spatial b", mw_models[6:7],"Spatial Both")

for (i in 1:length(mw_models)){
  final_tot_frame$model[which(final_tot_frame$model==mw_models[i])]<-rep(mw_models2[i],length(final_tot_frame$model[which(final_tot_frame$model==mw_models[i])]))
}

library(forcats)
final_tot_frame<-final_tot_frame %>% mutate(model=fct_relevel(model,mw_models2))

big_hists_all<-ggplot(data=subset(final_tot_frame,model %in% mw_models2[c(1:5,8)]))+
  geom_histogram(aes(x=est),alpha=0.2,col="black")+
  theme_bw()+
  ylab("Count")+xlab("Estimated Index")+
  facet_grid(rows=vars(model),cols=vars(method))

up_quant<-function(x){
  quantile(x,0.975)
}

low_quant<-function(x){
  quantile(x,0.025)
}

final_ests<-aggregate(est~model+method,data=final_tot_frame,FUN=mean)
final_up_ci<-aggregate(est~model+method,data=final_tot_frame,FUN=up_quant)
final_low_ci<-aggregate(est~model+method,data=final_tot_frame,FUN=low_quant)

final_est_ci<-data.frame(final_ests,up=final_up_ci$est,low=final_low_ci$est)

final_est_plot<-ggplot(data=subset(final_est_ci,model %in% mw_models2[c(1:5,8)]))+
  geom_point(aes(x=model,y=est,col=model,shape=model))+
  geom_errorbar(aes(x=model,col=model,ymin=low,ymax=up))+
  coord_cartesian(ylim=c(6600,7900))+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),axis.title.x = element_blank())+
  scale_color_viridis_d(name="Model")+
  scale_shape(name="Model")+
  ylab("Bootstrapped Index")+
  facet_wrap(~method)

#Midpoint - Individual
subset(final_est_ci,method==methods[1])$est-subset(final_est_ci,method==methods[2])$est
#Midpoint - Mean Height
subset(final_est_ci,method==methods[1])$est-subset(final_est_ci,method==methods[3])$est
#Individual - Mean Height
subset(final_est_ci,method==methods[2])$est-subset(final_est_ci,method==methods[3])$est

final_est_ci$ci_size<-final_est_ci$up-final_est_ci$low

final_est_ci$scaled_ci_size<-final_est_ci$ci_size/final_est_ci$est









