

#########################################################
############    MW-SH Modelling #########################
#########################################################

#Will focus on just spatial model, so will only pick 1 year to start
#Focusing on GB 

library(TMB)
library(sf)
library(raster)
library(mgcv)
library(ggplot2)
library(dplyr)

up_quant<-function(x){
  quantile(x,0.975)
}

low_quant<-function(x){
  quantile(x,0.025)
}

tot_area<-4252.073

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

#Models
# compile("./spatial_depth_glmm_offshore.cpp")
dyn.load(dynlib("./glmm_offshore_fixed_b"))
dyn.load(dynlib("./glmm_offshore_fixed_b_lognormal"))
dyn.load(dynlib("./depth_glmm_no_space"))
dyn.load(dynlib("./spatial_depth_glmm_direct_cov"))
dyn.load(dynlib("./spatial_depth_glmm_offshore"))

midpoints<-seq(0.975,1.675,by=0.05)
mw_models<-c("Offshore Current",
             "Offshore Current LN",
             "Inshore Depth LN GLM",
             "Spatial Depth GLMM",
             "Spatial Offshore LN GLMM",
             "Spatial Both")
#Fitting the first model

tmb_data_start<-list(weight=c(sub_mw$wmw),
                     varmat=as.matrix(cbind(1, c(log(sub_mw$sh/100)))),
                     depth=c(sub_mw$depth),
                     ind_loc=c(sub_mw$id_loc2-1),
                     locations=st_coordinates(gb_2023_surv))

parameters_start = list(
  beta = rep(0, 2),
  beta_depth = rep(0, 2),
  beta_s = rep(0,nrow(gb_2023_surv)),
  log_rho = log(10),
  log_sig = log(0.1),
  log_nu = log(1),
  log_phi = log(0.1))

maps_start <- list(
  # beta = factor(rep(NA, nvar)),
  # beta = factor(c(2,NA)),
  beta = factor(c(NA,2)),
  # beta = factor(c(3,2)),
  beta_depth=factor(c(1,NA)),
  log_nu=factor(NA)
)

#Spatial depth LMM
obj_start = MakeADFun(data=tmb_data_start,
                      parameters=parameters_start,
                      map=maps_start,
                      random=c("beta_s"),
                      DLL="spatial_depth_glmm_direct_cov",
                      silent = F)

Opt_start<-optimx::optimr(obj_start$par,obj_start$fn,obj_start$gr,
                          control=list(maxit=100000,maxeval=100000),
                          method="nlminb")

rep_start<-sdreport(obj_start,bias.correct=F)

Report_start<-obj_start$report()

#Did checkConsistency including simulating the random effects
# check<-checkConsistency(obj_start)


nboot<-1000
nvar<-2
n_loc<-nrow(gb_2023_surv)

n_sim<-10
subs<-200/n_sim
sim_index<-rep(NA,n_sim)
calc_tot_nonboot<-list()
calc_high_nonboot<-list()
calc_low_nonboot<-list()
final_est_ci<-list()
set.seed(190)
for (sub in 5:subs){
  for (sim in 1:n_sim){
    sim_data<-list(weight=rep(1,nrow(long_heights)),
                   varmat=as.matrix(cbind(1, c(log(long_heights$heights/100)))),
                   depth=long_heights$depth,
                   ind_loc=long_heights$tow_id-1,
                   locations=st_coordinates(gb_2023_surv))
    sim_par<-list(
      beta=c(0,obj_start$env$last.par.best[3]),
      beta_depth=c(obj_start$env$last.par.best[238],0),
      beta_s=c(obj_start$env$last.par.best[4:237]),
      log_rho=obj_start$env$last.par.best[1],
      log_sig=obj_start$env$last.par.best[2],
      log_nu=log(1),
      log_phi=obj_start$env$last.par.best[239]
    )
    
    sim_maps <- list(
      # beta = factor(rep(NA, nvar)),
      # beta = factor(c(2,NA)),
      beta = factor(c(NA,2)),
      # beta = factor(c(3,2)),
      beta_depth=factor(c(1,NA)),
      log_nu=factor(NA),
      beta_s=factor(rep(NA,nrow(gb_2023_surv)))
    )
    sim_obj<-MakeADFun(data=sim_data,
                       parameters=sim_par,
                       map=sim_maps,
                       random=c("beta_s"),
                       DLL="spatial_depth_glmm_direct_cov",
                       silent = F)
    
    full_sim_data<-sim_obj$simulate()
    
    agg_sims<-data.frame(weights=full_sim_data$weight,tow_id=sim_data$ind_loc+1)
    agg_sims<-aggregate(weights~tow_id,data=agg_sims,FUN=sum)
    agg_sims<-left_join(data.frame(tow_id=1:nrow(gb_2023_surv)),agg_sims)
    agg_sims$weights[is.na(agg_sims$weights)]<-0
    
    #Simulated index
    sim_index[sim]<-mean(agg_sims$weights/(actual_dists$dist*2.4384/10^6))*tot_area/1000000
    
    #Simulating the subset we interested in:
    sim_data2<-tmb_data_start
    
    sim_obj2<-MakeADFun(data=sim_data2,
                        parameters=sim_par,
                        map=sim_maps,
                        random=c("beta_s"),
                        DLL="spatial_depth_glmm_direct_cov",
                        silent = F)
    sim_sub_mw_data<-sim_obj2$simulate()
    
    new_sub_mw<-sub_mw
    new_sub_mw$wmw<-sim_sub_mw_data$weight
    
    #First, directly obtaining the CIs
    parameters3 = list(
      beta = rep(0, nvar),
      beta_depth = rep(0, nvar),
      log_phi = log(0.1))
    
    maps3 <- list(
      # beta = factor(rep(NA, nvar)),
      # beta = factor(c(2,NA)),
      beta = factor(c(NA,2)),
      # beta = factor(c(3,2)),
      beta_depth=factor(c(1,NA))
    )
    
    parameters4 = list(
      beta = rep(0, nvar),
      beta_depth = rep(0, nvar),
      beta_s = rep(0,nrow(gb_2023_surv)),
      log_rho = log(10),
      log_sig = log(0.1),
      log_nu = log(1),
      log_phi = log(0.1))
    
    maps4 <- list(
      # beta = factor(rep(NA, nvar)),
      # beta = factor(c(2,NA)),
      beta = factor(c(NA,2)),
      # beta = factor(c(3,2)),
      beta_depth=factor(c(1,NA)),
      log_nu=factor(NA)
    )
    
    parameters6 = list(
      beta = rep(0, nvar),
      beta_depth = rep(0, nvar),
      beta_s = rep(0,nrow(gb_2023_surv)),
      log_rho = log(10),
      log_sig = log(0.1),
      log_nu = log(1),
      log_phi = log(0.1))
    
    maps6 <- list(
      # beta = factor(rep(NA, nvar)),
      # beta = factor(c(2,NA)),
      beta = factor(c(NA,2)),
      # beta = factor(c(3,2)),
      beta_depth=factor(c(1,NA)),
      log_nu=factor(NA)
    )
    
    tmb_data3<-list(weight=c(new_sub_mw$wmw,rep(NA,length(midpoints)*nrow(gb_2023_surv))),
                    varmat=as.matrix(cbind(1, c(log(new_sub_mw$sh/100),log(rep(midpoints,nrow(gb_2023_surv)))))),
                    depth=c(new_sub_mw$depth,rep(gb_2023_surv$depth,each=length(midpoints))),
                    # ind_loc=mesh$idx$loc-1,
                    ind_loc=c(new_sub_mw$id_loc2-1,rep(0:(nrow(gb_2023_surv)-1),each=length(midpoints))))
    
    tmb_data4<-list(weight=c(new_sub_mw$wmw,rep(NA,length(midpoints)*nrow(gb_2023_surv))),
                    varmat=as.matrix(cbind(1, c(log(new_sub_mw$sh/100),log(rep(midpoints,nrow(gb_2023_surv)))))),
                    depth=c(new_sub_mw$depth,rep(gb_2023_surv$depth,each=length(midpoints))),
                    ind_loc=c(new_sub_mw$id_loc2-1,rep(0:(nrow(gb_2023_surv)-1),each=length(midpoints))),
                    locations=st_coordinates(gb_2023_surv))
    
    tmb_data6<-list(weight=c(new_sub_mw$wmw,rep(NA,length(midpoints)*nrow(gb_2023_surv))),
                    varmat=as.matrix(cbind(1, c(log(new_sub_mw$sh/100),log(rep(midpoints,nrow(gb_2023_surv)))))),
                    depth=c(new_sub_mw$depth,rep(gb_2023_surv$depth,each=length(midpoints))),
                    ind_loc=c(new_sub_mw$id_loc2-1,rep(0:(nrow(gb_2023_surv)-1),each=length(midpoints))),
                    locations=st_coordinates(gb_2023_surv))
    
    nobs = nrow(tmb_data3$varmat)
    
    #Depth LM
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
    
    #Spatial depth LMM
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
    
    #Spatial Depth lmm (offshore formulation)
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
    
    tmb_data9<-list(weight=sub_mw$wmw,
                    varmat=as.matrix(cbind(1, log(sub_mw$sh/100))),
                    depth=sub_mw$depth,
                    ind_loc=sub_mw$id_loc2-1,
                    locations=st_coordinates(gb_2023_surv))
    
    parameters9 = list(
      beta = rep(0, nvar),
      beta_depth = rep(0, nvar),
      beta_s = rep(0,nrow(tmb_data9$locations)),
      beta_b_s = rep(0,nrow(tmb_data9$locations)),
      log_rho = log(10),
      log_sig = log(0.1),
      log_nu = log(1),
      log_rho_b = log(10),
      log_sig_b = log(0.1),
      log_nu_b = log(1),
      log_phi = log(0.1))
    
    maps9 <- list(
      beta = factor(c(NA,2)),
      beta_depth=factor(c(1,NA)),
      log_nu=factor(NA),
      log_nu_b=factor(NA)
    )
    
    non_r<-names(parameters)[-which(names(parameters) %in% c("beta_s","beta_b_s"))]
    
    obj9 = MakeADFun(data=tmb_data9, 
                     parameters=parameters9,
                     map=maps9,
                     random=c("beta_s","beta_b_s"),
                     DLL="spatial_both",
                     silent = F)
    
    Opt9<-optimx::optimr(obj9$par,obj9$fn,obj9$gr,
                         control=list(maxit=100000,maxeval=100000),
                         method="nlminb")
    
    rep9<-sdreport(obj9,bias.correct=F)
    
    Report9<-obj9$report()
    
    pred_midpoints<-data.frame(preds=c(exp(Report3$mu[-c(1:nrow(new_sub_mw))]),exp(Report4$mu[-c(1:nrow(new_sub_mw))]),
                                       exp(Report6$mu[-c(1:nrow(new_sub_mw))]),exp(Report9$mu[-c(1:nrow(new_sub_mw))])),
                               high=c(exp(Report3$mu[-c(1:nrow(new_sub_mw))]+2*rep3$sd[which(names(rep3$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                      exp(Report4$mu[-c(1:nrow(new_sub_mw))]+2*rep4$sd[which(names(rep4$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                      exp(Report6$mu[-c(1:nrow(new_sub_mw))]+2*rep6$sd[which(names(rep6$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                      exp(Report9$mu[-c(1:nrow(new_sub_mw))]+2*rep9$sd[which(names(rep9$value)=="mu")][-c(1:nrow(new_sub_mw))])),
                               low=c(exp(Report3$mu[-c(1:nrow(new_sub_mw))]-2*rep3$sd[which(names(rep3$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                     exp(Report4$mu[-c(1:nrow(new_sub_mw))]-2*rep4$sd[which(names(rep4$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                     exp(Report6$mu[-c(1:nrow(new_sub_mw))]-2*rep6$sd[which(names(rep6$value)=="mu")][-c(1:nrow(new_sub_mw))]),
                                     exp(Report9$mu[-c(1:nrow(new_sub_mw))]-2*rep9$sd[which(names(rep9$value)=="mu")][-c(1:nrow(new_sub_mw))])),
                               model=c(rep(mw_models[3],length(Report3$mu[-c(1:nrow(new_sub_mw))])),
                                       rep(mw_models[4],length(Report4$mu[-c(1:nrow(new_sub_mw))])),
                                       rep(mw_models[5],length(Report6$mu[-c(1:nrow(new_sub_mw))])),
                                       rep(mw_models[6],length(Report9$mu[-c(1:nrow(new_sub_mw))]))),
                               loc=c(tmb_data3$ind_loc[-c(1:nrow(new_sub_mw))]+1,tmb_data4$ind_loc[-c(1:nrow(new_sub_mw))]+1,
                                     tmb_data6$ind_loc[-c(1:nrow(new_sub_mw))]+1,tmb_data9$ind_loc[-c(1:nrow(new_sub_mw))]+1))
    
    new_long_heights<-long_heights
    new_long_heights$cut_heights<-cut(new_long_heights$heights/100,breaks=c(midpoints-0.025,1.7))
    agg_long_heights<-aggregate(heights~cut_heights+tow_id,data=new_long_heights,FUN=length)
    agg_long_heights$mid<-midpoints[as.integer(agg_long_heights$cut_heights)]
    join_frame<-data.frame(tow_id=rep(1:234,each=length(midpoints)),mid=rep(midpoints,234))
    agg_long_heights<-left_join(join_frame,agg_long_heights,by=c("tow_id","mid"))
    agg_long_heights$heights[which(is.na(agg_long_heights$heights))]<-0
    colnames(agg_long_heights)[4]<-"total"
    
    temp_agg_long_heights<-data.frame(agg_long_heights,insh=rep(NA,nrow(agg_long_heights)),
                                      spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                      spat_both=rep(NA,nrow(agg_long_heights)))
    
    temp_high_long_heights<-data.frame(agg_long_heights,insh=rep(NA,nrow(agg_long_heights)),
                                       spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                       spat_both=rep(NA,nrow(agg_long_heights)))
    
    temp_low_long_heights<-data.frame(agg_long_heights,insh=rep(NA,nrow(agg_long_heights)),
                                      spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                      spat_both=rep(NA,nrow(agg_long_heights)))
    
    for (model in 3:6){
      for (i in 1:nrow(gb_2023_surv)){
        temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),model+2]<-temp_agg_long_heights[which(temp_agg_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$preds
        temp_high_long_heights[which(temp_high_long_heights$tow_id==i),model+2]<-temp_high_long_heights[which(temp_high_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$high
        temp_low_long_heights[which(temp_low_long_heights$tow_id==i),model+2]<-temp_low_long_heights[which(temp_low_long_heights$tow_id==i),]$total*pred_midpoints[which(pred_midpoints$model==mw_models[model] & pred_midpoints$loc==i),]$low
      }
    }
    
    agg_midpoints_scallops<-data.frame(tow_id=aggregate(insh~tow_id,data=temp_agg_long_heights,FUN=sum)$tow_id,
                                       insh=aggregate(insh~tow_id,data=temp_agg_long_heights,FUN=sum)$insh,
                                       spat=aggregate(spat~tow_id,data=temp_agg_long_heights,FUN=sum)$spat,
                                       spat_off=aggregate(spat_off~tow_id,data=temp_agg_long_heights,FUN=sum)$spat_off,
                                       spat_both=aggregate(spat_both~tow_id,data=temp_agg_long_heights,FUN=sum)$spat_both)
    
    high_midpoints_scallops<-data.frame(tow_id=aggregate(insh~tow_id,data=temp_high_long_heights,FUN=sum)$tow_id,
                                        insh=aggregate(insh~tow_id,data=temp_high_long_heights,FUN=sum)$insh,
                                        spat=aggregate(spat~tow_id,data=temp_high_long_heights,FUN=sum)$spat,
                                        spat_off=aggregate(spat_off~tow_id,data=temp_high_long_heights,FUN=sum)$spat_off,
                                        spat_both=aggregate(spat_both~tow_id,data=temp_high_long_heights,FUN=sum)$spat_both)
    
    low_midpoints_scallops<-data.frame(tow_id=aggregate(insh~tow_id,data=temp_low_long_heights,FUN=sum)$tow_id,
                                       insh=aggregate(insh~tow_id,data=temp_low_long_heights,FUN=sum)$insh,
                                       spat=aggregate(spat~tow_id,data=temp_low_long_heights,FUN=sum)$spat,
                                       spat_off=aggregate(spat_off~tow_id,data=temp_low_long_heights,FUN=sum)$spat_off,
                                       spat_both=aggregate(spat_both~tow_id,data=temp_low_long_heights,FUN=sum)$spat_both)
    
    
    all_sums<-data.frame(tow=rep(agg_midpoints_scallops$tow_id,length(mw_models[3:6])),
                         val=c(agg_midpoints_scallops$insh,agg_midpoints_scallops$spat,
                               agg_midpoints_scallops$spat_off,agg_midpoints_scallops$spat_both),
                         model=rep(mw_models[3:6],each=nrow(agg_midpoints_scallops)),
                         geometry=rep(gb_2023_surv$geometry,length(mw_models[3:6]))) %>% st_as_sf()
    all_sums$val<-all_sums$val/1000
    
    all_high_sums<-data.frame(tow=rep(high_midpoints_scallops$tow_id,length(mw_models[3:6])),
                              val=c(high_midpoints_scallops$insh,high_midpoints_scallops$spat,
                                    high_midpoints_scallops$spat_off,high_midpoints_scallops$spat_both),
                              model=rep(mw_models[3:6],each=nrow(high_midpoints_scallops)),
                              geometry=rep(gb_2023_surv$geometry,length(mw_models[3:6]))) %>% st_as_sf()
    all_high_sums$val<-all_high_sums$val/1000
    
    all_low_sums<-data.frame(tow=rep(low_midpoints_scallops$tow_id,length(mw_models[3:6])),
                             val=c(low_midpoints_scallops$insh,low_midpoints_scallops$spat,
                                   low_midpoints_scallops$spat_off,low_midpoints_scallops$spat_both),
                             model=rep(mw_models[3:6],each=nrow(low_midpoints_scallops)),
                             geometry=rep(gb_2023_surv$geometry,length(mw_models[3:6]))) %>% st_as_sf()
    all_low_sums$val<-all_low_sums$val/1000
    
    calc_tot_nonboot[[sim]]<-c(NA,NA,NA,NA)
    calc_high_nonboot[[sim]]<-c(NA,NA,NA,NA)
    calc_low_nonboot[[sim]]<-c(NA,NA,NA,NA)
    
    for (model in 1:4){
      calc_tot_nonboot[[sim]][model]<-mean((st_drop_geometry(all_sums)[which(all_sums$model==mw_models[model+2]),]$val/(actual_dists$dist*2.4384/10^6))*tot_area)/1000
      calc_high_nonboot[[sim]][model]<-mean((st_drop_geometry(all_high_sums)[which(all_high_sums$model==mw_models[model+2]),]$val/(actual_dists$dist*2.4384/10^6))*tot_area)/1000
      calc_low_nonboot[[sim]][model]<-mean((st_drop_geometry(all_low_sums)[which(all_low_sums$model==mw_models[model+2]),]$val/(actual_dists$dist*2.4384/10^6))*tot_area)/1000
    }
    
    tot_frame<-list()
    for(boot in 1:nboot){
      
      #Resampling meat weight
      boot_sub_mw<-new_sub_mw[sample(1:nrow(new_sub_mw),nrow(new_sub_mw),replace=T),]
      
      #Offshore
      
      #doing it in TMB check if get same
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
      
      #Lognormal offshore
      
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
      
      
      ##Depth glm
      
      tmb_data3<-list(weight=boot_sub_mw$wmw,
                      varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                      depth=boot_sub_mw$depth,
                      # ind_loc=mesh$idx$loc-1,
                      ind_loc=boot_sub_mw$id_loc-1)
      
      nvar = ncol(tmb_data3$varmat)
      nobs = nrow(tmb_data3$varmat)
      
      parameters3 = list(
        beta = Report3$beta,
        beta_depth = Report3$beta_depth,
        log_phi = summary(rep3)[3,1])
      
      maps3 <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA))
      )
      
      boot_obj3 = MakeADFun(data=tmb_data3,
                            parameters=parameters3,
                            map=maps3,
                            random=c(),
                            DLL="depth_glmm_no_space",
                            silent = F)
      
      boot_Opt3<-optimx::optimr(boot_obj3$par,boot_obj3$fn,boot_obj3$gr,
                                control=list(maxit=100000,maxeval=100000),
                                method="nlminb")
      
      # rep3<-sdreport(obj3,bias.correct=F)
      
      boot_Report3<-boot_obj3$report()
      
      #Spatial depth GLMM
      
      tmb_data4<-list(weight=boot_sub_mw$wmw,
                      varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                      depth=boot_sub_mw$depth,
                      ind_loc=boot_sub_mw$id_loc2-1,
                      locations=st_coordinates(gb_2023_surv))
      
      parameters4 = list(
        beta = Report4$beta,
        beta_depth = Report4$beta_depth,
        beta_s = Report4$beta_s,
        log_rho = summary(rep4)[1,1],
        log_sig = summary(rep4)[2,1],
        log_nu = log(1),
        log_phi = summary(rep4)[5,1])
      
      maps4 <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA)
      )
      
      boot_obj4 = MakeADFun(data=tmb_data4,
                            parameters=parameters4,
                            map=maps4,
                            random=c("beta_s"),
                            DLL="spatial_depth_glmm_direct_cov",
                            silent = F)
      
      boot_Opt4<-optimx::optimr(boot_obj4$par,boot_obj4$fn,boot_obj4$gr,
                                control=list(maxit=100000,maxeval=100000),
                                method="nlminb")
      
      # rep4<-sdreport(obj4,bias.correct=F)
      
      boot_Report4<-boot_obj4$report()
      
      boot_sub_mw$size_bin<-as.integer(cut(boot_sub_mw$sh,breaks=seq(60,160,by=5),
                                           include.lowest=T))-1
      boot_sub_mw$size_bin_name<-cut(boot_sub_mw$sh,breaks=seq(60,170,by=5),
                                     include.lowest=T)
      
      #Spatial Depth glmm (offshore formulation)
      
      tmb_data6<-list(weight=boot_sub_mw$wmw,
                      varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                      depth=boot_sub_mw$depth,
                      ind_loc=boot_sub_mw$id_loc2-1,
                      locations=st_coordinates(gb_2023_surv))
      
      parameters6 = list(
        beta = Report6$beta,
        beta_depth = Report6$beta_depth,
        beta_s = Report6$beta_s,
        log_rho = summary(rep6)[1,1],
        log_sig = summary(rep6)[2,1],
        log_nu = log(1),
        log_phi = summary(rep6)[5,1])
      
      maps6 <- list(
        # beta = factor(rep(NA, nvar)),
        # beta = factor(c(2,NA)),
        beta = factor(c(NA,2)),
        # beta = factor(c(3,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA)
      )
      
      boot_obj6 = MakeADFun(data=tmb_data6,
                            parameters=parameters6,
                            map=maps6,
                            random=c("beta_s"),
                            DLL="spatial_depth_glmm_offshore",
                            silent = F)
      
      boot_Opt6<-optimx::optimr(boot_obj6$par,boot_obj6$fn,boot_obj6$gr,
                                control=list(maxit=100000,maxeval=100000),
                                method="nlminb")
      
      # rep6<-sdreport(obj6,bias.correct=F)
      
      boot_Report6<-boot_obj6$report()
      
      tmb_data9<-list(weight=boot_sub_mw$wmw,
                      varmat=as.matrix(cbind(1, log(boot_sub_mw$sh/100))),
                      depth=boot_sub_mw$depth,
                      ind_loc=boot_sub_mw$id_loc2-1,
                      locations=st_coordinates(gb_2023_surv))
      
      parameters9 = list(
        beta = rep(0, nvar),
        beta_depth = rep(0, nvar),
        beta_s = rep(0,nrow(tmb_data9$locations)),
        beta_b_s = rep(0,nrow(tmb_data9$locations)),
        log_rho = log(10),
        log_sig = log(0.1),
        log_nu = log(1),
        log_rho_b = log(10),
        log_sig_b = log(0.1),
        log_nu_b = log(1),
        log_phi = log(0.1))
      
      maps9 <- list(
        beta = factor(c(NA,2)),
        beta_depth=factor(c(1,NA)),
        log_nu=factor(NA),
        log_nu_b=factor(NA)
      )
      
      non_r<-names(parameters)[-which(names(parameters) %in% c("beta_s","beta_b_s"))]
      
      boot_obj9 = MakeADFun(data=tmb_data9, 
                            parameters=parameters9,
                            map=maps9,
                            random=c("beta_s","beta_b_s"),
                            DLL="spatial_both",
                            silent = F)
      
      boot_Opt9<-optimx::optimr(boot_obj9$par,boot_obj9$fn,boot_obj9$gr,
                                control=list(maxit=100000,maxeval=100000),
                                method="nlminb")
      
      boot_rep9<-sdreport(boot_obj9,bias.correct=F)
      
      boot_Report9<-boot_obj9$report()
      
      match_gam_gb<-rep(NA,234)
      for (i in 1:length(match_gam_gb)){
        match_gam_gb[i]<-which(gb_2023_surv$ID==surv.dat$ID[i])
      }
      
      rematch_gam<-rep(NA,nrow(long_heights))
      for (i in 1:nrow(long_heights)){
        rematch_gam[i]<-which(match_gam_gb==long_heights$tow_id[i])
      }
      
      pred_midpoints<-data.frame(model=rep(mw_models,length(midpoints)*n_loc),
                                 midpoints=rep(midpoints,each=length(mw_models)*nrow(gb_2023_surv)),
                                 depth=rep(gb_2023_surv$depth,length(mw_models)*length(midpoints)),
                                 loc=rep(1:234,length(mw_models)*length(midpoints)),
                                 gam_loc=rep(match_gam_gb,length(mw_models)*length(midpoints)),
                                 geometry=rep(gb_2023_surv$geometry,length(mw_models)*length(midpoints)),
                                 pred_mid=rep(NA,length(rep(mw_models,each=length(midpoints)*n_loc)))) %>% st_as_sf()
      
      for (i in 1:nrow(pred_midpoints)){
        if (pred_midpoints$model[i]==mw_models[1]) pred_midpoints$pred_mid[i]<-(surv.dat$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
        if (pred_midpoints$model[i]==mw_models[2]) pred_midpoints$pred_mid[i]<-(surv.dat2$CFh[pred_midpoints$gam_loc])[i]*pred_midpoints$midpoints[i]^3
        if (pred_midpoints$model[i]==mw_models[3]) pred_midpoints$pred_mid[i]<-exp((boot_Report3$beta_depth[1])*log(pred_midpoints$depth[i])+(boot_Report3$beta[2])*log(pred_midpoints$midpoints[i]))
        if (pred_midpoints$model[i]==mw_models[4]) pred_midpoints$pred_mid[i]<-exp((boot_Report4$beta_depth[1])*log(pred_midpoints$depth[i])+(boot_Report4$beta_s[pred_midpoints$loc])[i]+(boot_Report4$beta[2])*log(pred_midpoints$midpoints[i]))
        if (pred_midpoints$model[i]==mw_models[5]) pred_midpoints$pred_mid[i]<-exp((boot_Report6$beta_depth[1])*log(pred_midpoints$depth[i])+((boot_Report6$beta_s[pred_midpoints$loc])[i]+(boot_Report6$beta[2]))*log(pred_midpoints$midpoints[i]))
        if (pred_midpoints$model[i]==mw_models[5]) pred_midpoints$pred_mid[i]<-exp((boot_Report9$beta_depth[1])*log(pred_midpoints$depth[i])+(boot_Report9$beta_s[pred_midpoints$loc])[i]+((boot_Report9$beta_b_s[pred_midpoints$loc])[i]+(boot_Report9$beta[2]))*log(pred_midpoints$midpoints[i]))
      }
      
      long_heights$cut_heights<-cut(long_heights$heights/100,breaks=c(midpoints-0.025,1.675))
      agg_long_heights<-aggregate(heights~cut_heights+tow_id,data=long_heights,FUN=length)
      agg_long_heights$mid<-midpoints[as.integer(agg_long_heights$cut_heights)]
      join_frame<-data.frame(tow_id=rep(1:234,each=length(midpoints)),mid=rep(midpoints,234))
      agg_long_heights<-left_join(join_frame,agg_long_heights,by=c("tow_id","mid"))
      agg_long_heights$heights[which(is.na(agg_long_heights$heights))]<-0
      colnames(agg_long_heights)[4]<-"total"
      
      temp_agg_long_heights<-data.frame(agg_long_heights,off=rep(NA,nrow(agg_long_heights)),
                                        off_ln=rep(NA,nrow(agg_long_heights)),insh=rep(NA,nrow(agg_long_heights)),
                                        spat=rep(NA,nrow(agg_long_heights)),spat_off=rep(NA,nrow(agg_long_heights)),
                                        spat_both=rep(NA,nrow(agg_long_heights)))
      
      for (model in 1:length(mw_models)){
        for (i in 1:nrow(gb_2023_surv)){
          if (model %in% c(1,2)){
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
                                         spat_both=aggregate(spat_off~tow_id,data=temp_agg_long_heights,FUN=sum)$spat_both)
      
      all_sums<-data.frame(tow=rep(agg_midpoints_scallops$tow_id,length(mw_models)),
                           val=c(agg_midpoints_scallops$off,agg_midpoints_scallops$off_ln,
                                 agg_midpoints_scallops$insh,agg_midpoints_scallops$spat,
                                 agg_midpoints_scallops$spat_off,agg_midpoints_scallops$spat_both),
                           model=rep(mw_models,each=nrow(agg_midpoints_scallops)),
                           geometry=rep(gb_2023_surv$geometry,length(mw_models))) %>% st_as_sf()
      all_sums$val<-all_sums$val/1000
      
      calc_tot<-rep(NA,length(mw_models))
      
      for (model in 1:length(mw_models)){
        calc_tot[model]<-mean((st_drop_geometry(all_sums)[which(all_sums$model==mw_models[model]),]$val/(actual_dists$dist*2.4384/10^6))*tot_area)/1000
      }
      
      tot_frame[[boot]]<-data.frame(model=mw_models,
                                    est=calc_tot)
      
    }
    
    final_tot_frame<-tot_frame[[1]]
    for (i in 2:nboot){
      final_tot_frame<-rbind(final_tot_frame,tot_frame[[i]])
    }
    
    final_ests<-aggregate(est~model,data=final_tot_frame,FUN=mean)
    final_up_ci<-aggregate(est~model,data=final_tot_frame,FUN=up_quant)
    final_low_ci<-aggregate(est~model,data=final_tot_frame,FUN=low_quant)
    
    final_est_ci[[sim]]<-data.frame(final_ests,up=final_up_ci$est,low=final_low_ci$est)
    
  }
}






