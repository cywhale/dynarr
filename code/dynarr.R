#' dynarr: A package model seedling species dynamics in regeneration processes
#'
#' The dynarr package reconstruct the community assembly of seedlings (recruits) under species-specific and chance deaths.
#' By random shuffling total number of births and deaths, simultae how the assemblage can "recover" and maintain species coexistence
#'
#' Still in coding, debugging, functions in progressing...
#'
#' @section dynarr functions:
#' increment_observ
#' num_distribx
#' sp_distribx
#' sp_wdistribx
#' rassembly_param
#' rassembly_param
#' rassembly_preset
#' rassembly_dyna

#' Create species dynamics (along time-lines) for births, deaths, and increment from raw-data
#'
#' @param sdl A data.table with species dynamics, including columns:...
#' @param varname A character vector: indicate column names in sdl
#' @return A list include data.table for unified sdl, and long data.table with time-lines of species dynamics
#' @examples
#' library(NanhsiDT)
#' library(data.table)
#' sdl<- copy(compo_sdl_dyna13)
#' increment_observ(sdl)
#' @rdname increment_observ
#' @export
increment_observ <- function(sdl,varname=c("spcode","recruit","recruit_time","dead_time")) {

  sdlx <- sdl %>% data.table() %>%
    .[,`:=`(spcode=paste0(get(varname[1])),
            recruit=as.integer(get(varname[2])),
            recruit_time=as.IDate(get(varname[3])),
            dead_time=as.IDate(get(varname[4]))#,
            #indv=1L
            )]

  nrec <- sdlx[recruit==1L,] %>% .[,indv:=1L] %>%
    data.table::dcast(recruit_time~spcode,fun.aggregate=sum, value.var="indv") %>%
    setnames(1,"times")
  drec <- sdlx[!is.na(dead_time),] %>% .[,indv:=-1L] %>%
    data.table::dcast(dead_time~spcode,fun.aggregate=sum, value.var="indv") %>%
    setnames(1,"times")

  #for(j in colnames(drec)[2:ncol(drec)]) set(drec,j=j, value=ifelse(drec[[j]]==0,0,-drec[[j]])) ## will minus the dead count

  if (nrow(sdlx[recruit==0L,])>0) {
    xini <- sdlx[recruit==0L,] %>% .[,indv:=1L] %>% # old seedling existed before first survey
      data.table::dcast(recruit_time~spcode,fun.aggregate=sum, value.var="indv") %>%
      setnames(1,"times")

    timex<- sort(unique(c(xini$times, nrec$times, drec$times)))
    xcur <- merge(xini,data.table(times=timex),by="times",all.y=T)
    for(j in colnames(xcur)[2:ncol(xcur)]) set(xcur,j=j, value=ifelse(is.na(xcur[[j]]),0,xcur[[j]]))

    xcur <- ladd_merge(list(xcur, nrec, drec), by="times")
  } else { ## if no old seedlings..

    timex<- sort(unique(c(nrec$times, drec$times)))
    xcur <- merge(nrec,data.table(times=timex),by="times",all.y=T)
    for(j in colnames(xcur)[2:ncol(xcur)]) set(xcur,j=j, value=ifelse(is.na(xcur[[j]]),0,xcur[[j]]))

    xcur <- ladd_merge(list(xcur, drec), by="times")
  }
  for(j in colnames(xcur)[2:ncol(xcur)]) set(xcur,j=j, value=ifelse(is.na(xcur[[j]]),0,xcur[[j]]))

  ccsum<- apply(xcur[,-1,with=F],2,cumsum)
  ctsum<- apply(ccsum,1,sum)
  #check
  #ctsum[length(ctsum)]
  #420
  #nrow(sdlx)-nrow(sdlx[!is.na(dead_time),])
  #420 !! match..

  xcurx<- cbind(xcur[,.(times)], data.table(ccsum), data.table(total=ctsum)) %>%
    data.table::melt(id.vars="times",measure.vars=colnames(.)[2:ncol(.)],
                     variable.name="spcode") %>% .[,class:="Survival"]
  nrecx<- cbind(nrec,data.table(total=apply(nrec[,-1,with=F],1,sum))) %>%
    data.table::melt(id.vars="times",measure.vars=colnames(.)[2:ncol(.)],
                     variable.name="spcode") %>% .[,class:="Recruit"]
  drecx<- cbind(drec,data.table(total=apply(drec[,-1,with=F],1,sum))) %>%
    data.table::melt(id.vars="times",measure.vars=colnames(.)[2:ncol(.)],
                     variable.name="spcode") %>% .[,class:="Dead"] %>%
    .[,value:=abs(value)]

  return(list(sdlx=sdlx, incdt=rbindlist(list(xcurx,nrecx,drecx))))
}


############################# Get distribution from discrete bins
num_distribx <- function(x, bins=1000L, seed=NA, min_occur=NA_integer_, max_occur=NA_integer_) {
  x <- na.omit(x)
  ntdx <- approxfun(density(x))

  if (!is.na(seed)) {set.seed(seed)}

  minx <- ifelse(is.na(min_occur),range(x)[1], as.integer(min_occur))
  maxx <- ifelse(is.na(max_occur),range(x)[2], as.integer(max_occur))

  tt <- sample(minx:maxx, bins, replace=TRUE)
  tt1<- ntdx(tt)
  ttx<- as.integer(tt1*round(1/min(tt1),0))

  newntd <- do.call(rep,args=list(x=tt,times=ttx)) %>%
    .[sample(seq_along(.),length(.),replace=F)] %>%
    .[sample(seq_along(.),length(.),replace=F)] ## some continuous rep() too long, may not random just one-time sample
  return(newntd)
}

sp_distribx <- function(data, splist, varname=c("spcode"), basebins=22L,
                        seed=NA, min_occur=NA_integer_, smooth_chance=FALSE) { ## chance death will let non-occurrence-death species have chance as 0.5*(min of set), otherwise, they will never died
  nsp <- data %>% data.table() %>%
    .[,`:=`(spcode=paste0(get(varname[1])))] %>%
    .[,.N, by=.(spcode)] %>% setnames(2,"spn")

  minx <- ifelse(is.na(min_occur), 1L, as.integer(min_occur))
  basen<- ifelse(is.na(basebins) | basebins<2L, 2L, basebins)

  nsp <-nsp[,spn:=spn*(as.integer(basen)-1L)] %>%
    rbind(data.table(spcode=splist[!splist %in% nsp$spcode]) %>% .[,spn:=minx]) %>%
    #.[,spcode:=paste0(spcode)] %>%
    setorder(spcode)
  #all.equal(nsp$spcode, spxt)
  #TRUE

  if(smooth_chance==TRUE) {
    nsp[,spn:=as.integer(round(10*sqrt(spn),0))]
  }

  if (!is.na(seed)) {set.seed(seed)}

  newsp <- do.call(rep,args=list(x=seq_along(splist),times=nsp$spn)) %>%
    .[sample(seq_along(.),length(.),replace=F)] %>%
    .[sample(seq_along(.),length(.),replace=F)] ## some continuous rep() too long, may not random just one-time sample
  #table(newsp) ## check
  return(newsp)
}

sp_wdistribx <- function(x, weights, maxbins=1000L, seed=NA) {

  nw <- weights/sum(weights)

  if (!is.na(seed)) set.seed(seed)

  ttx<- as.integer(round(nw*(1/min(nw)),0))
  ttx<- ttx*as.integer(maxbins/sum(ttx))

  newntd <- do.call(rep,args=list(x=seq_along(x),times=ttx)) %>%
    .[sample(seq_along(.),length(.),replace=F)] %>%
    .[sample(seq_along(.),length(.),replace=F)] ## some continuous rep() too long, may not random just one-time sample

  return(newntd)
}

rassembly_cfg <- function(file) {
  cfg1 <- fread(file,sep="\t",header = FALSE) %>% setnames(1:3,c("param","value","comment"))
  return(cfg1[,1:2])
}

rassembly_param <- function(var,extset,cfgset,defval) {
  #if(!is.null(extset)) {
  if(!is.na(extset)) {
    return(extset)
  }
  #}
  if(nrow(cfgset)>0) {
    if(any(cfgset$param %in% var)) {
      vt <- cfgset[param==var,]$value
      if(!is.na(vt) & is.numeric(vt)) {
        return(vt)
      }
    }
  }
  return(defval)
}

rassembly_preset <- function(sdl,varname=c("spcode","recruit","recruit_time","dead_time","d_cause"),
                             ctrl_dcause=c(), cfgfile=NA_character_,...) {
  dx <- increment_observ(sdl, varname=varname)
  dx$sdlx[,`:=`(d_cause=get(varname[5]))]

  otherd <- setdiff(sort(unique(dx$sdlx$d_cause)), c(0,ctrl_dcause))

  xddh <- dx$sdlx[!is.na(dead_time) & d_cause %in% ctrl_dcause,]   ##pool for seedlings damaged by herbivory
  xddn <- dx$sdlx[!is.na(dead_time) & d_cause %in% otherd,] ##pool for seedlings dammaged by other causes

  ##dynamic model configuration #######################################################################
  if(!is.na(cfgfile)) {
    fcfg <- rassembly_cfg(cfgfile)
  } else {
    fcfg <- data.table()
  }

  #param <- c("MaxDeadPeak","Appa_ctl"))
  #par <- lapply(param, function(x,fcfg,defv,...) {
  #  ifelse(hasArg(substitute(x)), Appa_ctl, NA)
  #})
  function(x,...) {
    ifelse(hasArg(substitute(x)), get(deparse(substitute(x))), NA)
  }

  MaxDeadPeak <- rassembly_param("MaxDeadPeak",ifelse(hasArg(MaxDeadPeak), MaxDeadPeak, NA), fcfg,
                                 defval=as.integer(nrow(dx$sdlx[is.na(dead_time)])*0.2)) ## too big dead number at one time at simulation will crash all system. Not reseaonble for dead caused by herbivory
  Appa_ctl <- rassembly_param("Appa_ctl",ifelse(hasArg(Appa_ctl), Appa_ctl, NA), fcfg, 0.1)  #when,r=3, I'll turn off Flunc_flag ## set NA let apparent competition not occur # set value as minimum occurred rate
  Flunc_flag <- rassembly_param("Flunc_flag",ifelse(hasArg(Flunc_flag), Flunc_flag, NA), fcfg, 1L)
  Min_interN <- rassembly_param("Min_interN",ifelse(hasArg(Min_interN), Min_interN, NA), fcfg, 3L) ## at least N >= Min_interN, interaction will act, otherwise, ineraction on the sp is zero
  Min_interSP<- rassembly_param("Min_interSP",ifelse(hasArg(Min_interSP), Min_interSP, NA), fcfg, 2L) ## at least specise >= Min_interSP, interaction will act, otherwise, ineraction will stop

  seed <- rassembly_param("seed",ifelse(hasArg(seed), seed, NA), fcfg, NA) ## for reproducible results, you can change it anyway
  B <- rassembly_param("B",ifelse(hasArg(B), B, NA), fcfg, 1000L)
  CriticalLoss <- rassembly_param("CriticalLoss",ifelse(hasArg(CriticalLoss), CriticalLoss, NA), fcfg, 0.5)  ## Critical Loss rate of species number to check (0.5 is Half Loss)
  Nct <- rassembly_param("Nct",ifelse(hasArg(Nct), Nct, NA), fcfg, 0.5) ## negative compensatory effect #make mean not larger than 0.35 #see 01_monthly_birth_death.R
  VarNDD_flag <- rassembly_param("VarNDD_flag",ifelse(hasArg(VarNDD_flag), VarNDD_flag, NA), fcfg, 1L) ## only test when Appa_ctl[1] (r=1), turn-off variety in density-denpendent deaths, fixed at Nct
  ctrl_odd <- rassembly_param("ctrl_odd",ifelse(hasArg(ctrl_odd), ctrl_odd, NA), fcfg, 0.67) ## control of d_cause, percentage of herbivory dead-cause
  dtriBins <- rassembly_param("dtriBins",ifelse(hasArg(dtriBins), dtriBins, NA), fcfg, 1000L)## control of d_cause, percentage of herbivory dead-cause
  rare_chance <- rassembly_param("rare_chance",ifelse(hasArg(rare_chance), rare_chance, NA), fcfg, 0.59) #0.53 ## give rare species more possibility to have a recruit or occationality of death

  ###############################################################################################
  ## get distributions of dead number and recruit nmuber
  dtrib <- num_distribx(dx$incdt[class=="Dead" & spcode=="total",]$value,
                        bins=dtriBins, seed=seed,
                        min_occur = 0L, max_occur = MaxDeadPeak)
  rtrib <- num_distribx(dx$incdt[class=="Recruit" & spcode=="total",]$value,
                        bins=dtriBins, seed=seed, min_occur = 0L)

  splist<- unique(dx$sdlx[,.(spcode)]) %>% unlist(use.name=F) %>% paste0 %>% sort()

  return(list(sdlx=dx$sdlx, incdt=dx$incdt, ctrldx=xddh, otherdx=xddn,
              dtrib=dtrib, rtrib=rtrib, splist=splist,
              par=list(MaxDeadPeak=MaxDeadPeak,
                       Appa_ctl=Appa_ctl,
                       Flunc_flag=Flunc_flag,
                       Min_interN=Min_interN,
                       Min_interSP=Min_interSP,
                       seed=seed, B=B, CriticalLoss=CriticalLoss,
                       Nct =Nct, VarNDD_flag=VarNDD_flag,
                       ctrl_odd=ctrl_odd,
                       dtriBins=dtriBins,
                       rare_chance=rare_chance)))
}

rassembly_dyna <- function (dyna_set) {
  dnum <- sample(dyna_set$dtrib, dyna_set$par$B, replace=FALSE) ## dead number at each time of simulations
  rnum <- sample(dyna_set$rtrib, dyna_set$par$B, replace=FALSE) ## recruit number at each time of simu

  dhnum<- round(dnum*dyna_set$par$ctrl_odd, 0) ## observed value is 0.67, roughly 0.67:0.33 = 2:1 (D1+D4) vs (D2+D3+D5) of dead causes
  dnnum<- dnum-dhnum

  spxt <- dyna_set$splist
  newsp<- sp_distribx(dyna_set$sdlx[recruit==1L,], spxt)#, seed=seed)    ## get distributions of species ratio of recruitments
  dhsp <- sp_distribx(dyna_set$ctrldx, spxt)#, seed=seed)    ## get distributions of species ratio of dead by (D1+D4)
  dnsp <- sp_distribx(dyna_set$otherdx,spxt,  #seed=seed,
                      smooth_chance=TRUE) #min_occur = TRUE, ..of dead by (D2+D3+D5)

  tcur <- dyna_set$sdlx[is.na(dead_time),.N,by=.(spcode)] %>% setnames(2,"N0") %>% ## time zero, initial value
    rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N0:=0L]) %>% .[,spcode:=paste0(spcode)] %>%
    setorder(spcode) %>% .[,nt:=N0] %>%
    .[,sdseed:=ifelse(N0==0,dyna_set$par$rare_chance,#0.9*
                      log(1+N0))] #(1-1/22 ~0.96 0.96*0.53~ 0.5 have prob=(1-0.96) to round to 1)
  ######################################## sdseed is for generating random noises for dead/recruit of each sp based on their current population size
  ######################################## we want to increase small variations in rare sp, because of their rare observationsfor dead causes or recruit events.

  appsp <- sp_wdistribx(x=spxt, weights=tcur$sdseed) ## for apparent competition

  nv <- matrix(0,nrow=length(spxt),ncol=dyna_set$par$B) ## total numbers variation in species, for debugging
  nv[,1] <-tcur$N0

  #curspn <- c(nrow(sdlx[is.na(dead_time),.N,by=.(spcode)]))
  spratio <- rep(NA_real_,dyna_set$par$B)
  spappct <- rep(0,length(spxt))
  spratio[1]<-nrow(dyna_set$sdlx[is.na(dead_time),.N,by=.(spcode)])/length(spxt) #curspn[1])  ### lOss rate against whole poor or initial pool
  #halfit <- NA_integer_;
  #sploss <- c(); ## for debugging
  flag<- FALSE
  ndd <- rep(0L, length(spxt))

  for(i in 1:(dyna_set$par$B-1L)) {
    ## NDD is the effect in next iteration (one-lagged)
    ## but random noise and apprant cometition follow current or follow previoust-time increment distribution?
    ## random fluctation can be random on/off, not always occur if <0
    ## if >0 it's occassional recruitment for sp, especially useful for rare sp, otherwise rare sp would not change
    ## tcur$sdseed is more slow, not so big-contrast, not sharp distrib between abund and rare sp
    ## NDD, negative density-dependency
    if (i>1L) {
      if (dyna_set$par$VarNDD_flag==1L) {
        ndr <- runif(length(spxt),0,dyna_set$par$Nct)
      } else {
        ndr <- rep(dyna_set$par$Nct*0.5, length(spxt))
      }
      #ndr <- ndr + runif(length(spxt),NDD_Old*(Nct-ndr),(Nct-ndr))*tcur$nt/max(tcur$nt)
      ndd <- -1L*as.integer(round(ndr*ddi$N,0)) #Neg Density-Dependency
      ndd[ddi$N<=1L | tcur$nt<dyna_set$par$Min_interN] <- 0L  ## only one recruit, too few indv, no density effect

      ## random noise or follow previoust-time increment distribution?
      ## random fluctation can be random on/off, not always occur if <0
      ## if >0 it's occassional recruitment for sp, especially useful for rare sp, otherwise rare sp would not change
      ## tcur$sdseed is more slow, not so big-contrast, not sharp distrib between abund and rare sp
      #if (i>1) {
      #ddi[which(ddi$N<=1 | tcur$nt<Min_interN), N:=0] ## only one recruit, too few indv, no density effect

      #ndd <- ndd + fluc ## now may increase or decrease, but may cause bug in next step, so move to later

      #maxd<-round(rnum[i-1]*Nct*0.5,0)-sum(abs(ndd))-sum(abs(fluc[fluc<=0])) #If NDD+random dead > mean NCT dead, should stop
      maxd<- round(sum(ddi$N)*dyna_set$par$Nct*0.5,0)-sum(abs(ndd))
      if (!is.na(dyna_set$par$Appa_ctl) & maxd>=2L) { ## apparent competition must work on at least two sp, and thus two deads
        dselt<- which(ddi$N>=1L & #(ddi$N>0 | fluc>0) &
                      tcur$nt>=dyna_set$par$Min_interN) ## work on rare sp! #those population increase had apparent competition with each others
        if (length(dselt)>=dyna_set$par$Min_interSP) {
          appd <- as.integer(round(runif(1,dyna_set$par$Appa_ctl*maxd,maxd),0))
          selapp<-spxt[sample(appsp[appsp %in% dselt],appd,replace=TRUE)]

          #if (length(dselt)>=Min_interSP & appd>=2) { ########## apparent competetion occur if sp_num >=2, sum(dead)>=2
          if (length(unique(selapp))>=dyna_set$par$Min_interSP & appd>=2) { ## apparent competetion occur if sp_num >=2, sum(dead)>=2
            appalpha<- min(tcur[spcode %in% selapp,]$sdseed)/sum(tcur[spcode %in% selapp,]$sdseed)
            #sdseed already log(nt+1)!! and set a minimum rule (rare chance) to prevent near zero value
            #ddalpha <- sapply(selapp, function(x) {sum(log(tcur[spcode %in% selapp[selapp!=x],]$nt+1))/(1+appalpha*sum(log(tcur[spcode %in% selapp,]$nt+1)))}) %>%
            ddalpha <- sapply(selapp, function(x) {sum(tcur[spcode %in% selapp[selapp!=x],]$sdseed)/(1+appalpha*sum(tcur[spcode %in% selapp,]$sdseed))}) %>%
              data.table(spcode=names(.)) %>% setnames(1,"da") %>%
              .[,{.(da=sum(da))}, by=spcode] ## ratios in repeated species in selapp should sum together first.

            #appdt<- data.table(spcode=selapp, N=as.integer(round(-1*appd*ddalpha/sum(ddalpha),0))) %>% #bug
            appdt<- ddalpha %>% .[,{.(N=as.integer(round(-1*appd*da/sum(.$da),0)))},by=spcode] %>%
              #.[,{.(N=sum(N))},by=.(spcode)] %>%
              rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L]) %>%
              setorder(spcode)

            selappt <- which(appdt$N!=0L)
            if (any(selappt)) {
              spappct[selappt] <- spappct[selappt] + rep(1L, length(selappt))
            }

            #Old version, rare sp hurted least
            #appdt<- data.table(spcode=spxt[sample(appsp[appsp %in% dselt],appd,replace=TRUE)]) %>% ##use appsp, weighted by tcur$sdseed #use newsp or 1:length(spxt)??
            #.[,-1*(.N), by=.(spcode)] %>% setnames(2,"N") %>%
            #rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L]) %>%
            #setorder(spcode)

            #appdt[which(tcur$nt<=1), N:=0] ## if indiv <=1 no bio-interaction with other species
            ndd <- ndd + appdt$N
          } #else {
        } #else {
      } #else {
      #ndd <- ndd + fluc
    }

    ## D1+D4 (herbivory ratio death)
    #  ddi <- xddh[sample(1:nrow(xddh), dhnum[i], replace=FALSE),.(spcode)] %>% .[,.N, by=.(spcode)] %>%
    #    rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L])
    if(any(which(tcur$nt==0L))) {
      delsp <- match(tcur[nt==0L,]$spcode, spxt)
      dselt <- dhsp[!dhsp %in% delsp]
    } else {
      dselt <- dhsp
    }
    if (length(dselt)>dhnum[i]) {
      ddi <- data.table(spcode=spxt[sample(dselt, dhnum[i], replace=FALSE)]) %>% .[,.N, by=.(spcode)] %>%
        rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L])

      tcur[ddi, nt:=nt-N, on=.(spcode)]

      iter <- 1L
      while (any(which(tcur$nt<0))) { ## negative values mean we had some excess dead on some species
        delsp <- match(tcur[nt<0,]$spcode, spxt)
        dselt <- dselt[!dselt %in% delsp]
        tt <- sum(abs(tcur[nt<0,]$nt))
        tcur[nt<0,nt:=0L]

        if (length(dselt) < tt) break
        ddi <- data.table(spcode=spxt[sample(dselt, tt, replace=FALSE)]) %>% .[,.N, by=.(spcode)] %>%
          rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L])

        tcur[ddi, nt:=nt-N, on=.(spcode)]
        iter <- iter+1L
        if (iter>=10L) {cat("Routine to find excess dead not stop? ", iter, " at i = ",i, "\n"); break}
      }
    }

    if(any(which(tcur$nt==0))) {
      delsp <- match(tcur[nt==0,]$spcode, spxt)
      dselt <- dnsp[!dnsp %in% delsp]
    } else {
      dselt <- dnsp
    }

    ## D2+D3+D5 other_cause ratio death
    if (length(dselt)>dnnum[i]) {
      ddi <- data.table(spcode=spxt[sample(dselt, dnnum[i], replace=FALSE)]) %>% .[,.N, by=.(spcode)] %>%
        rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L])

      tcur[ddi, nt:=nt-N, on=.(spcode)]

      iter <- 1L
      while (any(which(tcur$nt<0))) { ## negative values mean we had some excess dead on some species
        delsp <- match(tcur[nt<0,]$spcode, spxt)
        dselt <- dselt[!dselt %in% delsp]
        tt <- sum(abs(tcur[nt<0,]$nt))
        tcur[nt<0,nt:=0L]

        #if (length(dselt) < tt | length(unique(dselt)) < Min_interSP) break
        if (length(dselt) < tt) break
        ddi <- data.table(spcode=spxt[sample(dselt, tt, replace=FALSE)]) %>% .[,.N, by=.(spcode)] %>%
          rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L])

        tcur[ddi, nt:=nt-N, on=.(spcode)]
        iter <- iter+1L
        if (iter>=10L) {cat("Routine to find excess dead not stop? ", iter, " at i = ",i, "\n"); break}
      }
    }

    ## recruitment (for next turn)
    ddi <- data.table(spcode=spxt[sample(newsp, rnum[i], replace=FALSE)]) %>% .[,.N, by=.(spcode)] %>%
      rbind(data.table(spcode=spxt[!spxt %in% .$spcode]) %>% .[,N:=0L]) %>%
      setorder(spcode) %>%
      .[,N:=N+ndd]

    tcur[ddi, nt:=nt+N, on=.(spcode)] ## nv[,1] already is N0, so nt is actually for next turn

    ## randon noise + population effect will have a negative compensatory effect in next turn
    ## if not add here, but add in nt:=nt+N+snoise, will cause rare sp to increase, because they also rare dead events
    fluc<-  #random on/off fluctuation
      as.integer(round(runif(length(spxt),0,1),0)*              ## random turn-off, otherwise too big variation
                 round(runif(length(spxt),-1,1)*tcur$sdseed,0)) ## random fluctuation

    if (dyna_set$par$Flunc_flag==1L) { ## Turn on/off Flunc_flag, not make 0, but let the fluctuations syn with recruitment burst (otherwise, almost no rare sp born)
      fluc[(ddi$N>=1L & fluc>0L & tcur$nt>=dyna_set$par$Min_interN) | (tcur$nt<abs(fluc) & fluc<0L)] <- 0L ############# dont let continuous fluctuation
    } else {
      ## want to make recruitment sync (nt>=Min_interN for those rarely recruit, to let him occasionally had births)
      fluc[(ddi$N<1L & fluc>0L & tcur$nt>=dyna_set$par$Min_interN) | #### but indeed make abundant sp with fluc, but rare sp no fluc
             (tcur$nt<abs(fluc) & fluc<0L)] <- 0L
    }

    ddi[,N:=N+fluc] ### for next iteration's NDD used
    tcur[,nt:=nt+fluc]

    curspn <- length(spxt)-length(which(tcur$nt<=0))
    if (curspn<=0) {
      #break
      cat("TOTAL 100% loss of species at time: ", i, " and iteration: ", m, " and Param: ", k,"\n")
    }
    spratio[i+1] <- curspn/length(spxt) #curspn[1])

    if (spratio[i+1]<=dyna_set$par$CriticalLoss) { # curspn[1]
      #cat("Iteration: ",i, " and current species number: ", curspn, "\n")
      if (flag==FALSE) {
        cat(paste0((1-dyna_set$par$CriticalLoss)*100.0),"% loss of species at time: ", i, "\n")# and iteration: ", m, " and Param: ", k,"\n")

        halfit <- i #halfit[k,m]<- i
        flag <- TRUE
      }
    }
    #cat("Iteration: ",i, " and current species number: ", curspn[i+1], "\n")
    #for debugging
    #tt <- which(unlist(nv[,i])>0 & unlist(tcur$nt)<=0)
    #if (any(tt)) {
    #  cat("Iteration: ",i, " and cuddirrent species number: ", curspn[i+1], "\n")
    #  cat("New loss species: ",paste(spxt[tt],collapse = ","),"\n")
    #  sploss <- c(sploss,paste(spxt[tt],collapse = ","))
    #} else {
    #  sploss <- c(sploss,NA_character_)
    #}

    tcur[nt<0,nt:=0L] ## dont leave negative value
    tcur[,sdseed:=ifelse(nt==0,dyna_set$par$rare_chance,log(1+nt))] ## adjust SD according to current population
    appsp <- sp_wdistribx(x=spxt, weights=tcur$sdseed)
    nv[,i+1] <- tcur$nt
  }
  #debugging
  #cat(paste0((1-CriticalLoss)*100.0),"% loss iteration is: ", halfit, "\n")

  #if (Plot==1L) {
  #  par(mfcol=c(3,1))
  #  plot(1:dyna_set$par$B,spratio,type="l")
  #  chkt <- cbind(cbind(1:dyna_set$par$B),cbind(t(nv[c(8,16,18,20,31),])))
  #  matplot(chkt[,1],chkt[,2:6], type="l", col=1:5)
  #  legend(0.8*dyna_set$par$B, 0.8*range(nv[c(8,16,18,20,31),])[2],
  #         legend=spxt[c(8,16,18,20,31)],col=1:5,lty=rep(1,5))

  #  chkt <- cbind(cbind(1:dyna_set$par$B),cbind(t(nv[c(2,21),])))
  #  matplot(chkt[,1],chkt[,2:3], type="l", col=1:2)
  #  legend(0.8*dyna_set$par$B, 0.8*range(nv[c(2,21),])[2],
  #         legend=spxt[c(2,21)],col=1:2,lty=rep(1,2))
  #}

  #for (j in seq_along(Chk_SP)) {
  #  tt <- get(paste("n",Chk_SP[j],k,r,sep="_"), envir = yk)
  #  tt[m,] <- nv[match(Chk_SP[j], spxt),]
  #  assign(paste("n",Chk_SP[j],k,r,sep="_"),tt, envir = yk)
  #}
  ## store species abundance at last 100 runs
  #tt <- get(paste("spabund",k,r,sep="_"), envir = yk)
  #tt[,m] <- rowMeans(nv[,(B-99L):B],na.rm=T)
  #assign(paste("spabund",k,r,sep="_"),tt, envir = yk)
  return(list(nv=nv,spratio=spratio,halfit=halfit,spappct=spappct,tcur=tcur))
}


