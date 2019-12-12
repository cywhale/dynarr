#### This file is to plot a increse trend of ungulate herbivory outside fence
#### data come from R/seedling/survival/00x_init_herbivory01.R L1-L404 
#### save(ihx,nht, times, file="data_init_herbivory01.RData")
#### use random shuffle on fenced/unfenced individuals to boostrap this odds ratio

library(data.table)
library(magrittr)
library(ggplot2)
load("data_init_herbivory01.RData")

#nht can be derived from iht, check the original code
#nht
eval_oddx <- function(dt, times, returnx="odd_ratio") {
  require(data.table)
  require(magrittr)
  count_cy.02 <- function(x) { ## only count postive item
    return (length(which((!is.na(x))&(x>0))))
  }
  
  nht <- dt %>% .[,.(
    nh1=count_cy.02(hurt1)/length(which(!is.na(recruit_time) & 
                                        (as.integer(recruit_time) <= 
                                           as.integer(as.IDate(times[2]))) & 
                                        (is.na(dead_time) |
                                           as.integer(dead_time) >= 
                                           as.integer(as.IDate(times[2]))))),
    nh2=count_cy.02(hurt2)/length(which(!is.na(recruit_time) & 
                                        (as.integer(recruit_time) <= 
                                           as.integer(as.IDate(times[3]))) & 
                                        (is.na(dead_time) |
                                           as.integer(dead_time) >= 
                                           as.integer(as.IDate(times[3]))))),
    nh3=count_cy.02(hurt3)/length(which(!is.na(recruit_time) & 
                                        (as.integer(recruit_time) <= 
                                           as.integer(as.IDate(times[4]))) & 
                                        (is.na(dead_time) |
                                           as.integer(dead_time) >= 
                                           as.integer(as.IDate(times[4])))))
  ), by=.(fence)]

  if (returnx!="odd_ratio") {
    return(nht)
  }

  return(
    #odds ratio
    unlist(nht[fence=="N",2:4,with=F])/unlist(nht[fence=="F",2:4,with=F])
  )
}

################################################# Original case (observed)

dx <- ihx[,.(hurt1,hurt2,hurt3,recruit_time,dead_time,fence)]
tt = eval_oddx(dx, times, returnx="")
### use only .N
#   fence       nh1       nh2        nh3
#1:     N 0.1936620 0.1690141 0.16549296
#2:     F 0.2347826 0.1434783 0.07391304
### use alive as total num
#fence       nh1       nh2        nh3
#1:     N 0.2208835 0.1825095 0.17279412
#2:     F 0.2887701 0.1601942 0.07834101

all.equal(tt, nht)

tt1 <- eval_oddx(dx, times)
###nh1       nh2       nh3 
###0.7370965 1.1749049 2.2936581
class(tt1)
# [1] "numeric"
#################################################
################################################# bootstrapping (no parallel in this case)
library(NanhsiDT)
B = 1000
L = nrow(dx)
fout <- matrix(0,nrow= B, ncol=3) ## 1000 times randomization for fenced
nout <- matrix(0,nrow= B, ncol=3) ## 1000 times randomization for unfenced
randidx <- matrix(0, nrow=L,ncol=B)
randidx <- shuffleM_idx(c(1:L), B, seed = 123L, keep_origin = TRUE, index_only = FALSE)

for(i in 1:B) {
  dxt <- copy(dx) %>% .[,fence:=fence[randidx[,i]]]
  tt1 <- eval_oddx(dxt, times, returnx = "")
  fout[i,] <- unlist(tt1[2,2:4], use.names = FALSE)
  nout[i,] <- unlist(tt1[1,2:4], use.names = FALSE)
}

fx <- as.data.table(nout/fout) %>% setnames(1:3, c("grp1","grp2","grp3")) %>%
  .[, `:=`(fence="fenced", type="random")] %>% .[1, type:="observed"] 

fx[, id:=.I]
dfm <- melt(fx, id.vars = c("id", "fence", "type"), 
            measure.vars = c("grp1", "grp2", "grp3"),
            variable.name = "group", value.name = "odds_ratio") %>% #damage_ratio
  .[,stage:=fifelse(group=="grp1", times[2],
            fifelse(group=="grp2", times[3], times[4]))]

################################################# plotting output

gh1 <- ggplot(dfm[type!="observed",], aes(stage, odds_ratio)) +
  geom_jitter(col="darkgrey", alpha = 0.5) +
  stat_boxplot(col="black", alpha=0.5, width=0.3) + 
  geom_point(data=dfm[type=="observed",], aes(stage, odds_ratio), 
             col="red", size=4, shape=17) + xlab(NULL)

gh1
