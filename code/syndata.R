source('mcmc_func.R')
sim.GPqp <- function(par,t){
    ind <- grep('sigma.red',names(par))
    sigma.red <- par[ind[1]]
    ind <- grep('\\dl\\d',names(par))
    l <- par[ind[1]]
    ind <- grep('\\dPgp\\d',names(par))
    Pgp <- par[ind[1]]
    ind <- grep('\\dlp\\d',names(par))
    lp <- par[ind[1]]
####covariance
    cov.red.sim = rednoise.cov(sigma.red,l,tt=t,tol=1e-10,Pgp=Pgp,lp=lp,type='qp')
    rv.gp <- mvrnorm(1,rep(0,length(t)),cov.red.sim)
    return(rv.gp)
}
#noise.type <- 'white'
noise.types <- c('W0P','W2P', 'MA0P', 'MA2P', 'GP0P', 'GP2P')
version <- 4
for(noise.type in noise.types){
###load data
    car <- read.table('../data/real/GJ699_carmenes.dat')
    keck <- read.table('../data/real/GJ699_keck.dat')
    harps <- read.table('../data/real/GJ699_harps.dat')
    out <- rbind(car,keck,harps)
    trv.all <- out[,1]
    ind <- sort(trv.all,index.return=TRUE)$ix
    trv.all <- trv.all[ind]
    RV.all <- out[ind,2]
    eRV.all <- out[ind,3]
###simulate noise using the jitters determined for MA(1) for combined data
    if(grepl('W',noise.type)){
        par.opt <- read.table('../data/par/Wpar.txt',header=TRUE,check.names=FALSE)
    }
    if(grepl('MA',noise.type)){
        par.opt <- read.table('../data/par/MApar.txt',header=TRUE,check.names=FALSE)
    }
    if(grepl('GP',noise.type)){
        par.opt <- read.table('../data/par/GPpar.txt',header=TRUE,check.names=FALSE)
    }
    par.opt <- as.matrix(par.opt)[1,]
    prior.type <- 'e0'
    arma.type <- 'abs'
    gp.type <- 'qp'
    Np <- 2
    js <- par.opt[c('1s11','1s12','1s13')]
    jitter.car <- rnorm(nrow(car),0,js[1])
    jitter.harps <- rnorm(nrow(harps),0,js[2])
    jitter.keck <- rnorm(nrow(keck),0,js[3])
    sim.car <- jitter.car
    sim.harps <- jitter.harps
    sim.keck <- jitter.keck
    RV.car <- RV.kepler(pars.kep=par.opt,tt=car[,1]%%2400000,kep.only=TRUE)[[1]]
    RV.harps <- RV.kepler(pars.kep=par.opt,tt=harps[,1]%%2400000,kep.only=TRUE)[[1]]
    RV.keck <- RV.kepler(pars.kep=par.opt,tt=keck[,1]%%2400000,kep.only=TRUE)[[1]]
    if(grepl('2P',noise.type)){
        sim.car <- sim.car +RV.car
        sim.harps<- sim.harps+RV.harps
        sim.keck<- sim.keck+RV.keck
    }
    if(grepl('MA',noise.type)){
                                        #    par.data[[k1]]$RV2 <- par.data[[k1]]$RV-RV.noise[[k1]]
        sim.car <- sim.car+arma(par.opt,rv.kep=RV.car,j3=1,x=car[,2],t=car[,1]%%2400000)
        sim.harps <- sim.harps+arma(par.opt,rv.kep=RV.harps,j3=2,x=harps[,2],t=harps[,1]%%2400000)
        sim.keck <- sim.keck+arma(par.opt,rv.kep=RV.keck,j3=3,x=keck[,2],t=keck[,1]%%2400000)
    }
    if(grepl('GP',noise.type)){
        rv.gp <- sim.GPqp(par=par.opt,t=trv.all)
        ind1 <- match(car[,2],RV.all)
        ind2 <- match(harps[,2],RV.all)
        ind3 <- (1:length(RV.all))[-c(ind1,ind2)]
        sim.car <- sim.car+rv.gp[ind1]
        sim.harps <- sim.harps+rv.gp[ind2]
        sim.keck <- sim.keck+rv.gp[ind3]
    }
    f1 <- paste0('../data/synthetic/GJ699_car',noise.type,version,'.dat')
    cat('output file:\n',f1,'\n')
    write.table(cbind(car[,1],sim.car,car[,3]),file=f1,quote=FALSE,row.names=FALSE,col.names=c('JD','RV','eRV'))

    f2 <- paste0('../data/synthetic/GJ699_harps',noise.type,version,'.dat')
    cat('output file:\n',f2,'\n')
    write.table(cbind(harps[,1],sim.harps,harps[,3]),file=f2,quote=FALSE,row.names=FALSE,col.names=c('JD','RV','eRV'))

    f3 <- paste0('../data/synthetic/GJ699_keck',noise.type,version,'.dat')
    cat('output file:\n',f3,'\n')
    write.table(cbind(keck[,1],sim.keck,keck[,3]),file=f3,quote=FALSE,row.names=FALSE,col.names=c('JD','RV','eRV'))
}
