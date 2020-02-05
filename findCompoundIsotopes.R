#!/usr/bin/R
## Written By: H. Paul Benton
## Made on: Febuary 18, 2019
## Edited Feb, 4, 2020
## For: Sulphure code
## Version 1.1.0
options(stringsAsFactors = F)

require(signal) || stop("Cannot load signal package needed.")
require(caTools) || stop("Could not load caTools package which is needed.")
library(stringr)
library(magrittr)
require(xcms) || stop("Couldn't find package xcms\n")

load("cmpList-sulphure.Rda")
cmpList<-cmpListSulphure
# adductSelection="standard"
# minScanTime=3; ppmDevN=7.5; snthresh=2; maxoFilt=1000; plot.prog=F; atom<-"S"

ppmDev<-function (Mr, ppmE = 5) {
    error <- (ppmE/10^6) * Mr
    deviation <- c(Mr + error, Mr - error)
    return(range(deviation))
}

ppm<-function(mm, tm){
    abs(((mm-tm)/tm)*10^6)
}

runSulphurFinder<-function(atom="C", cmpList, minScanTime=3, ppmDevN=20, snthresh=2, maxoFilt=1000, adductSelection=c("all", "standard"), plot.prog=F){
    fl<-list.files(pattern="mzXML|mzData|mzML",
                   ignore.case = T,
                   recursive=T,
                   full.names = T)
    ## run the compound list getter here
    start.time<-Sys.time()
    pkl<-list()
    for(i in 1:length(fl)){
        cat(paste0("Reading ", basename(fl[i]),"\n" ))
        object<-readMSData(fl[i], mode="onDisk", msLevel. = 1)
        cat(paste0("Finding peaks....\n"))
        # pkl[[i]]<-tryCatch({
        #     findCompound(object, cmpList, inputAtom=atom, minScanTime=minScanTime, 
        #                  ppmDevN=ppmDevN, snthresh, maxoFilt, 
        #                  adductSelection, plot.prog)
        #     }, warning = function(w){
        #         warning(paste0("WARNING on file ", fl[i], " -- ", w, "\n"))
        #     }, error = function(e) {
        #         cat(paste0("Skipping file ", fl[i], "\n"))
        #         warning(paste0("ERROR on file ", fl[i], " -- ", e, "\n"))
        #     }
        # )
        pkl[[i]]<-findCompound(object, cmpList, inputAtom=atom, minScanTime=minScanTime, 
                               ppmDevN=ppmDevN, snthresh, maxoFilt, 
                               adductSelection, plot.prog)
    }
    end.time<-Sys.time()
    end.time-start.time
    return(pkl)
}

findCompoundIsotopes<-function(object, pkl, pklEIC,
                               maxLook=1, deltaIso, minScanTime, snthresh, 
                               maxoFilt, cmp.list, ppmDevN=20, addID){
    
    maxLook<-as.integer(maxLook)
    deltaIso<-as.numeric(deltaIso)
    mat.iso<-matrix(nrow=3,ncol=as.numeric(maxLook))
    rownames(mat.iso)<-c("min", "max", "med")
    for(i in 1:maxLook){
        mat.iso[c(1,2),i]<-range(ppmDev({as.numeric(pkl[1,"mz"])+{i*deltaIso}}, ppmDevN))
        mat.iso[3,i]<-{as.numeric(pkl[1,"mz"])+{i*deltaIso}}
    }
    ## this makes a matrix that we know where to look for every possible isotope for this compound.
    
    ans<-sapply(1:ncol(mat.iso), function(isoCol, mat.iso, pkl, object, pklEIC, minScanTime, snthresh, maxoFilt, cmp.list, addID){
        iso.mzr<-mat.iso[,isoCol]
        
        eic.iso<-getCompound(object, iso.mzr[c("min", "max")])
        eic.iso<-caTools::runmean(eic.iso, minScanTime)
        if(all(eic.iso == 0)){
            return(NULL)
        }
        ySmSums<-signal::sgolayfilt(eic.iso, p=minScanTime)
        ## get the scan data once for each maxLook isotope
        ResPeakIsoList<-list()
        addID.iso<-paste0(addID, ";", diff(cmp.list["mass"], iso.mzr["med"]))
        for(i in 1:nrow(pkl)){
            pk<-pkl[i,, drop=T]
            rtr<-as.numeric( c(pk["rtmin"], pk["rtmax"]) )
            scanr<-as.numeric(c(pk["scanMin"], pk["scanMax"]))
            
            ySmSums.iso<-ySmSums[scanr[1]:scanr[2]]
            max <- max(ySmSums.iso)
            maxy <- which.max(ySmSums.iso)
            noise <- mean(ySmSums.iso[ySmSums.iso > 0])
            sn <- ySmSums.iso[maxy]/noise
            
            if(ySmSums.iso[maxy] > 0 && ySmSums.iso[maxy] > snthresh*noise && ySmSums.iso[maxy] > 0){
                mzacc<- iso.mzr["med"] ## this is a lie
                deltamz<-ppm(iso.mzr["min"], iso.mzr["max"])
                
                into <- sum(ySmSums.iso)
                intb <- (sum(ySmSums.iso)-sum(rep(noise, length(scanr[1]:scanr[2])) ))
                maxo <- max(ySmSums.iso)
                if(maxo < maxoFilt){
                    return(NULL)
                }
                # stop("HERE !!!\n")
                rho<-cor(ySmSums.iso, pklEIC[seq(scanr[1], scanr[2])] )^2
                if(is.na(rho) || rho < 0.4){
                    return(NULL)
                }else{
                    ResPeakIsoList[[i]]<- cbind("mz"=mzacc, "deltamz"=deltamz, 
                                                "rtmed"=pk["rtmed"], "rtmax"=pk["rtmax"], "rtmin"=pk["rtmin"],
                                                "scanMin"=scanr[1], "scanMax"=scanr[2],
                                                "into"=into, "intb"=intb, "maxo"=maxo, 
                                                "sn"=sn, "adductID"=addID.iso, "corr"=rho, "isotopolog"=isoCol,
                                                "ID"=cmp.list["ID"], "Name"=cmp.list["Biocyc"], 
                                                "mass"=cmp.list["mass"], "formula"=cmp.list["formula"],
                                                "file"=basename(fileNames(object))
                                                )
                    cat(" . ")
                } 
            }else{
                next
            }
        }
        isoResMat<-do.call("rbind", ResPeakIsoList)
        return(isoResMat)
    },  mat.iso, pkl, object, pklEIC, minScanTime, snthresh, maxoFilt, cmp.list, addID)
    # stop("what...")
    if(is.null(ans)){
        return(NULL)
    } else if(is.list(ans) & is.null(ans[[1]])){
        return(NULL)
    } else if(ncol(ans) == 1){
        return(t(ans))
    }else{
        ResMat<-do.call("rbind", ans)
        rownames(ResMat)<-NULL
        return(ResMat)
    }
}

getCompound<-function(object, mzr){
    mz.vec<-mz(object)
    int.vec<-intensity(object)
    cmp<-sapply(1:length(mz.vec), function(i, mzL, intL, mzr){
        idx<-which(mzL[[i]] > mzr[1] & mzL[[i]]< mzr[2])
        return(sum(intL[[i]][idx]) )
    }, mz.vec, int.vec, mzr)
    return(cmp)
}

findCompound<-function(object, cmpList, inputAtom="C", minScanTime=5, 
                       ppmDevN=20, snthresh=2, maxoFilt=1000, adductSelection, plot.prog=FALSE){
    
    if(!pmatch("OnDiskMSnExp", class(object), nomatch = 0)){
        stop("Not a MSnbase object")
    }
    mode<-unique(polarity(object))
    if(length(mode) > 1) stop("We do not allow polarity switching currently\n")
    if(mode == "1"){
        mode<-"+"
        cat("Loading positive mode adducts\n")
        adduct_list<-read.table("adducts.csv", sep=",", header=TRUE)
        if(adductSelection == "standard") adduct_list<-adduct_list[c(1,2,18,24),]
    } else if (mode == "0"){
        mode<-"-"
        cat("Loading negative mode adducts\n")
        adduct_list<-read.table("adductsNEG.csv", sep=",", header=TRUE)
        if(adductSelection == "standard") adduct_list<-adduct_list[c(2,4,6,8),]
    } else {
        mode<-"+"
        cat("No found polarity defaulting to positive mode adducts\n")
        adduct_list<-read.table("adducts.csv", sep=",", header=TRUE)
        if(adductSelection == "standard") adduct_list<-adduct_list[c(1,2,18,24),]
    }
    
    cl <- makeCluster(35, out="")
    # clusterEvalQ(cl, { source("findCompoundIsotopes.R")})
    clusterEvalQ(cl, {
        require(signal) || stop("Cannot load signal package needed.")
        require(caTools) || stop("Could not load caTools package which is needed.")
        require(xcms) || stop("Couldn't find package xcms\n")
        library(stringr)
        library(magrittr)
    })
    clusterExport(cl, c("findCompoundIsotopes", "ppmDev", "ppm", "getCompound", "minRange",
                        "getAccIsoMass", "getAccurateMass", "findIsotopeDelta"))

    pkl.all<-parApply(cl, cmpList, 1, function(cmp.list, adduct.list, object, ppmDevN, mode,
                                                         maxoFilt, inputAtom, minScanTime, plot.prog){
    # pkl.all<-apply(cmpList,  1, function(cmp.list, adduct.list, object, ppmDevN, mode,
    #                                      maxoFilt, inputAtom, minScanTime, plot.prog) {
    cat(paste0("Compound:", cmp.list["Biocyc"], "\n"))
        pkl<-sapply(1:nrow(adduct_list), function(adductIDNo, cmp.list, object, ppmDevN, mode, adduct_list, 
                                                  maxoFilt, inputAtom, minScanTime, plot.prog){
            # mzr<-ppmDev(as.numeric(cmp.list["mass"])-adduct_list[adductID,"massdiff"], ppmDevN)
            adductID<-adduct_list[adductIDNo,]
            Mplus<-FALSE ## option to be added later
            if(any(adductID["name"] %in% "[M]-" || adductID["name"] %in% "[M]+") && Mplus == FALSE){
                return(NULL)
            }
            if(mode == "+"){
                mzM<-(((as.numeric(cmp.list["mass"])*adductID["nmol"])/adductID["charge"]) 
                      + adductID["massdiff"])
                mzr<-ppmDev(mzM , ppmDevN)
            } else if (mode == "-"){
                mzM<-(((as.numeric(cmp.list["mass"])*adductID["nmol"])/adductID["charge"]) 
                      - adductID["massdiff"])
                mzr<-ppmDev(mzM, ppmDevN)
            }
            ## this is the mass range that we'll be searching for
            
            eic<-getCompound(object, mzr)
            if(max(eic) < maxoFilt){
                return(NULL)
            }
            if(length(eic) <= minScanTime){
                return(NULL)
            }
            rtime<-unlist(rtime(object))
            
            if(plot.prog == T) plot(rtime, eic, type="l", main="Raw Spectra", xlab="RT", ylab="intensity")
            eic<-caTools::runmean(eic, minScanTime)
            
            if(sum(eic) == 0){
                return(NULL)
            }
            ySmSums<-signal::sgolayfilt(eic, p=minScanTime)
            ysum.sm<-ySmSums
            if(plot.prog == T) lines(rtime, ysum.sm, col="red", type="l")
            # ySmSums<-approx(ySmSums, n=max(header(object)[,"seqNum"]))$y
            
            ## make the output lists
            ResEICList<-ySmSums ## this is for the raw eic - isotope tracking
            ResPeakList<-list() ## this is for the raw peaklist
            # if(plot.prog == T) par(mfrow=c(2,4))
            for(q in 1:4){
                # ysum.sm, object, cmp.list, adduct.list, adductID, mode){
                ## find a max of 5 peaks !
                max <- max(ysum.sm)
                maxy <- which.max(ysum.sm)
                noise <- mean(ysum.sm[ysum.sm > 0])
                sn <- ysum.sm[maxy]/noise
                if (ysum.sm[maxy] > 0 && ysum.sm[maxy] > snthresh*noise && ysum.sm[maxy] > 0) {
                    # peakrange <- xcms:::descendZero(ysum.sm, maxy)
                    peakrange <- xcms:::descendMin(ysum.sm, maxy)
                    pwid <- (rtime[peakrange[2]] - rtime[peakrange[1]]) /
                        (peakrange[2] - peakrange[1])
                    
                    mzmed<-getAccurateMass(object, mzr, mzM, maxy)
                    mzacc<-mzmed["accurateMZ"]
                    deltamz<-mzmed["deltaMZ"]
                    
                    into <- sum(eic[peakrange[1]:peakrange[2]])
                    intb <- (sum(eic[peakrange[1]:peakrange[2]])-sum(rep(noise, length(peakrange[1]:peakrange[2])) ))
                    maxo <- max(eic[peakrange[1]:peakrange[2]])
                    if(maxo < maxoFilt){
                        next
                    }
                    if(diff(peakrange) < minScanTime){
                        next
                    }
                    
                    if(plot.prog == T) plot(ysum.sm, xlim=c({peakrange[1]-10}, {peakrange[2]+10}), type="l", main="Raw Spectra")
                    # ResEICList<-ySmSums ## this doesn't need to be a list as we're not taking peakrange
                    ## therefore jsut send the whole thing over once instead of multiple copies of the same thing !memeory!
                    ysum.sm[peakrange[1]:peakrange[2]] <- 0 
                    ySmSums[peakrange[1]:peakrange[2]] <- 0
                    rtmax<-if(is.na(rtime[peakrange[2]])) max(rtime) else rtime[peakrange[2]]
                    rtmin<-if(is.na(rtime[peakrange[1]])) min(rtime) else rtime[peakrange[1]]
                    ResPeakList[[q]]<- cbind("mz"=mzacc, "deltamz"=deltamz, 
                                             "rtmed"=rtime[maxy],
                                             "rtmax"=rtmax, "rtmin"=rtmin,
                                             "scanMin"=peakrange[1], "scanMax"=peakrange[2],
                                             "into"=into, "intb"=intb, "maxo"=maxo, 
                                             "sn"=sn, "adductID"=adductID[,"name"], "corr"=0, "isotopolog"=0,
                                             "ID"=cmp.list["ID"], "Name"=cmp.list["Biocyc"],
                                             "mass"=cmp.list["mass"], "formula"=cmp.list["formula"],
                                             "file"=basename(fileNames(object))
                                             )
                    cat("+")
                }
                if(plot.prog == T) dev.off()
            }
            cat(paste0("\tAdduct: ",adductID["name"],"\n"))
            if(length(ResPeakList) == 0){
                return(NULL)
            }
            
            ResMat<-do.call("rbind", ResPeakList)
            rownames(ResMat)<-NULL
            
            ## we need to do a calculation of the mass deviation that we're going to have here
            # inputAtom<-"S" ## this can be changed later to lookup other atoms
            ## mass deviation calculation
            if(!length(grep(inputAtom, cmp.list["formula"]) >0)){
                return(ResMat)
                ## If the compound does't have the input atom in it's formula 
                ## then return the current running list
            }
            deltaIso<-findIsotopeDelta(inputAtom)
            inputAtomNo<-str_extract(as.character(cmp.list["formula"]), 
                                     paste0("\\.*", inputAtom,"\\d*"))
            maxLook<-max(1, gsub(inputAtom, "", inputAtomNo), na.rm=T)
            # stop("check here...")
            ## isotope finding 
            ## now we can do a targeted peak detetor using the peaklist from above !
            isoPKL<-findCompoundIsotopes(object, pkl=as.data.frame(ResMat), 
                                         pklEIC=ResEICList,
                                         maxLook=maxLook, 
                                         deltaIso=deltaIso, 
                                         minScanTime=minScanTime, 
                                         snthresh=snthresh, 
                                         maxoFilt=maxoFilt, cmp.list, 
                                         ppmDevN = ppmDevN,
                                         addID=adductID[,"name"]
            )
            if(is.null(isoPKL)){
                return(ResMat) 
                ## then the isotope was below the thresholds or corr or soething else...
            }else{
                pkmat<-rbind(ResMat,isoPKL)
                return(pkmat)
            }
        }, cmp.list, object, ppmDevN, mode, adduct_list, maxoFilt, inputAtom, minScanTime, plot.prog)
        # pkl<-do.call("rbind", pkl)
        gc()
        return(pkl)
    }, adduct_list, object, ppmDevN, mode, maxoFilt, inputAtom, minScanTime, plot.prog)
    # pkl.all<-do.call("rbind", pkl.all)
    # rownames(pkl.all)<-NULL
    stopCluster(cl)
    return(pkl.all)
}

findIsotopeDelta<-function(inputAtom){
    Atoms<-data.frame("S"=c("S32"=31.97207100, "S34"=33.96786690), 
                      "C"=c("C12"=12, "C13"=13.003355),
                      "N"=c("N15"=15, "N16"=16),
                      "O"=c("O16"=16, "O18"=18)
                      )
    idx<-which(colnames(Atoms) == inputAtom)
    return(diff(Atoms[,idx]))
}

minRange<-function(rangeA, rangeB){
    minR<-min(min(as.numeric(rangeA)), min(as.numeric(rangeB) ))
    maxR<-min(as.numeric(max(rangeA)), max(as.numeric(rangeB) ))
    return(c(minR, maxR))
}

getAccIsoMass<-function(object, mzSearch, ymax){
    # pk<-peaks(object, ymax)
    pk<-filterAcquisitionNum(object, ymax)
    mzL<-unlist(mz(pk))
    intL<-unlist(intensity(pk))
    names(mzL)<-NULL
    
    # stop("error")
    mzrange<-which(mzL > mzSearch[1] && mzL < mzSearch[2])
    Int.max<-which.max(pkL[])
    acc.mz<-mzL[Int.max]
    delta.mz<-ppm(as.numeric(mzSearch[1]), as.numeric(mzSearch[2]))
    return(c("accurateMZ"=acc.mz, "deltaMZ"=delta.mz))
}

getAccurateMass<-function(object, mzRange, mzM, ymax){
    # pk<-peaks(object, ymax)
    pk<-MSnbase::filterAcquisitionNum(object, as.integer(ymax))
    pk<-unlist(mz(pk))
    names(pk)<-NULL
    # cat(paste0("mzSearch -",mzRange, "\n"))
    idx<-which.min(abs(pk - as.numeric(mzM)))
    # idx<-which(pk >= mzRange[1] && pk <= mzRange[2])
    acc.mz<-pk[idx]
    # cat(paste0("acc.mz -", acc.mz, "\n"))
    delta.mz<-ppm(acc.mz, as.numeric(mzM))
    
    # cat(paste0("delta.mz -",delta.mz, "\n"))
    return(c("accurateMZ"=acc.mz, "deltaMZ"=delta.mz))
}

#########################################################################
########################## extra functions ##############################
#########################################################################
rtNoise<-function (spec, rtime, gap = quantile(diff(spec), 0.9, na.rm=T)) {
    intmean <- mean(spec, na.rm=T)
    rtlen <- diff(range(c(rtime)))
    rtdiff <- diff(rtime)
    gaplen <- sum(rtdiff[rtdiff > gap])
    weighted.mean(c(intmean, min(spec)/2), c(1 -gaplen/rtlen, gaplen/rtlen))
}

getChromatogram<-function(object, lowMZ, highMZ){
    iter<-max(header(object)[,"seqNum"])
    eic<-sapply(1:iter, function(i, ob){
        pk<-peaks(object, i)
        dat<-sum(pk[which(pk[,1] > lowMZ && pk[,1] < highMZ) ,2])
        return(dat)
    }, ob)
    return(as.numeric(eic))
}
