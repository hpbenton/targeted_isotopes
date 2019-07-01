
foo<-do.call("rbind", pkl)
#### organising the shit stupid matrix thing
save(foo, file="foo_data.Rda")
res<-do.call("rbind", foo[1,1][[1]])

for(i in 1:ncol(foo)){
    for(j in 2:nrow(foo)){
        #if(is.list(foo[j,i])){
        for(k in 1:length(foo[j,i])){
            if(is.matrix(foo[j,i][[k]])){
                idx<-unique(which(foo[j,i][[k]] == "NULL", arr.ind = T)[,1])
                if(length(idx) > 0){
                    pb<-foo[j,i][[k]][-idx,]
                } else{
                    pb<-foo[j,i][[k]]
                    # cat(" . ")
                } ## end of checking for NULL
                
                for(l in 1:ncol(pb)){
                    res<-rbind(res, matrix(pb[,l], ncol=18))
                } ## fucking stupid !!
                
            } else if(is.list(foo[j,i][[k]])){
                res<-rbind(res, do.call("rbind", foo[j,i][[k]]))
                # cat(" - ")
            } ## end of list or matrix 
        } ## end of for k within list block
        #}## end of if list block
        cat(round({i/ncol(foo)}*100), "  ")
    } ## end of row block
    cat("\n")
}## end of column blocks


# {
#     row.names(res)<-NULL
#     ptions(stringsAsFactors = F)
#     write.table(res, file="res.tsv", sep="\t", row.names = F)
#     df<-as.data.frame(res)
#     rm(res)
#     res<-read.table("res.tsv", header=T, sep="\t")
#     res<-res[,-1]
#     colnames(res)<-colnames(df)
# }


res<-res[!(res[,"file"] %in% "NULL"),]
row.names(res)<-NULL
df<-as.data.frame(res)
for(i in 1:ncol(df)){
    df[,i]<-unlist(df[,i])
}
res<-df

### grouping the data by compound and Retention times.
FN<-unique(res[,"file"])
fin.Cnames<-c("mzmen", "ppmMax", "ppmMin", "rtmed", "rtmin", "rtmax", "scanMin", "scanMax", 
              "snMed", "snMin", "snMax", "corMed", "corMin", "corMax", "isotopolog",
              "adduct", "Name", "ID", "mass", "formula", FN) ## then file names for into
fin<-matrix(ncol=c(length(fin.Cnames) ))
uni.res<-unique(cbind(res[,c("Name", "ID")],  "mz"=as.integer(res[,"mz"])) )
colnames(fin)<-fin.Cnames

for(i in 1:nrow(uni.res)){
    idx<-which(res[,"ID"] == uni.res[i,"ID"] & as.integer(res[,"mz"]) == uni.res[i,"mz"])
    if(any(is.na(res[idx,"rtmed"]))){
        ## sort out an odd problem with the peak detector
        rtmed.vec<-as.numeric(res[idx,"rtmed"])
        medRT<-median(rtmed.vec, na.rm=T)
        rtmed.vec[is.na(rtmed.vec)]<-medRT
    }else{
        rtmed.vec<-as.numeric(res[idx,"rtmed"])
    }
    den<-density(rtmed.vec, bw= 20)
    
    while(median(den$y) > 0){ ## make this better
        g.idxMax<-which.max(den$y)
        g.max<-max(den$y)
        g.pr<-xcms:::descendMin(den$y, g.idxMax)
        if(diff(g.pr) < 1){
            break
        }
        rt.pr<-c("rtmin"=den$x[g.pr[1]], "rtmax"=den$x[g.pr[2]])
        
        ## zero out the plots
        den$y[g.pr[1]:g.pr[2]]<-0
        
        ## which peaks go into this group
        res.gr<-which(res[idx,"rtmed"] > rt.pr[1] & res[idx,"rtmed"] <= rt.pr[2])
        sel<-res[idx, ][res.gr,]
        selCh<-sel[,c("adductID", "ID", "Name", "formula", "file")]
        sel<-data.matrix(
            sel[,-match(c("adductID", "ID", "Name", "formula", "file"), colnames(res))]
            ) ## remove the character col from the data
        if(length(unique(selCh[,"file"])) < 2){
            break
        }
        
        ## assing the metaData into features
        mzmen<-weighted.mean(sel[,"mz"], sel[,"maxo"])
        rtmed<-weighted.mean(sel[,"rtmed"], sel[,"maxo"])
        ppmMax<-max(sel[,"deltamz"])
        ppmMin<-min(sel[,"deltamz"])
        scanMin<-min(sel[,"scanMin"])
        scanMax<-max(sel[,"scanMax"])
        snMed<-median(sel[,"sn"])
        snMin<-min(sel[,"sn"])
        snMax<-max(sel[,"sn"])
        corMed<-median(sel[,"corr"])
        corMin<-min(sel[,"corr"])
        corMax<-max(sel[,"corr"])
        isotop<-max(sel[,"isotopolog"])
        adduct<-selCh[1,"adductID"]
        CmpName<-selCh[1,"Name"]
        CmpID<-selCh[1,"ID"]
        mass<-sel[1,"mass"]
        formula<-selCh[1,"formula"]
        
        ## now assign the intensity to the groups
        intenRow<-numeric(length(FN))
        for(j in 1:length(FN)){
            gP.inx<-which(selCh[,"file"] == FN[j])
            into<-max(sel[gP.inx, "into"])
            intb<-max(sel[gP.inx, "intb"])
            # cat(" . ")
            intenRow[j]<-into
        }
        fin<-rbind(fin,
                   c(mzmen, ppmMax, ppmMin,
                     rtmed, rt.pr["rtmin"], rt.pr["rtmax"],
                     scanMin, scanMax, snMed, snMin, snMax,
                     corMed, corMin, corMax, isotop,
                     adduct, CmpName, CmpID, mass, formula, intenRow)
        )
    }
    cat(paste0(round({i/nrow(uni.res)}*100), "\r"))
}
cat("\n\n")


sulphurCompounds<-c("methionine"=149.212665, 
                    "cysteine"=121.159431,
                    "nActylCysteine"=163.196189,
                    "butanoic"=208.279989,
                    "S-methyl-L-methioninate"=163.239283,
                    "N-Formyl-Met-Leu-Phe"=437.554956,
                    "Sinigrin"=397.466269
)

checkIsotopes<-function(raw_data, mass, ppmError=25){
    mzr<-ppmDev(mass, ppmError)
    chr_raw <- chromatogram(raw_data, mz = mzr)
    
    png(file=paste0("CMD-", round(mass), ".png"))
    plot(chr_raw, col="blue")
    dev.off()
    
    S32<-31.97207100
    S34<-33.96786690 
    delta.sulphur<-S34-S32
    isoMass<-mass+delta.sulphur
    
    mzr<-ppmDev(isoMass, ppmError)
    chr_raw <- chromatogram(raw_data, mz = mzr)
    png(file=paste0("CMD-isotope", round(mass), ".png"))
    plot(chr_raw, col="red")
    dev.off()
}



