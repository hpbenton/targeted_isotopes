#Quick Run 

##Parameters 

Some needed parameters are 
adductSelection="standard"
minScanTime=3; ppmDevN=7.5; snthresh=2; maxoFilt=1000; plot.prog=F; atom<-"C"

You will also need a selection of compounds to target. This needs to have the following header names :

ID, name, mass, formula
## Code to run
```R
 out<-runSulphurFinder(atom="C", minScanTime=3, ppmDevN=20, snthresh=2, maxoFilt=1000, adductSelection=c("all", "standard") )
 
 ```
 
## ToDO
Currently the grouping code is not in a function
Change code into objects
Remove fixed quick run variables
Have code to check number of cores on the machine.
 
## to be cleaned up later 
