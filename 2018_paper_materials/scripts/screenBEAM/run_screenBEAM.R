library(ScreenBEAM)

args = commandArgs(trailingOnly=TRUE)

if( length(args) != 4 ){
    print(args)
    stop("Usage: run_screenBEAM.R countfile cases ctrls outfile")
} 
print(args[1])
print(strsplit(args[3],","))
print(strsplit(args[2],","))

###NGS data
r<-ScreenBEAM(

  ###input format
  input.file=args[1] # system.file("extdata", args[1], package = "ScreenBEAM")#tab-separted file
  ,
  control.samples=unlist(strsplit(args[3],",")) #column names of control samples
  ,
  case.samples=unlist(strsplit(args[2],",")) #column names of case/treated samples
  ,
  control.groupname='CTRL'#name your control group
  ,
  case.groupname='CASE'#name your case group
  ,

  ###data pre-processing
  data.type='NGS'#data type
  ,
  do.normalization=TRUE
  ,
  filterLowCount=TRUE
  ,
  filterBy = 'control'
  ,
  count.cutoff=4
  ,

  ###Bayesian computing
  nitt=1500,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=500#number of burnin in MCMC sampling, 5000 is default

  )

write.csv(r,file=file.path(args[4]),row.names=FALSE,na='')
