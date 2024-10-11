#Appendix S3. The R script (ShipleyJVS2013) to perform the simulations reported in Shipley, B. (2013).
#Measuring and interpreting trait-based selection vs. meta-community effects during local community
#assembly. Journal of Vegetation Science (In press).
cats_sim <- function(lambda=0.3, meta.community.u=0.3, random.weight=0.2, 
                     Nperm=500, Nsims=10, 
                     association="no",sd.metacommunity.x=1)
{
  # This script accompanies the paper: Shipley, B. (2013).Measuring and interpreting
  #trait-based selection vs.
  # meta-community effects during local community assembly. Journal of Vegetation
  #Science (in press).
  # Arguments:
  # lambda = strength of selection on (if positive) or against (if negative) the trait
  #in the local community.
  # meta.community.y = strength of selection on an unknown factor in the metacommunity.
  # association (logical) = is there a correlation between the values of the unknown
  #factor in the metacommunity
  # and the trait values?
  # random.weight = strength of the random demographic stochasticity.
  # Nperm, Nsims = number of ramdom permutations of each simulated run and the number of
  #simulated runs.
  # sd.metacommunity.x = standard deviation of the unknown factor in the metacommunity.
  #
  # This is a simulation involving 10 "species" and 1 "trait"
  #
  # Note that this script calls "maxent2" and "maxent.test2" since these replace the
  #original
  # functions in the FD library; when the FD library is updated then these new versions
  #will be
  # called maxent and maxent.test.
  source("src/maxent2.R")
  source("src/maxenttest2_annotated.R")  
  set.seed(19970606)
  # trait holds the trait values (from 1 to 10)
  trait<- c(10, 8, 6, 2, 1, 2, 3, 7, 9, 10)
  trait<-matrix(trait,nrow=1,ncol=10,
                dimnames=list("trait",c("A","B","C","D","E","F","G","H","I","J")))
  model.mean.null.given.uniform<-model.uniform.prior.plus.traits<-
    model.mean.null.given.prior<-
    model.prior.plus.traits<-cor.local.meta<-cor.trait.x<-rep(NA,Nsims)
  delta.R.trait<-delta.R.neutral<-rep(NA,Nsims)
  information.unique.to.local.trait.constraints<-information.unique.to.neutral.prior<-
    joint.information<-biologically.unexplained.information<-rep(NA,Nsims)
  for(i in 1:Nsims){
    # Now run Nsims simulations...
    # q is the prior (info on relative abundance in the meta-community)
    # x is the property of the species determining their meta-community relative
    # abundances; by default it is comes from a uniform random distribution.
    if(association=="no"){
      x<-runif(10,1,10)
      q<-exp(meta.community.u*x)}
    # if you want a positive correlation between meta and local then meta.community.u>0
    # if you want a negative correlation, then meta.community.u<0
    if(association=="yes"){
      x<-round(trait+rnorm(10,0,sd.metacommunity.x),1)
      x.min<-min(x)
      if(x.min<0)x<-x+abs(x.min)
      q<-exp(meta.community.u*x)}
    if(association!="yes" & association!="no")return("error in association argument")
    q<-q/sum(q)
    # the next line generates the observed relative abundances, which are a function of
    # the propagule input (true.metacommunity.ra), the trait selection and the random
    # component representing demographic stochastiticy.
    # rand is a random value (1 to 10) representing demographic stochasticity
    rand<-runif(10,min=1,max=10)
    ra<- q*exp(lambda*trait)*exp(random.weight*rand)
    ra<-ra/sum(ra)
    
    # OWN ADDITION
    ra <- as.matrix(1/(ra[1,1:10])/(sum(1/(ra[1,1:10])))) 
    ra <- c(10, 8, 6, 7, 6, 4, 3, 0, 0, 0 )
    ra <- as.matrix(t(ra/sum(ra)))
    
    q <- c(2, 3, 1, 2, 1, 4, 7, 9,  10, 10)
    q <- as.matrix(t(q/sum(q)))
    # END OWN ADDITION
    
    cor.local.meta[i]<-cor(as.vector(ra),as.vector(q))
    cor.trait.x[i]<-cor(as.vector(trait),as.vector(x))
    ra<-matrix(ra,nrow=1,ncol=10,
               dimnames=list("site",c("A","B","C","D","E","F","G","H","I","J")))
    # CWM is the single community-weighted trait mean...
    CWM<-ra%*%t(trait)
    CWM_meta <- q%*%t(trait)
    #..........................................................................
    # calculating model using the community-weighted trait plus a uniform prior...
    fit1<-maxent2(constr=CWM,states=trait,prior=rep(1/10,10),lambda=TRUE)
    temp1<-maxent.test2(model=fit1,obs=ra,nperm=Nperm,quick=FALSE)
    text(0.001,0.001,"uniform prior")
    # an ESTIMATE (using 1000 permutations) of model bias (model.mean.null.given.uniform)
    model.mean.null.given.uniform[i]<-temp1$mean.KLR2.null
    model.bias<-model.mean.null.given.uniform[i]
    # model fit using the CWM and only the uniform prior...
    model.uniform.prior.plus.traits[i]<-temp1$KLR2.prior.plus.traits
    # if the model using the maximally uninformative prior plus traits is less than the
    # model bias, then correct...
    if(model.uniform.prior.plus.traits[i]<model.bias)model.uniform.prior.plus.traits[i]<-
      model.bias
    #...
    # calculating model using the CWM plus neutral prior
    fit2<-maxent2(constr=CWM,states=trait,prior=q,lambda=TRUE)
    fitted.lambdas<-fit2$lambda
    temp2<-maxent.test2(model=fit2,obs=ra,nperm=Nperm,quick=FALSE)
    text(0.001,0.001,"neutral prior")
    model.mean.null.given.prior[i]<-temp2$mean.KLR2.null
    # if this null, given the neutral prior, is less than the mean given the maximally
    # uniformative prior, then correct...
    if(model.mean.null.given.prior[i]<model.bias)model.mean.null.given.prior[i]<-
      model.bias
    # fit using the CWM and the neutral prior...
    model.prior.plus.traits[i]<-temp2$KLR2.prior.plus.traits
    # if this model, given the neutral prior and the traits, is less than the fit of the model using
    # the maximally uniformative prior and the traits, then correct...
    if(model.prior.plus.traits[i]<model.uniform.prior.plus.traits[i])model.prior.plus.traits[i]<-model.uniform.prior.plus.traits[i]
    #...............................................................................
    just.neutral1<-model.prior.plus.traits[i]-model.uniform.prior.plus.traits[i]
    just.neutral2<-model.mean.null.given.prior[i]-model.mean.null.given.uniform[i]
    just.traits1<-model.uniform.prior.plus.traits[i]-model.mean.null.given.uniform[i]
    just.traits2<-model.prior.plus.traits[i]-model.mean.null.given.prior[i]
    information.unique.to.local.trait.constraints[i]<-just.traits2
    information.unique.to.neutral.prior[i]<-just.neutral1
    delta.R.trait[i]<-just.traits1
    delta.R.neutral[i]<-just.neutral2
    T1<-just.traits1-just.traits2
    T2<-just.neutral2-just.neutral1
    if(abs(T1-T2)>1e-6)return(
      list(message="problem test",T1=T1,T2=T2,traits1=just.traits1,traits2=just.traits2,neutral1=just.neutral1,
           neutral2=just.neutral2,
           diff.traits=abs(just.traits2-just.traits1),diff.neutral=abs(just.neutral1-
                                                                         just.neutral2)))
    joint.information[i]<-T1
    biologically.unexplained.information[i]<-(1-model.prior.plus.traits[i])
  }
  #output.....
  list(
    model.mean.null.given.uniform=round(c(mean(model.mean.null.given.uniform),
                                          sd(model.mean.null.given.uniform)/sqrt(Nsims)),4),
    model.uniform.prior.plus.traits=round(c(mean(model.uniform.prior.plus.traits),
                                            sd(model.uniform.prior.plus.traits)/sqrt(Nsims)),4),
    model.mean.null.given.prior=round(c(mean(model.mean.null.given.prior),sd(model.mean.null.given.prior)/sqrt(Nsims)),4),
    model.prior.plus.traits=round(c(mean(model.prior.plus.traits),sd(model.prior.plus.traits)/sqrt(Nsims)),4),
    delta.R.trait=round(c(mean(delta.R.trait),sd(delta.R.trait)),4),
    delta.R.neutral=round(c(mean(delta.R.neutral),sd(delta.R.neutral)),4),
    information.unique.to.local.trait.constraints=round(c(mean(information.unique.to.local.trait.constraints),
                                                          sd(information.unique.to.local.trait.constraints)/sqrt(Nsims)),4),
    information.unique.to.neutral.prior=round(c(mean(information.unique.to.neutral.prior),
                                                sd(information.unique.to.neutral.prior)/sqrt(Nsims)),4),
    joint.information=round(c(mean(joint.information),sd(joint.information)/sqrt(Nsims)),4
    ),
    biologically.unexplained.information=round(c(mean(biologically.unexplained.information
    ),sd(biologically.unexplained.information)/sqrt(Nsims)),4),
    cor.local.metacommunity.ra=round(c(mean(cor.local.meta),sd(cor.local.meta)/sqrt(Nsims)
    ),4),
    cor.trait.x=round(c(mean(cor.trait.x),sd(cor.trait.x)/sqrt(Nsims)),4),
    observed.ra=ra,observed.metacommunity.ra=q,
    trait = trait
  )
  
}
