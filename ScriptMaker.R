scriptmaker<-function(Title,
                      Matrix, # expects a csv file, with row names in first column and column names in first row.
                      Species, # expects a csv file as per directions. 
                      Years=1,
                      Timesteps, # Number of population records
                      TimeUnit, # Unit of time between each record
                      Data, # list of tab delimiated .txt files, with species in the header, with the population dynamic data. Not needed if using Fake Data
                      GlobalParameters, # .csv of global parameter means and bounds
                      AssimEf=2, # Assimilation efficency
                      MakeFake = F, # T/F). Incompatible with multiyear at present.
                      SeperateParams = T, # If set to false, consumers and producers share a single aT, aJ and aR, which will need to be added manualy to the code outut if you want to make fake data
                      AbsoluteError=1.0e-12, # Allowed error of the ODEsolver stepper function
                      RelativeError=1.0e-6, # Allowed error of the ODEsolver stepper function
                      Pset=-1, # Set global control of dynamic parameters  0=randomly choose start between bounds, -1=fix start at what is in table, 1= all parameters fixed (can then go and manually loosen a few)
                      STARTTYPE=1, # 1=  will start all the start points at the first data value, bounds set by table (unlikely to be valid for multiyear!)
                                   # 2= will start all the start points set in the parameter file (unlikely to be valid for multiyear!)
                                   # 3=  will start all the start points at the first data value, bounds given % leeway either side. 
                      StartBounds=5, # % leeway for bounds when using start type = 3. 
                      StartFixed=1, # as filzbach
                      RejectNegs = F, #  Set to 'True' to include a part of the liklihood that will cause immediate rejection of any parameters that result in negative populations or unrealitically large populations.  # This is most commonly caused by solver steps being too small, so be careful! 
                      MaxForReject= 1000, # Level if above the  population goes above the likilhood is stopped
                      MCMC=c(1000,1000,0,0,1) # List of number of iterations for each of Filzbach stages (Burn-in, Bayesian Sampling, Maximim Likelihood adjustments, Maximum Likilihood search, Number of Chains)
  ){
  
#########################
#### Useful Quantities  
  M<-Matrix
  print(M)
  all<- rownames(M)
  TotS <- nrow(Species)
print(TotS)
  preds<-colnames(M)
  PnumAdj <- length(all)-length(preds)-1  # Number used to output index used in the c++ code in some of the ODE code
  print(Species) # Print to screen to CHECK THAT THERE ARE NOT ANY VALUES TOO FAR OUT!
  ## Get other parameters that need to be defined if generating fake data, also start point for parameter search
  fR <-subset(GlobalParameters, Parameter == 'fR') # fraction of maximum producer growth rate realised
  K  <-subset(GlobalParameters, Parameter == 'K') # total carrying capacity of phytoplanton
  sdZ<-subset(GlobalParameters, Parameter == 'sdZ') # Measurement error for zooplankton on a lognormal scale
  sdP<-subset(GlobalParameters, Parameter == 'sdP') #  Measurement error for phytoplankton on a lognormal scale
  fJ <-subset(GlobalParameters, Parameter == 'fJ') # fraction of realised maximum ingestion rate
  d  <-subset(GlobalParameters, Parameter == 'd') # prey defense parameter 
  W  <-subset(GlobalParameters, Parameter == 'W') # total prey biomass at which the functional response is 0.5 
  q  <-subset(GlobalParameters, Parameter == 'q') # Type II vs Type III
 
###
Headers <- paste('#include "preamble.h" \n#ifdef ',Title,
                     ' \n#include <stdlib.h> \n#include <math.h> \n#include <string.h>
#include "filzbach.h"
#include <iostream> 
#include <boost/array.hpp> 
#include <boost/numeric/odeint.hpp>
using namespace std;
using namespace boost::numeric::odeint;
void pause(){PAUSE}  
void fake_data(); 
void read_data(); 
void setup_parameters();
void fit_model();
void final_output();\n')
    
    ParamDefs<- c()
    aTparams <- c()
    aJparams<- c()
    aRparams<- c()
    starts<-c()
    InverseMasses<-c()
    

if(SeperateParams==T){
    if(MakeFake==T){ 
      for(i in 1:TotS){
        
        if(Species$Type[i] == 'Producer'){
          aRparams<-paste0(aRparams,'double aR',all[i],' = ',Species$aR[i],';\n ')
        }
        if(Species$Type[i] == 'Consumer'){
          aTparams<-paste0(aTparams,'double aT',all[i],' = ',Species$aT[i],';\n ')
          aJparams<-paste0(aJparams,'double aJ',all[i],' = ',Species$aJ[i],';\n ')
        }
        starts<-paste0(starts,'double start', all[i], ' = ', Species$Start[i],'; \n')
        InverseMasses<-paste0(InverseMasses,'const double im',all[i],'= ', signif(Species$Masses[i]^(-0.25),3),'; \n')
      } 
    }
    if(MakeFake==F){ 
      for(i in 1:TotS){
        
        if(Species$Type[i] == 'Producer'){
          aRparams<-paste0(aRparams,'double aR',all[i],';\n ')
        }
        if(Species$Type[i] == 'Consumer'){
          aTparams<-paste0(aTparams,'double aT',all[i],';\n ')
          aJparams<-paste0(aJparams,'double aJ',all[i],';\n ')
        }
        if(Years==1){ starts<-paste0(starts,'double start', all[i],'; \n')}
        if(Years>1){for(y in 1:Years){starts<-paste0(starts,'double start',y, all[i],'; \n')}
        }
        InverseMasses<-paste0(InverseMasses,'const double im',all[i],'= ', signif(Species$Masses[i]^(-0.25),3),'; \n')
      } 
    } 
}

if(SeperateParams==F){
  if(MakeFake==T){ 
        aRparams<-paste0(aRparams,'double aR = ***;\n ')
        aTparams<-paste0(aTparams,'double aT = ***;\n ')
        aJparams<-paste0(aJparams,'double aJ = ***;\n ')
      
        
  for(i in 1:TotS){
      starts<-paste0(starts,'double start', all[i], ' = ', Species$Start[i],'; \n')
      InverseMasses<-paste0(InverseMasses,'const double im',all[i],'= ', signif(Species$Masses[i]^(-0.25),3),'; \n')
    } 
  } 
if(MakeFake==F){ 
           aRparams<-paste0(aRparams,'double aR;\n ')
           aTparams<-paste0(aTparams,'double aT;\n ')
           aJparams<-paste0(aJparams,'double aJ;\n ')
  
      for(i in 1:TotS){
      if(Years==1){ starts<-paste0(starts,'double start', all[i],'; \n')}
      if(Years>1){for(y in 1:Years){starts<-paste0(starts,'double start',y, all[i],'; \n')}
      }
      InverseMasses<-paste0(InverseMasses,'const double im',all[i],'= ', signif(Species$Masses[i]^(-0.25),3),'; \n')
    } 
  } 
}

    Qs<-c()
    for(i in 0:(TotS-1)){Qs<-paste0(Qs,'double Q',i,'; double W',i,';\n')}
    
    Others<- paste(
      'double fJ = ',fJ$TruthForFakeData,';
      double fR = ',fR$TruthForFakeData,';
      double d = ',d$TruthForFakeData,';
      double K = ',K$TruthForFakeData,';  
      double sdZ = ',sdZ$TruthForFakeData,'; 
      double sdP = ',sdP$TruthForFakeData,';
      double W = ',W$TruthForFakeData,';
      double q = ',q$TruthForFakeData,';
      double PhytoTotal; \n',Qs,'\n'
    )
    
    ParamDefs<-paste0(aTparams,aJparams, aRparams, starts,Others,
                      '// Note that these are mass ^ -0.25 \n',InverseMasses,
                      '\n 
                      
                      typedef boost::array< double , ',TotS,' > state_type; // Create internal array the size of what needs to be tracked
                      typedef runge_kutta_cash_karp54< state_type > error_stepper_type; // Type of stepper used see http://headmyshoulder.github.io/odeint-v2/index.html docs for details
                      typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
                      double abs_err = ',AbsoluteError,' , rel_err = ',RelativeError,' , a_x = 1.0 , a_dxdt = 1.0;
                      controlled_stepper_type controlled_stepper(
                      default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );
                      int timerow = 0; \n')
    
    ################################
    ## Call ODE Function
    ###############################
    
    Type<- Species$Type
    all<- rownames(Matrix)
    preds<-colnames(Matrix)
    PnumAdj <- length(all)-length(preds)-1  # Number used to output index used in the c++ code in some of the ODE code
   
    TotalCode<-c()  
    DefinitionCode<-c()
    PhytoTotal<-'PhytoTotal = 0' # Gets addes to a little later
    for(species in 1:length(all)){DefinitionCode <- paste0(DefinitionCode,'W',species-1,'=x[',species-1,']/W;', 'Q',species-1,'=pow(W',species-1,',q);')} ## Up here calculate all the necessary biomass powers
    

###
if(SeperateParams==F){SPall<-rep('',length(all) )} # insert a space rather than anything useful if seperate parameters are not needed at certain places
if(SeperateParams==T){SPall<-all} # Otherwise carry on as normal
###


    for(i in 1:length(all)) # Cycle through all the species detailed in rows
    {
      IntroCode<-c(); ProdCode<-c();RespCode<-c();LossCode<-c();GainCode<-c(); ## Reset all segments for the new species
      
      ## First Create Title and LHS
      IntroCode<- paste0('/*  ', all[i],' ')
      if(nchar(IntroCode) >10){print('Your names are pretty long, consider shortening them')
      } else{while(nchar(IntroCode) < 11){IntroCode<-paste0(IntroCode," ")}}
      IntroCode<- paste0(IntroCode,' */ dxdt[',i-1,'] =') 
      
      ## Second see if it is a producer and intrinsically gains biomass or a consumer and loses mass through respiration
      if(Type[i]=='Producer'){
        print(paste(all[i], 'is a producer'))
        ProdCode <- paste0('fR*aR', SPall[i],'*im',all[i],'*x[',i-1,']*(1 -(PhytoTotal/K)) \n                          ')
        PhytoTotal<- paste0(PhytoTotal, '+x[',i-1,']')
      }
      if(Type[i]=='Consumer'){
        print(paste(all[i], 'is a consumer'))
        RespCode <- paste0('-(aT',SPall[i],'*im',all[i],'*x[',i-1,']) // Loss from Respiration \n                          ')
      }
      
      #Third, if instead it is a consumer, see what it gains from consuming others
      ## Is it a consumer?? See whether it is worth going through the loop by seeing if it has a column)
      if(i > length(all)-length(preds) )
      {
        column<- i-(length(all)-length(preds))
        print(paste(all[i], "has", sum(M[,column]), "prey"))
        numprey<- sum(M[,column]) # the number of prey that the predator has
        Pabr<-all[i] # Get the standard abbreviation of the predator in question
        
        GainCode<-paste0(GainCode,'fJ*aJ', SPall[i],'*im',Pabr,'*x[', i-1,'] *(0' ) 
        block<-c()
        for(prey in 1:length(all)) # Cycle through all of the potential prey to calculate the functional responses of each prey
               {if(M[prey,column]==1) {block<-paste0(block,"+Q", prey-1)}}
        GainCode<- paste0(GainCode , block , ')/(1 + d*x[' , column+PnumAdj , '] +(0' , block ,'))//  Gain from predation \n                          ')
      }else{print(paste(all[i], "does not have a predation column"))}
      
      ## Fourth see if it has any predators that cause it to lose biomass
      if (sum(M[i,])>0)
         { print(paste(all[i], "has", sum(M[i,]), "predator(s)"))
           numpreds<- sum(M[i,]) # the number of predators that the prey has
             for (pred in 1:length(preds)) # Cycle through all of the potential predators 
                     {
                    if(M[i,pred]==1){
                             Pabr<-preds[pred] # Get the standard abbreviation of the predator in question
                             # First part # ingestion rate * biomass of predators
                             if(SeperateParams==T){ LossCode<-paste0(LossCode,' -(fJ*aJ', Pabr,'*im',Pabr,'*',AssimEf,'*x[', pred+PnumAdj,']' )}  
                             if(SeperateParams==F){ LossCode<-paste0(LossCode,' -(fJ*aJ*im',Pabr,'*',AssimEf,'*x[', pred+PnumAdj,']' )}  
                             # Functional Response
                             top<-paste0(' * Q', i-1 )
                             bottom<- paste0('1+d*x[',pred+PnumAdj,']+(')
                              for(z in 1:length(all)) # Cycle through all the prey of the predator: z = index of all the prey species that the predator eats
                                       {if(M[z,pred]==1) {bottom <- paste0(bottom, "Q",z-1,'+')}}
                          LossCode<-paste0(LossCode,paste0(top,'/(',bottom,'0))) // Impact of ',Pabr,'\n                          ')) #Adding the 0 just means there doesn't have to be any mucking around with counters to switch off inserting the +'s
                                     }
                  }
         }else{ print(paste(all[i], "has no predators"))}
      TotalCode<-paste(TotalCode, IntroCode,ProdCode, GainCode,RespCode,LossCode, ';\n')
    }  
    
    PhytoTotal<- paste0(PhytoTotal, '; \n')
    ODEs<-paste0('\n',
                 ' void ODE( const state_type &x , state_type &dxdt , double t )
{',DefinitionCode,PhytoTotal,TotalCode, '\n } \n ')
    
    ##########################
    
    if(MakeFake == T)  # write_fakeHy()
    {
      rec<-c();for(i in 1:TotS){rec<-paste0(rec,'double rec',all[i],'; ')}
      write_fakeHy<- paste0(rec,'\n
                            void write_fakeHy( const state_type &x , const double t )
{
                            ')
      recdraw<-c()
      for(i in 1:TotS){
        if(Species$Type[i]=='Producer')
        {recdraw<-paste0(recdraw,'double rec',all[i],' = lognormal_draw(x[',i-1,'],sdP);\n')} 
        if(Species$Type[i]=='Consumer')
        {recdraw<-paste0(recdraw,'double rec',all[i],' = lognormal_draw(x[',i-1,'],sdZ);\n')} 
      }
      write_fakeHy<-paste0(write_fakeHy,'\n',recdraw)
      
      cout<-'cout << t'
      for(i in 1:TotS){cout<-paste0(cout,"<<'\\t'<< x[",i-1,']')} # double backslash escapes the escape....
      write_fakeHy <- paste0(write_fakeHy,cout,'<< endl; \n')
      
      for(i in 1:TotS){write_fakeHy<-paste0(write_fakeHy,'table_writevalue("fakedata", "',all[i],'", timerow,rec',all[i],'); \n' )  }
      
      write_fakeHy<-paste0(write_fakeHy, 'timerow++; 
} \n')
}
    
    if (MakeFake == T) # fake_data()
    {
      fakedata<- 'void fake_data()
{table_create("fakedata"); \n'
      
      for(i in 1:TotS){fakedata<-paste0(fakedata,'table_addcolumn("fakedata","',all[i],'");\n')}
      
      fakedata<- paste0(fakedata,'state_type x = {')
      
      for(i in 1:(TotS-1)){fakedata<-paste0(fakedata,'start',all[i],',')}
      fakedata<-paste0(fakedata,'start',all[TotS],'}; \n')
      
      fakedata<-paste0(fakedata,'integrate_const(controlled_stepper,ODE,x,0.0,', Timesteps*TimeUnit ,',', TimeUnit , ',write_fakeHy ); \n')
      fakedata<-paste0(fakedata,'table_output("fakedata","./workspace/',Title,'fake.txt"); \n return;} \n')
      
}
    
    #############
    #read_data()
    ###############
    readdata<- 'void read_data()
{'
    if(MakeFake==T)
    {readdata<-paste0(readdata,'table_read("mydata","./workspace/',Title,'fake.txt",',TotS,');\n' )}
    if(MakeFake==F && Years==1)
    { 
      readdata<-paste0(readdata,'table_read("mydata","./workspace/',Data,'",',TotS,');\n' )
      for(i in 1:TotS)
      {readdata<-paste0(readdata, 'table_addcolumn("mydata", "found', all[i],'"); \n')  }}
    
    if(MakeFake==F && Years>1)
    { for(y in  1:Years){readdata<-paste0(readdata,'table_read("mydata',y,'","./workspace/',Data[y],'",',TotS,');\n' )
                         for(i in 1:TotS)
                         {readdata<-paste0(readdata, 'table_addcolumn("mydata',y,'", "found', all[i],'"); \n')  }}}
    
    readdata<-paste0(readdata,'return;} \n \n')
    #######################
    # setup_parameters()  #
    #######################
    SP<- 'void setup_parameters()
{   /* min - max - initial - norm/lnorm -  fixed? - display?  */ \n'
    

if(SeperateParams==T){
    for (i in 1:TotS){if(Species$Type[i]=='Producer'){SP<- paste0(SP,'parameter_create("aR',all[i],'",',Species$aRLower[i],',',Species$aRUpper[i],',',Species$aR[i],',0,',Pset,',1); \n')} }
    for (i in 1:TotS){if(Species$Type[i]=='Consumer'){SP<- paste0(SP,'parameter_create("aJ',all[i],'",',Species$aJLower[i],',',Species$aJUpper[i],',',Species$aJ[i],',0,',Pset,',1); \n')} }
    for (i in 1:TotS){if(Species$Type[i]=='Consumer'){SP<- paste0(SP,'parameter_create("aT',all[i],'",',Species$aTLower[i],',',Species$aTUpper[i],',',Species$aT[i],',0,',Pset,',1); \n')} }
}
if(SeperateParams==F){
SP<- paste0(SP,'parameter_create("aR",Min,Max,Start,0,',Pset,',1); \n','parameter_create("aJ",Min,Max,Start,0,',Pset,',1); \n','parameter_create("aT",Min,Max,Start,0,',Pset,',1); \n')
}

    if(Years==1){
      if(STARTTYPE==1){
        st<-read.delim(Data)[1,]
        for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',all[i],'",',Species$StartLower[i],',',Species$StartUpper[i],',',st[i],',0,', StartFixed,',1);\n')}
      }
      if(STARTTYPE==2){
        for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',all[i],'",',Species$StartLower[i],',',Species$StartUpper[i],',',Species$Start[i],',0,', StartFixed,',1);\n')}
      }    
      if(STARTTYPE==3){
        st<-read.delim(Data)[1,]
        for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',all[i],'",',st[i]*(1-(StartBounds*0.01)),',',st[i]*(1+(StartBounds*0.01)),',',st[i],',0,', StartFixed,',1);\n')  }
      }      
    }else{
      for(y in 1:Years){
        st<-read.delim(Data[y])[1,]       
        
        if(STARTTYPE==1){
          for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',y,all[i],'",',Species$StartLower[i],',',Species$StartUpper[i],',',st[i],',0,', StartFixed,',1);\n')}
        }
        if(STARTTYPE==2){
          for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',y,all[i],'",',Species$StartLower[i],',',Species$StartUpper[i],',',Species$Start[i],',0,', StartFixed,',1);\n')}
        }    
        if(STARTTYPE==3){
          for (i in 1:TotS){SP<- paste0(SP,'parameter_create("start',y,all[i],'",',st[i]*(1-(StartBounds*0.01)),',',st[i]*(1+(StartBounds*0.01)),',',st[i],',0,', StartFixed,',1);\n')}
        } 
      }
    }
    
    ####################################################################################################################################
    SP<- paste0(SP,'parameter_create("fR','",',fR$LowerBound,',',fR$UpperBound,',',fR$Start,',0,',1,',1);\n')
    SP<- paste0(SP,'parameter_create("K','",',K$LowerBound,',',K$UpperBound,',',K$Start,',0,',Pset,',1);\n')
    SP<- paste0(SP,'parameter_create("sdZ','",',sdZ$LowerBound,',',sdZ$UpperBound,',',sdZ$Start,',0,',Pset,',1);\n')
    SP<- paste0(SP,'parameter_create("sdP','",',sdP$LowerBound,',',sdP$UpperBound,',',sdP$Start,',0,',Pset,',1);\n')
    SP<- paste0(SP,'parameter_create("fJ','",',fJ$LowerBound,',',fJ$UpperBound,',',fJ$Start,',0,',1,',1);\n')
    SP<- paste0(SP,'parameter_create("d','",',d$LowerBound,',',d$UpperBound,',',d$Start,',0,',Pset,',1);\n')
    SP<- paste0(SP,'parameter_create("W','",',W$LowerBound,',',W$UpperBound,',',W$Start,',0,',Pset,',1);\n')
    SP<- paste0(SP,'parameter_create("q','",',q$LowerBound,',',q$UpperBound,',',q$Start,',0,',Pset,',1);\n')
    
    ### TODO Add priors in here
    
    SP<- paste0(SP, 'parameter_showall();
                PAUSE
                return;
}  \n ')
##########################
    # Write_Hy()'s
    #################
    
    
    if(Years==1){
      write_Hy<- paste0('void write_Hy( const state_type &x , const double t )
{ \n')
      if(RejectNegs==T){for(i in 1:TotS){write_Hy<-paste0(write_Hy,'if(x[',i-1,']<0.0|| x[',i-1,']>',MaxForReject,'){set_metr_ltotnew(-99999.0);force_reject();return;};\n')}}
      
      for(i in 1:TotS){write_Hy<-paste0(write_Hy,'table_writevalue_multichain("mydata", "found',all[i],'", timerow,x[',i-1,']); \n' )  }
      write_Hy<-paste0(write_Hy, 'timerow++;
} \n')
  
}else{write_Hy <-c()
      for(y in 1:Years){  
        write_Hy<- paste0(write_Hy,'void write_Hy',y,'( const state_type &x , const double t )
{\n')
        if(RejectNegs==T){for(i in 1:TotS){write_Hy<-paste0(write_Hy,'if(x[',i-1,']<0.0|| x[',i-1,']>',MaxForReject,'){set_metr_ltotnew(-99999.0);force_reject();return;};\n')}}
        
        for(i in 1:TotS){write_Hy<-paste0(write_Hy,'table_writevalue_multichain("mydata',y,'", "found',all[i],'", timerow,x[',i-1,']); \n' )  }
        write_Hy<-paste0(write_Hy, 'timerow++; 
} \n') 
} 
}
    
    #################
    # Likilihood()
    #################
    
    Like1<-c()
    
    for(i in 1:TotS)
    {Like1<-paste0(Like1, 'double prob', all[i], '; double real', all[i],'; double found', all[i],'; \n')}
    
    Like1 <- paste0(Like1, 'void likelihood()
{
set_metr_ltotnew(0.0); /* set sum over log-Likelihood to zero */
set_metr_number_ok(0);
int numdata;
timerow = 0;\n' )
    

if(SeperateParams==T){

    for (i in 1:TotS)
    {if(Species$Type[i]=='Producer'){Like1<- paste0(Like1,'aR',all[i],'=cv("aR',all[i],'");\n')} }
    for (i in 1:TotS)
    {if(Species$Type[i]=='Consumer'){Like1<- paste0(Like1,'aT',all[i],'=cv("aT',all[i],'");\n')} }
    for (i in 1:TotS)
    {if(Species$Type[i]=='Consumer'){Like1<- paste0(Like1,'aJ',all[i],'=cv("aJ',all[i],'");\n')} }  
}
if(SeperateParams==F){Like1<- paste0(Like1,'aR=cv("aR");\naJ=cv("aJ");\naT=cv("aT");\n')}
    
    
    if(Years==1){
      for (i in 1:TotS){Like1<- paste0(Like1,'start',all[i],'=cv("start',all[i],'");\n')} 
    }else{
      for(y in 1:Years){for (i in 1:TotS){Like1<- paste0(Like1,'start',y,all[i],'=cv("start',y,all[i],'");\n')}
      }}
    
    Like1<- paste0(Like1,'fR','=cv("fR','");\n')
    Like1<- paste0(Like1,'K','=cv("K','");\n')
    Like1<- paste0(Like1,'sdZ','=cv("sdZ','");\n')
    Like1<- paste0(Like1,'sdP','=cv("sdP','");\n') 
    Like1<- paste0(Like1,'fJ','=cv("fJ','");\n') 
    Like1<- paste0(Like1,'d','=cv("d','");\n') 
    Like1<- paste0(Like1,'W','=cv("W','");\n')
    Like1<- paste0(Like1,'q','=cv("q','");\n')
    
    Like2<-c()
    if(Years==1){
      Like2<-('  timerow = 0;
              state_type x = {')
      
      for(i in 1:(TotS-1)){Like2<-paste0(Like2,'start',all[i],',')}
      Like2<-paste0(Like2,'start',all[TotS],'}; \n')
      Like2<-paste0(Like2,'integrate_const(controlled_stepper,ODE,x,0.0,', Timesteps*TimeUnit ,',', TimeUnit , ',write_Hy ); \n')
      
      
      Like2<- paste0(Like2, 'numdata = ',Timesteps,';
                     for(int ii = 0; ii < numdata; ii++)
{ ')
      
      for(i in 1:TotS){Like2<-paste0(Like2, 'real',all[i],' = table_getvalue("mydata", "',all[i],'",ii);\n')}
      for(i in 1:TotS){Like2<-paste0(Like2, 'found',all[i],' = table_getvalue("mydata", "found',all[i],'",ii);\n')}
      for(i in 1:TotS){if(Species$Type[i]=='Producer')
      {Like2<-paste0(Like2, 'inc_metr_ltotnew(log(prob',all[i],'=lognormal_density(real',all[i],', found',all[i],',sdP)));\n')}}
      for(i in 1:TotS){if(Species$Type[i]=='Consumer') 
      {Like2<-paste0(Like2, 'inc_metr_ltotnew(log(prob',all[i],'=lognormal_density(real',all[i],', found',all[i],',sdZ)));\n')}}
      Like2<- paste0(Like2,'inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on \n }')
}else{for(y in 1:Years){
  Like2<-paste0(Like2,' \n  timerow = 0;\n 
                state_type x',y,' = {')
  
  for(i in 1:(TotS-1)){Like2<-paste0(Like2,'start',y,all[i],',')}
  Like2<-paste0(Like2,'start',y,all[TotS],'}; \n')
  Like2<-paste0(Like2,'integrate_const(controlled_stepper,ODE,x',y,',0.0,', Timesteps*TimeUnit ,',', TimeUnit , ',write_Hy',y,' ); \n')
  
  Like2<- paste0(Like2, 'numdata = ',Timesteps,';
                 for(int ii = 0; ii < numdata; ii++)
{ ')
  
  for(i in 1:TotS){Like2<-paste0(Like2, 'real',all[i],' = table_getvalue("mydata',y,'", "',all[i],'",ii);\n')}
  for(i in 1:TotS){Like2<-paste0(Like2, 'found',all[i],' = table_getvalue("mydata',y,'", "found',all[i],'",ii);\n')}
  for(i in 1:TotS){if(Species$Type[i]=='Producer')
  {Like2<-paste0(Like2, 'inc_metr_ltotnew(log(prob',all[i],'=lognormal_density(real',all[i],', found',all[i],',sdP)));\n')}}
  for(i in 1:TotS){if(Species$Type[i]=='Consumer') 
  {Like2<-paste0(Like2, 'inc_metr_ltotnew(log(prob',all[i],'=lognormal_density(real',all[i],', found',all[i],',sdZ)));\n')}}
  Like2<- paste0(Like2,'inc_metr_number_ok(1); // adds one to the sample size, used in BIC later on \n }')
}
                }
    
    Like2<- paste0(Like2,'
                   return;
} \n')

    ################
    # final_output() and main()
    ################
    
    if(Years==1){
      end<- paste0('void final_output()
{
                   /* create link to file */
                   char fname[100];
                   /* create file name for output -- outp is the path */
                   get_filzbach_path(fname, 100);
                   /* now your bit to add to the path */
                   strcat_s(fname,"_my_output.txt");
                   /* run the most likely set of parameters */
                   params_set_to_posterior_mean();
                   likelihood();
                   /* print internal table out to file with name taken from above */
                   table_output("mydata",fname);
                   return;
}
                   
                   /* ************************************************* */
                   /* The control function that calls everything else.  */
                   /* ************************************************* */
                   
                   int main()
{
                   atexit(pause);
                   // set the likelihood function pointer
                   pfn_likelihood = &likelihood;
                   initialize_filzbach();
                   name_analysis("', Title, '"); \n')
}else{
  ZZZZ<-c()  
  for(y in 1:Years){ZZZZ<-paste0(ZZZZ, 'char fname',y,'[100];
                                 get_filzbach_path(fname',y,', 100);
                                 strcat_s(fname',y,',"_my_outputYEAR',y,'.txt");
                                 table_output("mydata',y,'",fname',y,');' )}
  end<- paste0('void final_output()
{ /* run the most likely set of parameters */
               params_set_to_posterior_mean();
               likelihood();',ZZZZ,' \nreturn;
}
               
               /* ************************************************* */
               /* The control function that calls everything else.  */
               /* ************************************************* */
               
               int main()
{
               atexit(pause);
               // set the likelihood function pointer
               pfn_likelihood = &likelihood;
               initialize_filzbach();
               name_analysis("', Title, '"); \n')  
}

if(MakeFake==T){end<-paste0(end,'            fake_data(); \n')}

end<- paste0(end,
'            read_data();
             setup_parameters();
             set_chains(',MCMC[5],');
             runmcmc(',MCMC[1],',',MCMC[2],',',MCMC[3],',',MCMC[4],'); // Burn-in, Bayes, MLE adaption, MLE search
             final_output();
}
             #endif
             ')
############


if(MakeFake==T){GrandOutput<-paste0(Headers,ParamDefs,ODEs,write_fakeHy,fakedata,readdata,SP,write_Hy, Like1, Like2,end)}
if(MakeFake==F){GrandOutput<-paste0(Headers,ParamDefs,ODEs,readdata,SP,write_Hy, Like1, Like2,end)}


################
dir.create('Output Scripts', showWarnings = FALSE)
cat(GrandOutput,file=paste0('Output Scripts/',Title,"CppODEs.txt"))

cat('************************
All done, now check it!
************************')

}