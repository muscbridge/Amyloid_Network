PROC IMPORT OUT= WORK.NEW2 
            DATAFILE= "S:\CU_Projects\Paul\Nolin, Sara\Data\data_long.xlsx" 
            DBMS=EXCEL REPLACE;
     RANGE="data_long$"; 
     GETNAMES=YES;
     MIXED=NO;
     SCANTEXT=YES;
     USEDATE=YES;
     SCANTIME=YES;
RUN;

proc sort data=new2;
by studyid network_type network;
run;

proc mixed data=new2;
class studyid  ;
model amyloid =  / solution cl;
random int  / subject=studyid vcorr v;
ods output covparms=c1;
run;
data c2; set c1;
sqrt_estimate = sqrt(estimate);
run;

** the model above yields an estimate for between person variance of 0.005406,
corresponding to a SD of 0.0735229378, which will be used in the effect size calculations.;


** model 1;
proc mixed data=new2;
class studyid  network_type;
model amyloid = age network_type / solution cl outp=out2;
repeated  / subject=studyid type=un rcorr;
lsmeans network_type / diff cl ;
run;

** model 2;
proc mixed data=new2;
class studyid  network_type ;
model amyloid = age | network_type / solution cl;
repeated  / subject=studyid type=un rcorr;
run;

** model 3;
proc mixed data=new2;
class studyid network;
model amyloid = age | network / solution cl;
repeated  / subject=studyid type=un rcorr;
lsmeans network / diff cl adj=tukey;
run;

** model 4;
proc mixed data=new2;
class studyid network_type pet_read;
model amyloid = age network_type | pet_read / solution cl;
repeated  / subject=studyid type=un rcorr;
run;

** model 5;
proc mixed data=new2;
class studyid network pet_read;
model amyloid = age pet_read|network / solution cl;
repeated  / subject=studyid type=un rcorr;
lsmeans pet_read*network / cl diff adj=tukey;
run;
