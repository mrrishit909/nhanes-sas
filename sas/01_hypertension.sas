/*  Hypertension in US adults, NHANES August 2021 - August 2023.
    Part 1 reproduces CDC's NCHS Data Brief No. 511 with the full survey design (strata, PSUs, exam weights).
    Part 2 asks who is being missed: survey logistic models for being unaware of hypertension, and for being
    uncontrolled while on treatment.
    Ran in SAS Enterprise Guide (SAS 9.4) on USF's virtual desktop; the log is sas/01_hypertension.log. Lines starting
    CSV| in the log are the results; the repository's check.py recomputes them in Python.                        */

options nodate nonumber validvarname=v7;
%put NOTE: SAS &sysvlong on &sysscp &sysscpl;

/* 1. Read the five CDC transport files straight from the web */
%macro get(f);
  filename f_&f temp;
  proc http url="https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/&f..xpt" out=f_&f; run;
  libname x_&f xport "%sysfunc(pathname(f_&f))";
  proc sort data=x_&f..&f out=&f; by seqn; run;
%mend;
%get(DEMO_L) %get(BPXO_L) %get(BPQ_L) %get(BMX_L) %get(HIQ_L)

/* 2. One analysis file. Everyone stays in it: subgroups are domains, never subsets, so standard errors are right. */
data nh;
  merge demo_l(in=in_demo) bpxo_l bpq_l bmx_l hiq_l;
  by seqn;
  if in_demo;
  sbp = mean(of bpxosy1-bpxosy3);                     /* average of up to three oscillometric readings */
  dbp = mean(of bpxodi1-bpxodi3);
  adult = (ridageyr >= 18 and sbp > . and dbp > . and wtmec2yr > 0 and ridexprg ne 1);   /* as the data brief */
  meds = (bpq150 = 1);
  htn = (sbp >= 130 or dbp >= 80 or meds);            /* 2017 ACC/AHA threshold, as the data brief */
  if bpq020 in (1, 2) then aware = (bpq020 = 1);      /* don't know / refused stay missing */
  controlled = (sbp < 130 and dbp < 80);
  hyp = (adult and htn);
  treated = (hyp and meds);
  if hyp and aware ne . then unaware = 1 - aware;
  if treated then uncontrolled = 1 - controlled;
  length agegrp $5 sex $5 race $11 income $7 insured $7 bmi $7;
  if 18 <= ridageyr < 40 then agegrp = '18-39'; else if 40 <= ridageyr < 60 then agegrp = '40-59'; else agegrp = '60+';
  sex = ifc(riagendr = 1, 'Men', 'Women');
  select (ridreth3);
    when (1, 2) race = 'Hispanic';
    when (3)    race = 'White';
    when (4)    race = 'Black';
    when (6)    race = 'Asian';
    otherwise   race = 'Other';
  end;
  /* NHANES 2021-2023 has 15 strata and 30 PSUs: 15 design degrees of freedom, so the models keep categories few
     (12 parameters) and leave unknown values out rather than modelling them */
  if indfmpir > . then income = ifc(indfmpir < 2, '<2x', ifc(indfmpir < 4, '2-4x', '4x+'));
  if hiq011 in (1, 2) then insured = ifc(hiq011 = 1, 'yes', 'no');
  if bmxbmi > . then bmi = ifc(bmxbmi < 25, '<25', ifc(bmxbmi < 30, '25-30', '30+'));
run;

/* 3. Survey-weighted proportions, one domain request per call so every output table has the same shape */
%macro est(var, dom, by, out);
  ods output domain=&out;
  proc surveymeans data=nh nobs mean stderr clm nomcar;     /* NOBS: N= would be the population-size option */
    strata sdmvstra; cluster sdmvpsu; weight wtmec2yr;
    var &var;
    domain &dom%if %length(&by) %then *&by;;
  run;
  data &out; set &out; where &dom = 1; length level $11 measure $10; measure = "&var";
    %if %length(&by) = 0 %then level = 'All';
    %else level = catx(' ', %sysfunc(tranwrd(&by, *, %str(,))));;      /* sex*agegrp -> catx(' ', sex, agegrp) */
    keep measure level n mean stderr lowerclmean upperclmean;
  run;
%mend;
%est(htn, adult, , p0)        %est(htn, adult, sex, p1)        %est(htn, adult, agegrp, p2)    %est(htn, adult, sex*agegrp, p3)
%est(aware, hyp, , a0)        %est(aware, hyp, sex, a1)        %est(aware, hyp, agegrp, a2)
%est(meds, hyp, , t0)         %est(meds, hyp, sex, t1)         %est(meds, hyp, agegrp, t2)
%est(controlled, hyp, , c0)   %est(controlled, hyp, sex, c1)   %est(controlled, hyp, agegrp, c2)
data estimates; length measure $10 level $11; set p0 p1 p2 a0 a1 a2 t0 t1 t2 c0 c1 c2; run;

/* 4. Age-adjusted prevalence: direct standardisation to the 2000 US population (NCHS weights, ages 18+) */
data ageadj;
  set p2(in=all) p3;
  length group $5;
  group = ifc(all, 'All', scan(level, 1, ' '));
  age = ifc(all, level, scan(level, 2, ' '));
  w = ifn(age = '18-39', 0.420263, ifn(age = '40-59', 0.357202, 0.222535));
run;
proc sql;
  create table ageadj as select group, sum(w * mean) as adj_mean, sqrt(sum(w * w * stderr * stderr)) as adj_se
    from ageadj group by group;
quit;

/* 5. Who is being missed: unaware of their hypertension (all with it), and uncontrolled (those on treatment) */
%macro model(y, dom, out);
  ods output oddsratios=&out;
  proc surveylogistic data=nh nomcar;
    strata sdmvstra; cluster sdmvpsu; weight wtmec2yr;
    domain &dom;
    class agegrp(ref='18-39') sex(ref='Women') race(ref='White') income(ref='4x+') insured(ref='yes') bmi(ref='<25') / param=ref;
    model &y(event='1') = agegrp sex race income insured bmi;
  run;
  data &out; set &out; where &dom = 1; length model $12; model = "&y"; keep model effect oddsratioest lowercl uppercl; run;
%mend;
%model(unaware, hyp, or1)  %model(uncontrolled, treated, or2)
data odds; length model $12 effect $40; set or1 or2; run;

/* 6. Results to the log as CSV lines (copied into results/ in the repository) */
data _null_;
  length line $300;
  set estimates(in=e) ageadj(in=a) odds(in=o);
  if e then line = catx('|', 'CSV', 'estimate', measure, level, n, put(mean, 10.8), put(stderr, 10.8), put(lowerclmean, 10.8), put(upperclmean, 10.8));
  if a then line = catx('|', 'CSV', 'ageadj', group, put(adj_mean, 10.8), put(adj_se, 10.8));
  if o then line = catx('|', 'CSV', 'odds', model, effect, put(oddsratioest, 10.6), put(lowercl, 10.6), put(uppercl, 10.6));
  putlog line;
run;

/* 7. Charts */
data cascade;
  set estimates(where=(level in ('18-39', '40-59', '60+')));
  length stage $26;
  select (measure);
    when ('htn')        stage = '1 Has it (of all adults)';
    when ('aware')      stage = '2 Aware (of those with it)';
    when ('meds')       stage = '3 Treated';
    when ('controlled') stage = '4 Controlled';
  end;
  pct = 100 * mean; lo = 100 * lowerclmean; hi = 100 * upperclmean;
  format pct 5.1;
run;
proc sgplot data=cascade;
  title 'Hypertension care cascade by age, US adults, NHANES 2021-2023';
  vbarparm category=level response=pct / group=stage groupdisplay=cluster limitlower=lo limitupper=hi datalabel;
  xaxis label='Age'; yaxis label='Percent (95% CI)' max=100;
run;
proc sgpanel data=odds;
  title 'Who is missed: adjusted odds ratios (95% CI)';
  panelby model / columns=2 novarname;
  scatter y=effect x=oddsratioest / xerrorlower=lowercl xerrorupper=uppercl;
  refline 1 / axis=x;
  colaxis type=log label='Odds ratio (log scale)'; rowaxis display=(nolabel) discreteorder=data;
run;
title;
