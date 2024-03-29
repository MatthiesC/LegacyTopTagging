Taken from: https://gitlab.cern.ch/akhukhun/roccor/-/tree/master

Original README follows now:

# README



Muon momentum corrections for 2016, 2017, 2018 end-of-year (EOY) and legacy (UL) samples

Please read the presentation in the muon POG 
https://indico.cern.ch/event/981770/
or the "additional notes" below.

Inputs:
-----------------------------------
UL samples:
2016a:
Run2016{B..F}-21Feb2020(_ver2)_UL2016_HIPM
RunIISummer19UL16MiniAODAPV-106X_mcRun2_asymptotic_preVFP_v8

2016b:
Run2016{F..H}-21Feb2020_UL2016
RunIISummer19UL16MiniAOD-106X_mcRun2_asymptotic_v13

2017:
Run2017{B..F}-09Aug2019_UL2017
RunIISummer19UL17MiniAOD-106X_mc2017_realistic_v6

2018:
Run2018{A..D}-12Nov2019_UL2018
RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1
-----------------------------------
-----------------------------------

EOY:
RoccoR2016.txt:
Run2016{B..H}-17Jul2018
RunIISummer16MiniAODv3 

RoccoR2017.txt:
Run2017{B..F}-31Mar2018
RunIIFall17MiniAODv2

RoccoR2018.txt:
Run2018{A..C}-17Sep2018 + Run2018DPromptReco
RunIIAutumn18MiniAOD
-----------------------------------




Example: 
-------------------------------------------------------------------------------------
//for root macro
#include "RoccoR.cc"
RoccoR  rc("RoccoR2017UL.txt"); 

//Or in cmssw:
#include "RoccoR.h"
RoccoR rc; 
rc.init(edm::FileInPath("path/to/RoccoR2017UL.txt").fullPath()); 
------
//scale factors for momentum of each muon:
double dtSF = rc.kScaleDT(Q, pt, eta, phi, s=0, m=0); //data
double mcSF = rc.kSpreadMC(Q, pt, eta, phi, genPt, s=0, m=0); //(recommended), MC scale and resolution correction when matched gen muon is available
double mcSF = rc.kSmearMC(Q, pt, eta, phi, nl, u, s=0, m=0); //MC scale and extra smearing when matched gen muon is not available
-------------------------------------------------------------------------------------
Here: 
Q is charge
nl is trackerLayersWithMeasurement \\
u is a random number distributed uniformly between 0 and 1 (gRandom->Rndm());
s is error set (default is 0)
m is error member (default is 0, ranges from 0 to nmembers-1)
For MC, when switching to different error sets/members for a given muon, random number (u) should remain unchanged. 

-------------------------------------------------------------------------------------
Following variations are provided to estimate uncertainties. 
-------------------------------------------------------------------------------------
	 set        nmembers	comment
Default  0	    1		default, reference based on madgraph sample, with adhoc ewk (sw2eff and Z width) and Z pt (to match data) weights. 
Stat     1          100         pre-generated stat. replicas; 
Zpt      2          1           derived without reweighting reference pt to data. 
Ewk      3          1           derived without applying ad-hoc ewk weights 
deltaM   4          1           one representative set for alternative profile deltaM mass window 
Ewk2     5          1           reweight the reference from constant to s-dependent Z width 
-------------------------------------------------------------------------------------
For statistical replicas, std. dev. gives uncertainty.
For the rest, difference wrt the cental is assigned as syst. error.
Total uncertainty is calculated as a quadrature sum of all components. 
-------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------


Additional notes:

1) Normally, uncertainties are negligible compared to other uncertainties in the analysis. 
As a simple check, one can compare results with and without applying these corrections. 
If the effect on the analysis is small compared to other uncertainties, then muon calibration 
uncertainties can probably be neglected. 

As an additional option for a quick (single-variation) check, provided functions below return uncertainties, 
evaluated by propagating the variations described above for each muon.
----
double deltaDtSF = rc.kScaleDTerror(Q, pt, eta, phi);
double deltaMcSF = rc.kSpreadMCerror(Q, pt, eta, phi, genPt);
double deltaMcSF = rc.kSmearMCerror(Q, pt, eta, phi, nl, u);
----
Since there is no information here on correlations between charges or different eta/phi bins, these functions 
are not recommended to be used as uncertainties, but only to be used as an estimate of its upper bound to see if it's negligible 
(to do so, scale-factors should be varied by these delta's up or down for different charges and eta-phi regions in 
a most conservative way for a given analysis). 


2) If an analysis is only sensitive to data/mc difference, which is typically the case, it can be more convenient
to apply data systematic variations to MC. Something like:
   kSpreadMC(..., s, m) * kScaleDT(..., 0, 0) / kScaleDT(..., s, m)
while keeping data fixed with central corrections. 


3) Input signal samples use pdg Z mass value with fixed-width propagator parameterization. 
This introduces ~ 34 MeV shift in Z peak position. 
It affects the overall scale by dk~0.0004.

By default we chose the MC peak position since most analyses are not sensitive to this kind of shift
and they can still check their data/mc agreement with default MC. 

However if your analysis is sensitive to this level of shift, 
you may prefer to use this corresponding set as default (Ewk2 in the above table, which corresponds to s=5, m=0)
and exclude this set from the evaluation of systematic uncertainties. 

In case of questions or problems, please email aleko.khukhunaishvili at gmail.com 
