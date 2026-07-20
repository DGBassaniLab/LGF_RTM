
// PROJECT: 		Stunting Reversal vs RTM
// INVESTIGATORS: 	Diego Bassani, Huma Qamar, Daniel Roth, Kelly Watson
// INSTITUTION: 	The Hospital for Sick Children, Toronto, Canada
// TASK: 			Simulate population datasets and generate supporting Figures
// VERSION DATE: 	Jan 28, 2026

/* Table of contents: 
PART 1: Simulate dataset and generate Figure 1A and 1B
PART 2: Simulate dataset and generate Figure 2
PART 3: Simulate dataset and generate Figure 3
PART 4:	Simulate dataset and generate Figure 4A and 4B
*/

**************************************************
// PART 1: Simulate dataset and generate Figure 1
**************************************************

*****************************************
* Part 1A - without faltering on LAZ scale
***************************************** 
clear all
set obs 100000
set seed 998876

* Set Mean HAZ =-1 (constant across all ages)
matrix MD = (-1, -1, -1, -1, -1)
* Set SD =1  (constant across all ages)
matrix SDD = (1, 1, 1, 1, 1)

* Set correlation structure between means at each timepoint: 
** Using empirical correlations from a Bangladeshi trial cohort 
** Matrix for: 3,6,12,18,24-months
matrix C =      (1.00, 0.82, 0.80, 0.77, 0.75)
matrix C =   C\ (0.82, 1.00, 0.85, 0.82, 0.80)
matrix C =   C\ (0.80, 0.85, 1.00, 0.92, 0.89)
matrix C =   C\ (0.77, 0.82, 0.92, 1.00, 0.99)
matrix C =   C\ (0.75, 0.80, 0.89, 0.99, 1.00)

*Use the mean, SD, and correlation matrix to simulate a population cohort from 3 to 24 months
drawnorm haz3 haz6 haz12 haz18 haz24, means(MD) sds(SDD) corr(C)

*Verify that the resulting dataset has the expected mean, sd and correlations
summarize
pwcorr haz3 haz6 haz12 haz18 haz24

*Generate an indicator variable ('tag') for stunted at 3 months
gen partID=_n
gen tag = (haz3 <-2) // tag obs with '1' if HAZ <-2 at 3-months

** Tag obs. that 'reverse' (HAZ>=-2) at 6, 12, 18, or 24 mos
* 'Reversal' always calculated from 3-month baseline
gen rev6 = 1 if tag ==1 & haz6 >=-2
gen rev12 = 1 if tag ==1 & haz12 >=-2
gen rev18 = 1 if tag ==1 & haz18 >=-2
gen rev24 = 1 if tag ==1 & haz24 >=-2

** Calculate prop. of obs. that 'reverse' at all timepoints after 3-months
count if tag==1 // denominator used for stunting reversal proportion 
global obs=r(N)
count if rev6==1
local reversal6 = round((r(N)/$obs )*100, 1)
count if rev12==1
local reversal12 = round((r(N)/$obs )*100, 1)
count if rev18==1
local reversal18 = round((r(N)/$obs )*100, 1)
count if rev24==1
local reversal24 = round((r(N)/$obs )*100, 1)

/*
*Consider methods that may be used to estimate catch-up 
*	(using 3 to 6 month interval as example)
*--------------------------------------------------------------
*Mean HAZ at 3 and 6 months among children stunted at 3 months
sum haz3 haz6 if tag==1
quietly sum haz3 if tag==1
scalar haz3_mean_stunted=r(mean)
quietly sum haz6 if tag==1
scalar haz6_mean_stunted=r(mean)
quietly corr haz3 haz6 
scalar rho_haz3_haz6=r(rho)

*Estimate "RTM effect" per Linden 2013 (https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-13-119)
rtmci haz3 haz6, cutoff(-2)

*Conditional model: Regression of HAZ at 6 months on HAZ at 3 months
reg haz6 haz3

*Using coefficients from cond model above, predict mean haz6 for children-stunted-at-3 months (based on their mean HAZ at 3 months)
display _b[_cons]+(_b[haz3]*haz3_mean_stunted)

*Reduced form of Tim Cole's formula for catch-up: Z2 - rZ1 (as per Cameron et al. 2005)
di haz3_mean_stunted*rho_haz3_haz6
di haz6_mean_stunted-(haz3_mean_stunted*rho_haz3_haz6)

*Reduced form of the Cole formula by using correlation only in children-stunted-at-3 months group
corr haz3 haz6 if tag==1
scalar rho_haz3_haz6_stunted=r(rho)
di haz3_mean_stunted*rho_haz3_haz6_stunted
di haz6_mean_stunted-(haz3_mean_stunted*rho_haz3_haz6_stunted)

*--------------------------------------------------------------
*/

*DRAW FIGURE 1a
*Select 0.5% of sample 1% of  observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long haz, i(partID) j(age_cat)

* Plot the data
gen jitter = runiform(-0.5, 0.5)
gen agejitter = age_cat+jitter
gen mean_haz = -1

twoway (scatter haz agejitter if tag == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag == 1, mcolor(purple) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 6 & rev6 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 6 & rev6 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 12 & rev12 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 12 & rev12 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 18 & rev18 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 18 & rev18 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 24 & rev24 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 24 & rev24 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (line mean_haz age_cat,  sort lcolor(gs3%30) lpattern(-) lwidth(thin)) ///
	   (scatteri 4.75 3 4.75 6, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 4.75 3 4.75 6, ///
	   recast(dropline) base(4.6) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 4.25 3 4.25 12, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 4.25 3 4.25 12, ///
	   recast(dropline) base(4.1) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 3.75 3 3.75 18, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 3.75 3 3.75 18, ///
	   recast(dropline) base(3.6) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 3.25 3 3.25 24, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 3.25 3 3.25 24, ///
	   recast(dropline) base(3.1) lw(thin) mc(none) lc(black) lp(solid)), ///
	   legend(order(2 3 4) label(2 "HAZ<-2 at 3 months") label(3 "HAZ<-2 at 3 months but HAZ>=-2 at the age shown") ///
	   label(4 "HAZ<-2 at 3 months and at the age shown") pos(6) col(1) size(2) region(lstyle(none))) ///
       xlabel(3 6 12 18 24, nogrid) ylabel(-4(2)5, nogrid) ysc(range(-4 4.5)) ///
	   title("") ///
       xtitle("Age (months)") ytitle("HAZ") /// 
	   graphregion(color(white)) ///
	   text(-1.7 7.75 "`reversal6'%", color(red) size(3)) /// 
	   text(-1.7 13.75 "`reversal12'%", color(red) size(3)) /// 
	   text(-1.7 19.75 "`reversal18'%", color(red) size(3)) ///
	   text(-1.7 25.75 "`reversal24'%", color(red) size(3)) /// 
	   text(-0.84 0.75 "Mean", color(gs3%60) size(2)) ///
	   yline(-2, lwidth(thin) lpattern(solid) lcolor(purple%30)) ///
	   yscale(range(-5 5)) ///
	   xscale(range(0 26)) ///
	   text(5 4.5 "r = 0.82", color(black) size(2.5)) /// 
	   text(4.5 9 "r = 0.80", color(black) size(2.5)) ///
	   text(4 15 "r = 0.77", color(black) size(2.5)) ///
	   text(3.5 21 "r = 0.75", color(black) size(2.5)) ///
	   note("A", pos(10) size(medlarge)) ///
	   name(figure1a, replace)
	
graph save "StuntR_Fig1a.gph", replace
graph export "StuntR_Fig1a.png", as(png) replace width(10000)

*************************************** 
* Part 1B - with faltering on LAZ scale
***************************************  

clear
set obs 100000
set seed 998877

* Set Mean HAZ =-1, reducing -0.25 at each timepoint (loosely based on https://www.nature.com/articles/s41586-023-06418-5/figures/13):
matrix MD = (-1, -1.25, -1.5, -1.75, -2)
* Set SD =1 (same as non-faltering scenario)
matrix SDD = (1, 1, 1, 1, 1)

* Set correlation structure between means at each timepoint: 
** Using empirical correlations from a Bangladeshi trial cohort 
** Matrix for: 3,6,12,18,24-months
matrix C =      (1.00, 0.82, 0.80, 0.77, 0.75)
matrix C =   C\ (0.82, 1.00, 0.85, 0.82, 0.80)
matrix C =   C\ (0.80, 0.85, 1.00, 0.92, 0.89)
matrix C =   C\ (0.77, 0.82, 0.92, 1.00, 0.99)
matrix C =   C\ (0.75, 0.80, 0.89, 0.99, 1.00)

*** Simulate population dataset using empirical correlation structure
drawnorm haz3 haz6 haz12 haz18 haz24, means(MD) sds(SDD) corr(C)

*Verify that the resulting dataset has the expected mean, sd and correlations
summarize
pwcorr haz3 haz6 haz12 haz18 haz24

*Generate an indicator variable ('tag') for stunted at 3 months
gen partID=_n
gen tag = (haz3 <-2) // tag obs with '1' if HAZ <-2 at 3-months

** Tag obs. that 'reverse' (HAZ>=-2) at 6, 12, 18, or 24 mos
* 'Reversal' always calculated from 3-month baseline
gen rev6 = 1 if tag ==1 & haz6 >=-2
gen rev12 = 1 if tag ==1 & haz12 >=-2
gen rev18 = 1 if tag ==1 & haz18 >=-2
gen rev24 = 1 if tag ==1 & haz24 >=-2

** Calculate prop. of obs. that 'reverse' at all timepoints after 3-months
count if tag==1 // denominator used for stunting reversal proportion 
global obs=r(N)
count if rev6==1
local reversal6 = round((r(N)/$obs )*100, 1)
count if rev12==1
local reversal12 = round((r(N)/$obs )*100, 1)
count if rev18==1
local reversal18 = round((r(N)/$obs )*100, 1)
count if rev24==1
local reversal24 = round((r(N)/$obs )*100, 1)

// save mean haz for figure 1b
foreach age of numlist 3 6 12 18 24 {
	
su haz`age'
gen mean_haz`age' = r(mean) 
	
}
/*
*Methods to estimate catch-up (same as for figure 1a)
*	(using 3 to 6 month interval as example)
*--------------------------------------------------------------
*Mean HAZ at 3 and 6 months among children stunted at 3 months
sum haz3 haz6 if tag==1
quietly sum haz3 if tag==1
scalar haz3_mean_stunted=r(mean)
quietly sum haz6 if tag==1
scalar haz6_mean_stunted=r(mean)
quietly corr haz3 haz6 
scalar rho_haz3_haz6=r(rho)

*Compare to estimates from RTM effect macro
rtmci haz3 haz6, cutoff(-2)

*Conditional model: Regression of HAZ at 6 months on HAZ at 3 months
reg haz6 haz3

*Using coefficients from cond model above, predict mean haz6 for children-stunted-at-3 months (based on their mean HAZ at 3 months)
display _b[_cons]+(_b[haz3]*haz3_mean_stunted)

*Reduced form of Tim Cole's formula for catch-up: Z2 - rZ1 (as per Cameron et al. 2005)
di haz3_mean_stunted*rho_haz3_haz6
di haz6_mean_stunted-(haz3_mean_stunted*rho_haz3_haz6)

*Reduced form of the Cole formula by using correlation only in children-stunted-at-3 months group
corr haz3 haz6 if tag==1
scalar rho_haz3_haz6_stunted=r(rho)
di haz3_mean_stunted*rho_haz3_haz6_stunted
di haz6_mean_stunted-(haz3_mean_stunted*rho_haz3_haz6_stunted)

*--------------------------------------------------------------
*/

*DRAW FIGURE 1b
*sample 0.5% of the observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations.
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long haz mean_haz, i(partID) j(age_cat)

// Drop 2 observations with very low haz for the purpose of generating the figure (all calculations retained these values) 
drop if haz < -5.75

* Plot the data
gen jitter = runiform(-0.5, 0.5)
gen agejitter = age_cat+jitter

twoway (scatter haz agejitter if tag == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag == 1, mcolor(purple) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 6 & rev6 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 6 & rev6 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 12 & rev12 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 12 & rev12 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 18 & rev18 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 18 & rev18 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 24 & rev24 == 1, mcolor(red) msize(tiny) msymbol(circle)) ///
	   (scatter haz agejitter if age_cat == 24 & rev24 == . & tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
	   (line mean_haz age_cat,  sort lcolor(gs3%30) lpattern(-) lwidth(thin)) ///
	   (scatteri 4.75 3 4.75 6, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 4.75 3 4.75 6, ///
	   recast(dropline) base(4.6) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 4.25 3 4.25 12, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 4.25 3 4.25 12, ///
	   recast(dropline) base(4.1) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 3.75 3 3.75 18, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 3.75 3 3.75 18, ///
	   recast(dropline) base(3.6) lw(thin) mc(none) lc(black) lp(solid)) ///
	   (scatteri 3.25 3 3.25 24, recast(line) ///
	   lw(thin) mc(none) lc(black) lp(solid)) (scatteri 3.25 3 3.25 24, ///
	   recast(dropline) base(3.1) lw(thin) mc(none) lc(black) lp(solid)), ///
	   legend(order(2 3 4) label(2 "HAZ<-2 at 3 months") label(3 "HAZ<-2 at 3 months but HAZ>=-2 at the age shown") ///
	   label(4 "HAZ<-2 at 3 months and at the age shown") pos(6) col(1) size(2) region(lstyle(none))) ///
       xlabel(3 6 12 18 24, nogrid) ylabel(-4(2)4, nogrid) ysc(range(-4 4.5)) ///
	   title("") ///
       xtitle("Age (months)") ytitle("HAZ") /// 
	   graphregion(color(white)) ///
	   text(-1 7.75 "`reversal6'%", color(red) size(3)) /// 
	   text(-1 13.75 "`reversal12'%", color(red) size(3)) /// 
	   text(-1 19.75 "`reversal18'%", color(red) size(3)) ///
	   text(-1 25.75 "`reversal24'%", color(red) size(3)) ///
	   text(-0.84 0.75 "Mean", color(gs3%60) size(2)) /// 
	   yline(-2, lwidth(thin) lpattern(solid) lcolor(purple%30)) ///
	   yscale(range(-5 5)) ///
	   xscale(range(0 26)) ///
	   text(5 4.5 "r = 0.82", color(black) size(2.5)) /// 
	   text(4.5 9 "r = 0.80", color(black) size(2.5)) ///
	   text(4 15 "r = 0.77", color(black) size(2.5)) ///
	   text(3.5 21 "r = 0.75", color(black) size(2.5)) ///
	   note("B", pos(10) size(medlarge)) ///
	   name(figure1b, replace)
	   
	   
graph save "StuntR_Fig1b.gph", replace
graph export "StuntR_Fig1b.png", as(png) replace width(5000) height(3000)

/***** Combine panels A and B *******/
grc1leg "StuntR_Fig1a.gph" figure1b, col(1) graphregion(color(white))
gr display, xsize(8) ysize(12)
gr save "StuntR_Fig1_panel.gph",  replace
gr export "StuntR_Fig1_panel.png", as(png) height(10000) replace

**************************************************
// PART 2: Simulate dataset and generate Figure 2
**************************************************

clear all
set obs 100000
set seed 113432

* Set Mean HAZ =0 at each timepoint:
matrix M2 = (0,0,0,0,0)

* Set SD =1 at each timepoint:
matrix SD2 = (1,1,1,1,1)

* Set correlation structure between means at each timepoint: 
** Using empirical correlations from a Bangladeshi trial cohort 
** Matrix for: 3,6,12,18,24-months
matrix C =  (1.00, 0.82, 0.80, 0.77, 0.75)
matrix C =   C\ (0.82, 1.00, 0.85, 0.82, 0.80)
matrix C =   C\ (0.80, 0.85, 1.00, 0.92, 0.89)
matrix C =   C\ (0.77, 0.82, 0.92, 1.00, 0.99)
matrix C =   C\ (0.75, 0.80, 0.89, 0.99, 1.00)

*** Simulate population dataset using empirical correlation structure
drawnorm haz3 haz6 haz12 haz18 haz24, means(M2) sds(SD2) corr(C)

*** Verify that the resulting dataset has the expected mean, sd and correlations
summarize
pwcorr haz3 haz6 haz12 haz18 haz24

gen partID=_n

** Tag obs. that have a HAZ less than -1 ('tag') or greater than +1 ('tag2') at 3-months
gen tag = (haz3 < -1)
gen tag2 = (haz3 > 1)

** Tag obs. that have a HAZ greater than or equal to -1 or less than or equal to +1 at 3-months but >+1 or <-1 at 6 months ('tag3')
gen tag3 = ((haz3 >= -1 & haz3 <= 1) & (haz6 > 1 | haz6 < -1))

** Tag obs. that cross above threshold of -1 (up toward mean=0) at 6,12,18,or24-mos
gen up6 = 1 if tag ==1 & haz6 >=-1
gen up12 = 1 if tag ==1 & haz12 >=-1
gen up18 = 1 if tag ==1 & haz18 >=-1
gen up24 = 1 if tag ==1 & haz24 >=-1


** Calculate prop. of obs. that 'reverse' (move 'up') at all timepoints after 3-months
count if tag==1 // denominator used when calculating props.
global obs=r(N)
count if up6==1
local poschange6 = round((r(N)/$obs )*100,1)
count if up12==1
local poschange12 = round((r(N)/$obs )*100,1)
count if up18==1
local poschange18 = round((r(N)/$obs )*100,1)
count if up24==1
local poschange24 = round((r(N)/$obs )*100,1)

** Tag obs. that cross below threshold of +1 (down toward mean=0) at 6,12,18,or24-mos
gen down6 = 1 if tag2 ==1 & haz6 <=1
gen down12 = 1 if tag2 ==1 & haz12 <=1
gen down18 = 1 if tag2 ==1 & haz18 <=1
gen down24 = 1 if tag2 ==1 & haz24 <=1


** Calculate prop. of obs. that move 'down' at all timepoints after 3-months
count if tag2==1 // denominator used when calculating props.
global obs=r(N)
count if down6==1
local negchange6 = round((r(N)/$obs )*100,1)
count if down12==1
local negchange12 = round((r(N)/$obs )*100,1)
count if down18==1
local negchange18 = round((r(N)/$obs )*100,1)
count if down24==1
local negchange24 = round((r(N)/$obs )*100,1)

** Calculate prop. of obs. that move from >=-1 and =<+1 at 3 months to tails at 6-months
count if tag==0 & tag2==0 // denominator
global obs=r(N)
count if tag3==1 & haz6 > 1 // subset that moved from >=-1 and =<+1 at 3 months to >+1 at 6-months
local uptotail6 = round((r(N)/$obs)*100,1)
count if tag==0 & tag2==0 // denominator
global obs=r(N)
count if tag3==1 & haz6 < -1 // subset that moved from >=-1 and =<+1 at 3 months to <-1 at 6-months
local downtotail6 = round((r(N)/$obs)*100,1)

*sample 0.5% of the observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations.
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long haz, i(partID) j(age_cat)

* Plot the data
gen jitter = runiform(-0.5, 0.5)
gen agejitter = age_cat+jitter


twoway (scatter haz agejitter if tag == 0 & tag2 == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag == 1, mcolor(blue*0.80) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag2 == 1, mcolor(dkgreen*0.7) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag3 == 1 & (age_cat==3 | age_cat==6), mcolor(pink%40) msize(tiny) msymbol(Oh)), ///
	   legend(order(2 3 4) row(1) label(2 "HAZ < -1 at 3 months") label(3 "HAZ > +1 at 3 months") label(4 "-1 ≥ HAZ ≤ +1 at 3 months, but <-1 or >+1 at 6 months") pos(6) /// 
	   size(2.5) region(lcolor(black) lwidth(0.08) )) ///
	   xlabel(3 6 12 18 24, nogrid) ylabel(, nogrid) ///
       title("") ///
       xtitle("Age (months)") ytitle("HAZ") ///
	   text(0.16 0.12 "Mean", color(gs3%60) size(2)) /// 
	   yline(0, lwidth(thin) lcolor(gs3%30)) ///
	   graphregion(color(white)) ///
	   text(-0.7 7.3 "`poschange6'%", color(blue) size(2)) /// 
	   text(0.7 7.3 "`negchange6'%", color(dkgreen) size(2)) /// 
	   text(-0.7 13.3 "`poschange12'%", color(blue) size(2)) /// 
	   text(0.7 13.3 "`negchange12'%", color(dkgreen) size(2)) ///
	   text(-0.7 19.3 "`poschange18'%", color(blue) size(2)) /// 
	   text(0.7 19.3 "`negchange18'%", color(dkgreen) size(2)) ///
	   text(-0.7 25.3 "`poschange24'%", color(blue) size(2)) /// 
	   text(0.7 25.3 "`negchange24'%", color(dkgreen) size(2)) ///
	   text(1.4 7.3 "`uptotail6'%", color(pink) size(1.75)) ///
	   text(-1.4 7.3 "`downtotail6'%", color(pink) size(1.75)) ///
	   text(-1 -0.8 "-1", color(blue%60) size(2)) /// 
	   yline(-1, lwidth(thin) lcolor(blue%30)) ///
	   text(1 -0.8 "+1", color(dkgreen%60) size(2)) /// 
	   yline(1, lwidth(thin) lcolor(dkgreen%30)) ///
	   xscale(range(0 26))
	   
graph save "StuntR_Fig2.gph", replace
graph export "StuntR_Fig2.png", as(png) replace width(10000)


**************************************************
// PART 3: Simulate dataset and generate Figure 3
**************************************************

clear all
set obs 100000
set seed 08252107

* Set Mean HAZ at each timepoint:
matrix M3 = (-1,-1)

* Set SD of Mean HAZ at each timepoint:
matrix SD3 = (1,1)

* Set correlation structure between means at each timepoint: 
** Using empirical correlations from a Bangladeshi trial cohort 
** Matrix for: 3,24-months
matrix C =  (1.00, 0.75)
matrix C =   C\ (0.75, 1.00)

*** Simulate population dataset using empirical correlation structure
drawnorm haz3 haz24, means(M3) sds(SD3) corr(C)

*** Verify that the resulting dataset has the expected mean, sd and correlations
summarize
pwcorr haz3 haz24

*** Calculate correlation of HAZ3 with change from 3 to 24 months
gen delta3_24 = haz24 - haz3 
pwcorr haz3 haz24 delta3_24
drop delta3_24

gen partID=_n

* Tag obs. if HAZ at 3-mos less than -2 but greater than or equal to -2.5
gen tag = (haz3<-2 & haz3>=-2.5)

* Tag obs. if HAZ at 3-mos less than -2.5
gen tag2 = (haz3<-2.5)

* Calculate avg HAZ increase for obs. w/ HAZ b/w -2 & -2.5 at 3mos
gen diff = haz24 - haz3 if tag==1
summarize diff
local meandiff = round(r(mean), 0.01)
local SD = round(r(sd), 0.01)

* Calculate avg HAZ increase for obs. w/ HAZ less than -2.5 at 3mos
gen diff2 = haz24 - haz3 if tag2==1
summarize diff2
local meandiff2 = round(r(mean), 0.01)
local SD2 = round(r(sd), 0.01)

* Calculate 'reversal' prop. for obs. with HAZ b/w -2 & -2.5 at 3mos
count if tag==1 						// denominator
global obs=r(N)
gen bluerev = 1 if tag ==1 & haz24 >=-2 // blue corresponds to colour on plot
count if bluerev==1
local reversalblue = round((r(N)/$obs )*100, 1)

* Calculate 'reversal' prop. for obs. with HAZ less than -2.5 at 3mos
count if tag2==1 						// denominator
global obs=r(N)
gen redrev = 1 if tag2 ==1 & haz24 >=-2 // red corresponds to colour on plot
count if redrev==1
local reversalred = round((r(N)/$obs )*100, 1)

*sample 0.5% of the observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations.
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long haz, i(partID) j(age_cat)

* Plot the data
gen jitter = runiform(-1.5,1.5)
gen agejitter = age_cat+jitter
twoway (scatter haz agejitter if tag == 0 & tag2 == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
       (scatter haz agejitter if tag2 == 1, mcolor(red) msize(tiny) msymbol(circle)), ///
       legend(order (2 3) label(2 "HAZ<-2 and HAZ>=-2.5 at 3 months") label(3 "HAZ<-2.5 at 3 months") ///
	   position(6) col(2) size(2) region(lcolor(black) lwidth(0.08))) ///
       xlabel(3 24, nogrid) ylabel(, nogrid) ///
	   ylabel(-4(2)2) ///
	   xscale(range(-6 33)) ///
       title("") ///
       xtitle("Age (months)") ytitle("HAZ") ///
	   text(-0.88 -5.6 "Mean", color(gs3%60) size(2)) /// 
	   text(-1 -8 "-1", color(gs3%60) size(2)) /// 
	   yline(-1, lwidth(thin) lcolor(gs3%30)) ///
	   graphregion(color(white)) /// 
	   xsize(4) /// 
	   text(1 13.5 "Mean △HAZ (SD): +0`meandiff' (0`SD')", color(blue) size(2.5)) /// 
	   text(0.5 13 "Stunting 'reversal': `reversalblue'%", color(blue) size(2.5)) ///
	   text(0 13.5 "Mean △HAZ (SD): +0`meandiff2' (0`SD2')", color(red) size(2.5)) ///
	   text(-0.5 13 "Stunting 'reversal': `reversalred'%", color(red) size(2.5)) ///
	   yline(-2, lwidth(thin) lcolor(blue%30)) ///
	   text(-2.5 -8.2 "-2.5", color(red%60) size(2)) /// 
	   yline(-2.5, lwidth(thin) lcolor(red%30))
	   
graph save "StuntR_Fig3.gph", replace
graph export "StuntR_Fig3.png", as(png) replace height(7000)

*****************************************************************************

*****************************************************
// PART 4A: Simulate dataset and generate Figure 4A
*****************************************************

clear all
set obs 100000
set seed 80252687

* Set mean length for girls at 3mo and 24m (corresponding to HAZ = -1). 

matrix M4 = (57.7, 83.2)

* Set SD of Mean length at each timepoint:
matrix SD4 = (2.1051, 3.2267)

* Set correlation structure between means at each timepoint: 
** Using same empirical correlations as for LAZ

matrix C4 =  (1.00, 0.75)
matrix C4 =   C4\ (0.75, 1.00)

*** Simulate population dataset
drawnorm length3 length24, means(M4) sds(SD4) corr(C4)

*** Verify that the resulting dataset has the expected mean, sd and correlations
summarize, det
pwcorr length3 length24

*** Calculate correlation of length3 with change from 3 to 24 months
gen delta3_24 = length24 - length3 
sum delta3_24
pwcorr length3 length24 delta3_24
drop delta3_24

gen partID = _n

*Generate variables for 5th and 10th percentiles at 3 months
egen length3_p5 = pctile(length3), p(5)
egen length3_p10 = pctile(length3), p(10)

*Generate variables for 10th percentile at 24 months
egen length24_p10 = pctile(length24), p(10)

* Tag obs. if length at 3-mos less than 10th but greater than (or equal to) 5th
gen tag = (length3 < length3_p10 & length3 >= length3_p5)

* Tag obs. if length at 3-mos less than 5th percentile
gen tag2 = (length3 < length3_p5)

* Calculate avg length increase for obs. w/ length between 5th and 10th percentile at 3mos
gen diff = length24 - length3 if tag == 1
summarize diff
local meandiff = round(r(mean), 0.1)
local SD = round(r(sd), 0.1)

* Calculate avg length increase for obs. w/ length less than 5th percentile at 3mos
gen diff2 = length24 - length3 if tag2 == 1
summarize diff2
local meandiff2 = round(r(mean), 0.1)
local SD2 = round(r(sd), 0.1)

* Calculate 'catch-up' prop. for obs. with length b/w 5th and 10th percentiles at 3mos
*	In this context, catch-up means that at 24-mo, it is above the 10th percentile
count if tag == 1 										// denominator
global obs = r(N)
gen bluerev = 1 if tag == 1 & length24 >= length24_p10 	// blue corresponds to colour on plot
count if bluerev == 1
local reversalblue = round((r(N)/$obs )*100, 1)

* Calculate 'catchp-up' prop. for obs. with lenth less than 5th percentile at 3mos
count if tag2 == 1 				// denominator
global obs = r(N)
gen redrev = 1 if tag2 == 1 & length24 >= length24_p10 // red corresponds to colour on plot
count if redrev == 1
local reversalred = round((r(N)/$obs )*100, 1)

*sample 0.5% of the observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations.
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long length, i(partID) j(age_cat)

* Plot the data - Figure 4A

gen jitter = runiform(-3,3)
gen agejitter = age_cat+jitter
twoway (scatter length agejitter if tag == 0 & tag2 == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter length agejitter if tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
       (scatter length agejitter if tag2 == 1, mcolor(red) msize(tiny) msymbol(circle)), ///
       legend(order (2 3) label(2 "Length <10th & >=5th percentile at 3 months") /// 
	   label(3 "Length <5th percentile at 3 months") ///
	   position(6) col(2) size(2) region(lcolor(black) lwidth(0.08))) ///
       xlabel(3 24, nogrid) ylabel(, nogrid) ///
	   ylabel(50(10)90) ///
	   xscale(range(-6 33)) ///
       title("") ///
       xtitle("Age (months)") ///
	   ytitle("Length (cm)" "  ") ///
	   text(59 -3.5 "Mean (3 mo)", color(gs3%60) size(2)) ///
	   yline(57.7, lwidth(thin) lcolor(gs3%30)) ///
	   text(84 -3.5 "Mean (24 mo)", color(gs3%60) size(2)) /// 
	   text(80 -1.5 "10th percentile (24 mo)", color(purple%30) size(2)) ///
	   yline(79.1, lwidth(thin) lcolor(purple%30)) ///
	   yline(83.2, lwidth(thin) lcolor(gs3%30)) ///
	   graphregion(color(white)) /// 
	   xsize(4) /// 
	   text(73 13.5 "Mean △Length (SD): +`meandiff'.0 (`SD')", color(blue) size(2.5)) /// 
	   text(71 13 "Catch-up: `reversalblue'%", color(blue) size(2.5)) ///
	   text(67 13.5 "Mean △Length (SD): +`meandiff2' (`SD2')", color(red) size(2.5)) ///
	   text(65 13 "Catch-up: `reversalred'%", color(red) size(2.5))
	   
graph save "StuntR_Fig4A.gph", replace
graph export "StuntR_Fig4A.png", as(png) replace height(9000)

****************************************************
// PART 4B: Simulate dataset and generate Figure 4B
****************************************************

clear all
set obs 100000
set seed 80252687

* Set mean length for girls at 3mo and 24m (corresponding to HAZ = -1). 

matrix M4 = (57.7, 83.2)

* Set SD of Mean length at each timepoint:
matrix SD4 = (2.1051, 3.2267)

* Set correlation structure between means at each timepoint: 
** Using same empirical correlations as for LAZ

matrix C4 =  (1.00, 0.75)
matrix C4 =   C4\ (0.75, 1.00)

*** Simulate population dataset
drawnorm length3 length24, means(M4) sds(SD4) corr(C4)

*** Verify that the resulting dataset has the expected mean, sd and correlations
summarize, det
pwcorr length3 length24

gen partID = _n

*Generate height-for-age difference (HAD) values at 3 and 24 months using median lengths (HAZ=0)
gen HAD3 = length3 - 59.8
gen HAD24 = length24 - 86.4

**Summarize HAD3 and HAD24, delta-HAD and verify correlations
sum HAD3 HAD24, det
gen deltaHAD = HAD24 - HAD3
pwcorr HAD3 HAD24 deltaHAD
summ deltaHAD, det
drop deltaHAD

*For comparison: 
display 57.7 - 59.8		// mean HAD at 3 months
display 83.2 - 86.4		// mean HAD at 24 months
display -3.2 - -2.1		// mean delta HAD

*Generate variables for 5th and 10th percentiles at 3 months
egen HAD3_p5 = pctile(HAD3), p(5)
egen HAD3_p10 = pctile(HAD3), p(10)

* Tag obs. if HAD at 3-mos less than 10th but greater than (or equal to) 5th
gen tag = (HAD3 < HAD3_p10 & HAD3 >= HAD3_p5)

* Tag obs. if HAD at 3-mos less than 5th percentile
gen tag2 = (HAD3 < HAD3_p5)

* Calculate avg HAD increase for obs. w/ length between 5th and 10th percentile at 3mos
gen diff = HAD24 - HAD3 if tag == 1
summarize diff
local meandiff = round(r(mean), 0.1)
local SD = round(r(sd), 0.1)

* Calculate avg length increase for obs. w/ length less than 5th percentile at 3mos
gen diff2 = HAD24 - HAD3 if tag2 == 1
summarize diff2
local meandiff2 = round(r(mean), 0.1)
local SD2 = round(r(sd), 0.1)

* Calculate 'catch-up' prop. for obs. with length b/w 5th and 10th percentiles at 3mos
*	In this context, catch-up means that HAD is lower magnitude (closer to 0) at 24-mo than at 3-mo
count if tag == 1 										// denominator
global obs = r(N)
gen bluerev = 1 if tag == 1 & HAD24 > HAD3		// blue corresponds to colour on plot
count if bluerev == 1
local reversalblue = round((r(N)/$obs )*100, 1)

* Calculate 'catch-up' prop. for obs. with HAD less than 5th percentile at 3mos
count if tag2 == 1 									// denominator
global obs = r(N)
gen redrev = 1 if tag2 == 1 & HAD24 > HAD3		// red corresponds to colour on plot
count if redrev == 1
local reversalred = round((r(N)/$obs )*100, 1)

*sample 0.5% of the observations to produce cleaner plot, but estimates are still based on simulated dataset with 100,000 observations.
local sample = _N * 0.005
sample `sample', count

* Reshape the data to long format
reshape long HAD, i(partID) j(age_cat)

* Plot the data - Figure 4B

gen jitter = runiform(-3,3)
gen agejitter = age_cat+jitter
twoway (scatter HAD agejitter if tag == 0 & tag2 == 0, mcolor(gs13) msize(tiny) msymbol(circle)) ///
       (scatter HAD agejitter if tag == 1, mcolor(blue) msize(tiny) msymbol(circle)) ///
       (scatter HAD agejitter if tag2 == 1, mcolor(red) msize(tiny) msymbol(circle)), ///
       legend(order (2 3) label(2 "HAD <10th & >=5th percentile at 3 months") /// 
	   label(3 "HAD <5th percentile at 3 months") ///
	   position(6) col(2) size(2) region(lcolor(black) lwidth(0.08))) ///
       xlabel(3 24, nogrid) ylabel(, nogrid) ///
	   ylabel(-12(2)6) ///
	   xscale(range(-6 33)) ///
       title("") ///
       xtitle("Age (months)") ///
	   ytitle("Height-for-age difference (cm)" "  ") ///
	   text(-1.7 -3.5 "Mean (3 mo)", color(gs3%60) size(2)) ///
	   yline(-2.09, lwidth(thin) lcolor(gs3%30)) ///
	   text(-2.8 -3.5 "Mean (24 mo)", color(gs3%60) size(2)) /// 
	   yline(-3.19, lwidth(thin) lcolor(gs3%30)) ///
	   graphregion(color(white)) /// 
	   xsize(4) /// 
	   text(6 13.5 "Mean △HAD (SD): `meandiff' (`SD')", color(blue) size(2.5)) /// 
	   text(5 13 "Catch-up: `reversalblue'%", color(blue) size(2.5)) ///
	   text(3 13.5 "Mean △HAD (SD): `meandiff2' (`SD2')", color(red) size(2.5)) ///
	   text(2 13 "Catch-up: `reversalred'%", color(red) size(2.5)) ///
	   name(figure4B, replace)

graph save "StuntR_Fig4B.gph", replace
graph export "StuntR_Fig4B.png", as(png) replace width(3000) height(3000)


/***** Combine panels A and B *******/
grc1leg "StuntR_Fig4A.gph" figure4B, col(1) graphregion(color(white))
gr display, xsize(8) ysize(12)
gr save "StuntR_Fig4_panel.gph",  replace
gr export "StuntR_Fig4.png", as(png) height(9000) replace



*****************************************************************************
*****************************************************************************END
	  
