timer on 1
use ak80.dta, replace

/* Table V Columns (1) and (2) */
qui: eststo r1: reg  lwage education i.yob
qui: eststo r2: reg  lwage education i.yob, robust
qui: eststo r3: reg  lwage education i.yob in 1/100
qui: eststo r4: reg  lwage education i.yob in 1/100, robust

esttab r1 r2 r3 r4, b(%11.10f) se(%11.10f) drop(*yob*)
/* 0.0710810458 (0.0003390067) [0.0003814625] */
/* 0.1004438463 (0.0346482285) [0.0233448020] */

qui: eststo i1: ivregress 2sls lwage  (education=i.qob) i.yob in 1/100
qui: eststo i2: ivregress 2sls lwage  (education=i.qob) i.yob in 1/100, robust
qui: eststo i3: ivregress 2sls lwage  (education=i.yob##i.qob) i.yob
qui: eststo i4: ivregress 2sls lwage  (education=i.yob##i.qob) i.yob, robust

esttab i1 i2 i3 i4, b(%11.10f) se(%11.10f) drop(*yob*)
/*-0.1609028242 (0.3277804536) [0.4120184496] */
/* 0.0891154613 (0.0161098202) [0.0162120317] */


/* MBTSLS */
/* first get number of exog and endog variagbles */
qui: ivreg2 lwage (education=i.qob) i.yob in 1/100
local k1 = 1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
qui: eststo i1: ivreg2 lwage (education=i.qob) i.yob in 1/100, kclass(`k1')
qui: eststo i2: ivreg2 lwage (education=i.qob) i.yob in 1/100, kclass(`k1') robust
qui: ivreg2 lwage (education=i.yob##i.qob) i.yob
local k2 = 1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
qui: eststo i3: ivreg2 lwage (education=i.yob##i.qob) i.yob, kclass(`k2')
qui: eststo i4: ivreg2 lwage (education=i.yob##i.qob) i.yob, kclass(`k2') robust

esttab i1 i2 i3 i4, b(%11.10f) se(%11.10f) drop(*yob*)
/* 0.3505927902 (0.0504391651) [0.3894877272] */
/* 0.0937333665 (0.0180984698) [0.0204147326] */



qui: eststo i1: ivregress liml lwage  (education=i.qob) i.yob in 1/100
qui: eststo i2: ivregress liml lwage  (education=i.qob) i.yob in 1/100, robust
qui: eststo i3: ivregress liml lwage  (education=i.yob##i.qob) i.yob
qui: eststo i4: ivregress liml lwage  (education=i.yob##i.qob) i.yob, robust

esttab i1 i2 i3 i4, b(%11.10f) se(%11.10f) drop(*yob*)
/*-0.2054464334 (0.3788818685) [0.4981386013] */
/* 0.0928764165 (0.0177441446) [0.0196323640] */

/* Table VII */

qui: eststo i1: reg  lwage education i.yob i.sob
qui: eststo i2: reg  lwage education i.yob i.sob, robust
qui: eststo i3: ivregress 2sls lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob
qui: eststo i4: ivregress 2sls lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob, robust
qui: eststo i5: ivregress liml lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob
qui: eststo i6: ivregress liml lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob, robust
esttab i1 i2 i3 i4 i5 i6, b(%11.10f) se(%11.10f) drop(*yob* *sob*)
/* 0.0673389705 (0.0003464258) [0.0003883466] */
/* 0.0928180625 (0.0093013345) [0.0096641481] */
/* 0.1063979834 (0.0116383738) [0.0149803714] */

/* MBTSLS */
/* first get number of exog and endog variagbles */
qui: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob
local k3 =  1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
qui: eststo i1: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob, kclass(`k3')
qui: eststo i2: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob, kclass(`k3') robust
esttab i1 i2, b(%11.10f) se(%11.10f) drop(*yob* *sob*)
/* 0.1089429483 (0.0120411994) [0.0159979308] */

timer off 1
timer list
