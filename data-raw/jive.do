use ak80.dta, replace

/* Table 2 in AIK99 */
qui: eststo i1: ivregress 2sls lwage (educ=i.qob#i.yob) i.yob
qui: eststo i2: ivregress liml lwage (educ=i.qob#i.yob) i.yob
qui: eststo r1: ivregress 2sls lwage (educ=i.qob#i.yob) i.yob, robust
qui: eststo r2: ivregress liml lwage (educ=i.qob#i.yob) i.yob, robust
esttab i1 i2 r1 r2, b(%11.10f) se(%11.10f) drop(*yob*)

/* 0.0891154613 (0.0161098202) (0.0162120317)
   0.0928764165 (0.0177441446) (0.0196323640) */

qui: eststo i3: ivregress 2sls lwage (educ=i.qob#i.yob i.qob#i.sob) i.yob i.sob
qui: eststo i4: ivregress liml lwage (educ=i.qob#i.yob i.qob#i.sob) i.yob i.sob
esttab i3 i4, b(%11.10f) se(%11.10f) drop(*yob* *sob*)
/*    0.0928180625 (0.0093013345)
      0.1063979826 (0.0116383736) */


/* TODO Small sample */
qui: eststo i1: ivregress liml lwage (education=i.qob) i.yob in 1/50
