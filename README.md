hammy
=====

`hammy` is a small program for fitting homoscedastic additive models (HAMs), which are generalizations of linear regression models (but, unsurprisingly, not as general as [generalized additive models](https://en.wikipedia.org/wiki/Generalized_additive_model)). HAMs generalize linear regression models by allowing each predictor to have a nonlinear relationship with the dependent variable to be predicted. Instead of the linear regression model

> *y* = &alpha; + &beta;<sub>1</sub> *x*<sub>1</sub> + &beta;<sub>2</sub> *x*<sub>2</sub> + &hellip; + &beta;<sub>*k*</sub> *x*<sub>k</sub> + &epsilon;

the HAM assumes the relationship

> *y* = *f*<sub>1</sub>(*x*<sub>1</sub>) + *f*<sub>2</sub>(*x*<sub>2</sub>) + &hellip; + *f*<sub>k</sub>(*x*<sub>k</sub>) + &epsilon;

where *f*<sub>i</sub> is a smooth function of the predictor variable *x*<sub>i</sub>. To make the model identifiable, the HAM assumes (without loss of practical generality) the special case

> *y* = &lang;*y*&rang; + *f*<sub>1</sub>(*x*<sub>1</sub>) + *f*<sub>2</sub>(*x*<sub>2</sub>) + &hellip; + *f*<sub>k</sub>(*x*<sub>k</sub>) + &epsilon;

where &lang;*y*&rang; is *y*'s sample mean (so that the total contribution of the smooth functions is zero on average).

There are various ways to estimate the smooth functions.
`hammy` uses Nadaraya-Watson kernel-smoothing to estimate each one-dimensional smooth function *f*<sub>i</sub>.
`nwgks.c` and the standalone `nwgks` program implement that smoothing, and by carrying out that smoothing for multiple predictors `hammy` is able to estimate a complete additive model with homoscedastic errors.
