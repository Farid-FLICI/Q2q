# Q2q: Interpolating mortality rates at all ages

## Description 

Mortality Rates are usually published following an abridged description, i.e., by age groups 0, [1, 5[, [5, 10[, [10, 15[ and so on. For some applications, a detailed (single) ages description is required. Despite the huge number of the proposed methods in the literature, there is a limited number of methods ensuring a high performance at lower and higher ages in the same time. For example, the 6-terms `Lagrange` interpolation function is well adapted to mortality  interpolation at lower ages (with unequal intervals) but is not adapted to higher ages. On the other hand, the `Karup-King` method allows a good performance at higher ages but not adapted to lower ages.  The package `Q2q` allows combining both the two methods to allow interpolating mortality rates at all ages. First, it starts by implementing each method separately, then the resulted curves are joined based on a 5-age averaged error between the two curves.


## Installation

You can install the released version of Q2q from [github](https://github.com/Farid-FLICI/Q2q/blob/master/Q2q_0.1.0.tar.gz) with:

## Functions

The package `Q2q` provides two functions `getqxt()` and `getqx()`.

### getqxt()

The `getqxt()` functions allows interpolating age specific mortality rates $$ASMRs$$ starting from an abridged mortality surface. This later should provides the five ages mortality rates, usually noted as $$_nQ_{x,t}$$ in the literature, with $$x$$ representing the age and $t$ the year and $n$ the length of the age interval which is set to be $5$ except for age $0$ and the age group $1-5$.

The general formulation of the function is `getqxt(Qxt, nag, t)` with `Qxt` is the matrix of the five ages mortality rates with no header and identification column. The number of rows in the matrix should correspond to the number of age groups refereed as `nag` in the function. `t` corresponds to the number of years contained in the mortality matrix, and it should be equal to the number of columns in `Qxt`.

The function results principally in `qxt` which represents the matrix of the age-specific mortality rates for age $$x$$ and year $$t$$. This matrix had resulted from the junction of two initial matrices `qxtl` and `qxtk` obtained with the Lagrange and the Karup-King methods respectively. These two matrices are also provided as function returns. additionally, for each year $$t$$, the junction age is provided in a vector `$jonct_ages`. The functions also returns the survivorship matrix `lxt` and the matrix of theoretical deaths `qxt`.

### getqxt()

The `getqx()` function is a special case of `getqxt()` with `t=1`.

```{r}
getqx <- function(Qx,nag) {
                           getqxt(Qxt=Qx, nag, t=1)
                          }
```
## Examples

lets consider the Algerian mortality surface from 1977 to 2014, for males, by 5 ages groups.

```{r}
LT <- read.table("https://raw.githubusercontent.com/Farid-FLICI/farid-flici.github.io/y/completed_Mortality_Rates_Males_nQx.txt", dec=".", sep="\t", quote="" )

LT <- data.matrix(LT)
```

### Example 1

lets consider a case of an annual life table (Algeria, males, year 1995).

```{r echo=TRUE, include=TRUE}
Ax <- LT[ 2:18  , 18  ]
Ax
```

We plot the mortality curve.

```{r include=TRUE, fig.align="center",fig.width=10,fig.height=6}
plot(x=c(0,1, seq(5, (length(Ax)-2)*5, by=5)), y=log(Ax), type="b", xlab="Age", ylab="log(nQx)", main="Figure 1. Mortality Rates (Algeria, males, 1995)" ) 
```

We interpolate $$q_{x,t}$$ using `getqx()`

```{r}
library(Q2q)
interpolated_curve <-getqx(Qx=Ax, nag=17)
str(interpolated_curve)
```

plot the interpolated mortality curve using `plot_qx()`

```{r include=TRUE , fig.align="center",fig.width=10,fig.height=6}
plot(x=c(0:79), y=log(interpolated_curve$qx), type="b", xlab="Age", ylab="log(qx)", main="Figure 2. Interpolated Mortality Rates (Algeria, males, 1995")
```

```{r include=TRUE , fig.align="center",fig.width=10,fig.height=6}
plot(x=c(5:79), y=log(interpolated_curve$qxk[6:80]), type="p", cex=0.75, pch=1,xlab="Age", ylab="log(qx)", main="Figure 3. Interpolated Mortality Rates (Algeria, males, 1995) - Junction age", xlim=c(0,80))
lines(x=c(0:69), y=log(interpolated_curve$qxl), lwd=2)
lines(x=rep(interpolated_curve$jonct_age, 2), y=c(min(log(interpolated_curve$qxl)),max(log(interpolated_curve$qxk))), lty=2, lwd=1.25, col="red")
```

## Example 2 

Lets consider the following mortality surface with of a dimension 17 by 38, the Algerian mortality surface from 1977 to 2014 for males.

```{r echo=TRUE}
Bxt <- LT[2:18, 2:39]
head(Bxt[  ,1:5 ])
tail(Bxt[  ,1:5 ])
```

The surface of `log(nQxt)` can be plotted using 

```{r include=TRUE, fig.align="center" }
persp(x=c(0,1, seq(5,75, by=5)), y=c(1977:2014), z=log(Bxt), theta=-10, expand = 0.8 , phi=25, xlab="age", ylab="year", zlab="log (nQxt)", main="Figure 4. Five ages mortality surface, Algeria, males, 1977-2014")
```

Then lets deduce the single ages mortality surface using `getqxt()`

```{r include=TRUE }
interpolated_surface <- getqxt(Qxt=Bxt, nag=17, t=38)

str(interpolated_surface)
```

```{r include=TRUE ,fig.align="center"}
persp(x=c(0:79), y=c(1977:2014), z=log(interpolated_surface$qxt), theta=-10, expand = 0.8 , phi=25, xlab="age", ylab="year", zlab="log (qxt)", main="Figure 5. Detailled ages mortality surface, Algeria, males, 1977-2014")
```

