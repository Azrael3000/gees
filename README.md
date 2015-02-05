GEES
====

GPL Euler equation solver

-------------------------------

GEES (pronounced like cheese because it's yummy) solves the Euler equations in 1-D. It is intended as teaching tool and to produce reference solutions and the whole code is in Fortran90, although it should be understandable for people who don't know that language very well.

-------------------------------

**Implementation details:**

- Time-stepping: 4th order explicit Runge Kutta
- Flux calculation: MUSCL scheme (Kurganov and Tadmor central scheme)
- Flux limiter: Ospre
- Boundary conditions: Second order polynomial extrapolation
- Equation of state: Tait

-------------------------------

There is also a LaTeX document to be found in the *doc* folder that shows the gory mathematical details.

-------------------------------

Before putting the code on github it was (and still is) available on my [homepage](http://sci.amconception.de/index.php?nav=gees). I put the code up for review on [CodeReview](http://codereview.stackexchange.com/questions/10326).
