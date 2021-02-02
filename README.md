COVID - Ventilation
===================================
Ventilation simulation in a room with multiple active scalar fields
(temperature, humidity, CO2). Inserts bodies using point-in-convex hull
algorithms.

References
 * [obj_io: John Burkardt @ University of South Carolina](https://people.math.sc.edu/Burkardt/f_src/obj_io/obj_io.html)
 * [geometry: John Burkardt @ University of South Carolina](https://people.math.sc.edu/Burkardt/f_src/geometry/geometry.html)

Code development by Chong Shen Ng, Steven Chong, and Rui Yang and others.


Note
===================================
We assume the breath out volume of 0.5L, carrying injectmeanq for the whole volume of air

Hard coded part:
For addbodyobj.f90, injectedvol,injectmeanq,time_shift,breath_interval
are normalized by system height 3m, free fall vel 0.71m/s, and free fall time 4.25s
