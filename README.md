# chains

This is a prototype built to try out some thing about managing and analyzing scans
 of n-body simulations. Singulairty containers are used to provide a repeatable
 and portable environment for building and executing code. The idea is that 
this is particularly handy when dealing with myriad dependencies, like Python 
based analysis scripts tend to end up having.

The basic probem this system handles is the long-term stability of resonant 
chains of planets, in the style of:

> Matsumoto, Y., Nagasawa, M., & Ida, S. (2012). 
> The orbital stability of planets trapped in the first-order mean-motion resonances. 
> Icarus, 221(2), 624–631. http://doi.org/10.1016/j.icarus.2012.08.032

