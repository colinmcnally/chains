# chains

This is a prototype built to try out some thing about managing and analyzing scans
of n-body simulations. Singularity containers are used to provide a repeatable and
portable environment for building and executing code. The idea is that this is 
particularly handy when dealing with myriad dependencies, like Python based 
analysis scripts tend to end up having.

The basic probem this system handles is the long-term stability of resonant 
chains of planets, in the style of:

> Matsumoto, Y., Nagasawa, M., & Ida, S. (2012). 
> The orbital stability of planets trapped in the first-order mean-motion resonances. 
> Icarus, 221(2), 624–631. http://doi.org/10.1016/j.icarus.2012.08.032

The first set of results generated confirms the basic pattern in
Matsumoto et al. (2012), that although chains in mean motion resonance (MMR) have 
much longer stability times than nonresonant ones with similar spacing, there exists 
a critical length of chain beyond which the stability timescale rapidly decreases.
For planets of planet mass / star mass ratio q=1e-5 they give critcal chain lengths 
of 8 in 6:5 MMR and 4 in 7:6 MMR.

The data from these simulations on the orbit crossing time (proxy for stability time)
is:

![plot of collision times for chains of planets](./plots/tcross_q1em5_n3-10_p3-6.png)

where orange triangles indicate lower bounds.
