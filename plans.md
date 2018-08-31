

# What we've learned

I can probably also greatly reduce the total search by 
first sampling 2-3 runs across all parameters and only persueing
those that have at least some sign of bistability

simulations plans:

 - organize parameters
 - run W_m_σ_t, and W_m_σ_ϕ simulation (~10 hours)
 - run freq,scale,track simulation, 3st only, 3 runs (~130 hours)
 - then run hopefully a much smaller set over the parameters that showed
   reasonable bistability

 - we're going to try to switch to 0.7 before too many other changes

Next plan: move to 0.7, since MARCC is now 0.7 only this may also make stuff
even faster, if I'm lucky
  - as the first step to that, I'm moving AuditoryModel to a separate 0.7
    package. That is a small amount of code and should be a good test bed to get
    a sense for how much work it is going to be and whether it's really possible
    (the only other big question mark is the NMF code, which *seems* to be
    working on 0.7 but there might be a sublte bug given the issues with the
    test...  perhaps that will be fixed before I have to worry about it though).

we also want to fix the scales to make sure we have the right values (might have
to alter the parameters I actually use, so I have to be careful with this
change).
