can we combine an short-term acousting memory (something like the approach
described in McDermott et al 2011 (PNAS)) within this model, could that take on the 
role of 'attention'? A.k.a. enhancing a particular object representation. 

NOTE: the TC model is really a version of some sort of *simple* WM, i.e. each
component is an 'object' in this memory. Those components could be used
to identify peaks and use that to build some ongoing template of the 
principle object.

Implementation questiosn for MI

- do we determine inhibition as a function of overall input
  to layer or based on some estimate of the strength of
  a scale or rate relative to the others? (i.e. does
  overall activation of a particular scale matter or just that
  a scale is present?) is this well formulated, what
  would the "estimate" really be?
  
- do we use rate and/or scale MI?

- is there adapt/MI present in the spectrum?

- does/how does the object-level mask guide this MI (as proposed)?

- do we really want to use the mangitude, or would 
  accounting for phase be appropriate?
