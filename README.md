

# Bistability

# RESET!!!

okay, now we're going to think about how this works with the cortical model,
rather than the DNN

things I need to be working on (by priority)
X. plot response of L1 to varying delta F
2. bistability in layer 1 (i.e. via appropriate inhibition patterns in layer 1)
3. visualization of layer 2 and 3
X. fundamentals of bistable models, expanding to higher dimensions
X. figure out what the heck is up with sigmoid screwing up
   model behavior for L1
6. what about sigmoid messing up L2 and L3?
   - note: does it help if we use 0 instead of -1 in L3
  
# File organization

* **data**: all files generated by model simulations
* **models**: all of the source code to run model simulations
  * **rankin2015**: Replication of Rankin et al. 2015
  * **chakrabarty2017**: Replication of Deb's model results.
    * **model**: Auditory-model related code
	* **stimulus**: Stimulus-generation related code
	* **original**: some of the original scripts used by Deb, for reference
  * **little**: my current version of the model
* **plots**: all of the plots generated by model simulations
* **util**: utility files for managing the project (see below)

# Questsions for Deb

1. what is up with intertia_delta? why modify only the first 4 time points? (to emphasize the target sound, which helps separate it from the non-target)
2. why subtract the starting value for the buildup experiemnt,
   i know it gives the "right" answers but what is the theoretical
   justification for it? (short answer, there isn't really one, deb just wanted to get a result consistent with the human data)

# Problems/Questions 
1. I noticed that when you plot out the buildup by frequency, there are some
   strange re-orderings of the delta's at the most extreme frequencies, using my
   continuous response setup. Why? NOTE: This happens a little bit in the upper
   frequencies for deb's approach as well, but not nearly as badly.

# Older thoughts

This is a work-in-progress bistability modeling project. 

- [x] Step 1 - replicate the basic behavior described in [Rankin et al. 2015](https://doi.org/10.1371/journal.pcbi.1004555)
- [ ] Step 2 - replicate outcomes of Deb's (Chakrabarty) model
  - [x] asynchronous vs. synchronous
  - [x] buildup
  - [ ] [Micheyl et al 2013](https://doi.org/10.1121/1.4789866) data
        - [ ] get asynch working with new model
  - [ ] generate and save pdfs of results
        - [x] asynch
        - [x] buildup
        - [ ] micheyl
- [ ] Step 3 - introduce simple bi-stable components to Deb's model
  - [x] figure out how to generate a continuous response 
        from Deb's model
  - [x] refactor the model to ensure I can modify it with ease
  - [x] create a graph for the continuous response version of the build up experiment
        - [x] include a version that breaks down the results by frequency
        - [x] do the same for the old approach
        - [x] compare new approach using both ``L_2`` and ``L_1`` norm
  - [x] create a graph to visualize stimulus response in layer1
        - [x] simple - discriminability of a and b
		- [x] map MI frequency response of each unit, order by frequency and bandwidth
  - [x] introduce adaptation in layer 1
        - [x] how does first approach affect psychophysical benchmarks
              - [x] generate data for async (with and without adaptation)
              - [x] generate graph for async
              - [x] generate data for buildup (with and without adaptation)
              - [x] generate graph for buildup
        - [x] generate graph of adaptation effect on units
  - [x] introduce mutual inhibition in layer 1
        - [x] implement
        - [x] create plots
  - [x] introduce adaptation and MI to layer 2
        - [x] find a usuful way to plot
        - [x] blindingly implement and observe the results

Hold on: is it worth even generating these visualizations? I it would
be useful to determine if there is a way to actually get the right
behavior for this network.

Hypotheis: Bistability will not work in layer1 and layer2. In the rankin model,
there is a monotonic relationship between the presence of a percept and the
magnitude of values in the network. The dynamics of the bistability rely on
this.

I can test this by creating two simple models. One with the zero-means-something network, and one with aribtrary transformations to a vector space. This seems like the most useful test. 

In this respect bi-stability *should* be ablqe to work on layer 3.

if I can generate a simple AAE for layer1 that is more interpretable
that could also be helpful... This would at least help determine if this
approach is worth pursuing.

things to show/discuss with mounya

generate some graphs for three versions of a bistable model
- essential, 2 element network
- course coded, 10 element network
- distributed representation, 10 element network

waiting on deb's info for visualization

propose using a model with ELU's and sparisity prior (AAN's)

the key thing wee need is for units to some how "represent" a particular
grouping, and for those to compete with another particular grouping

<!--------------------------------------------------------------------------->

things to discuss with mounya:

topic 0 - summary of results, no bi-stability produced with
first choice of parameters? how are we going to pick parameters?

first appraoch - visualize to develop udnerstanding and manually adjust
learning? - something about natural scenes?

topic 1 - what do we need for bi-stability

grouping - I've asked, in principle what ingredients do we need for bi-stability
when we have more units: one way to make this work is to group
the inhibitions properly - this suggests we need layer3 to either
be the source of bi-stability OR for it to have top-down effects on the
lower-level units: alternatively there are potential ways for bi-stability
to occur on the subj-object level. Even in Cao et al (2016) we have to
"know" which units compete with one another.

sparsity/montonicity - I think another thing we need for bi-stability to "work"
is that utlimately inhibition is driven by higher activity in one region, so,
higher activity has to correspond to percept presence, which doesn't seem to be
the case in the current model.

topic 2 - visualization

interpretable models - I've also looked at a few ways we might visualize the units:
many of the approaches I've found would require re-training the model
to learn more interpretable representations: AAN's with ELUs and a sparsity prior,
or InfoGAN. Both have working implementations in TensorFlow, so they'd be "easy" to adopt: the challenges wouldbe application specific.

using such a model could kill two birds with one stone - visualization challenges
and the viability of bi-stability 

local, linear aproximations - however I could use something like LIME (or
similar) to examine what features are locally being used to determine
"discrimination" of a and b events at the different layers, and then examine how
these dynamics change as a function of the new parameters I'm introducing. I
could do this for subsets of neurons, at various layers, e.g. by cluster. I
could probably pick various parameters of the model to visualize on top of this
"pixel"-wise measure.

Next questions:

can we slow down buildup? (do certain frequencies)
