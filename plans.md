# Plans

## Research plans

What are the next obvious steps:

basic approach: justify subsequent innovations on the basis of new
data/stimuli to incorporate into the model

even more basic: examine new stimuli/conditions the model can already handle

another important question: what would improve the pilot data for my grant?
(largely this will involve extending the stimuli since the grant is going to
validate the model across different paradigms, but it could also be worth
actually collecting some pilot data using the paradigms that work weel)

### First steps

Tasks

- [ ] clean up existing stimuli (still a few things to work through)
- [ ] organize the different tasks, figure out what changes they require
- [ ] create a list of priorities and things to try

One ideal way to go about this: identify novel stimuli that appear to play on
the same challenges other paradigms of ABA play on. This way we can expand
the scope of ABA paradigms covered by the model *and* cover more types of
stimuli.

Based on our meeting:

we want to examine how many stimuli we can cover with an extended model

These are my notes about different sets of stimuli we're considering:

#### ABC stimulus

([Thomassen & Bendixin 2019](https://asa.scitation.org/doi/10.1121/1.4973806)). 
There are two challenges here,

1. We'll need a much larger number of scene interpretations; will slow down
inference by 27/4 = 6.75
2. It will be much more challenging to generate a decision for this stimulus; one approach could be a sort of template matching --- compare the present mask to a set of template masks representing different possibilities.

### Interleaved Speakers

([Gandras et al 2018](https://doi.org/10.1016/j.neuroscience.2017.07.021)) - A is a male speaker, B is a female spaker, and syllabels are interleaved.

This might actually be fairly straightforward to model. The key challenge would be 
making sure the features modeled were relatively reasonable.

### Varying temporal rates

(Van Noorden 1975): It's not yet clear to me whether this would work or not.
It's definitely worth testing; I'm guessing this would be a challenge
however, in part because the time constants are a poorly explored component
of the model.

I think what's needed here some inference over the temporal blurring (in addition to spectral blurring) of the object level analysis.

### Shepard tones

Ala Chambers '17, and others --- This could invovle just adding some short-term
temporal modulation features. Likely an increase in complexity of x10.

This could handle the basic ambiguous stimulus, but wouldn't get context. We'd need to model build-up to get that right (but it might be worth doing that)

There's also the stuff by Grave's --- look into that...


How would decision making work here?

I think this could be pretty simple; it would be about some readout of a mask
of the pitch trajectory (and or, the read out of pitch-change features)

### Ambiguous Brightness 

Ala Seidenberg '18 --- This would presumably be fairly similar to the work
for the Shepard tones.

### Context stimuli

These are Joel's stimuli --- similar challenge to the Chambers '17 data,
in terms of getting build-up to work.

### Tone disruption

These are my and Nathan's data --- also related is the data from Rankin '17.
This would likely invovle modeling some form of attention/buildup of some
kind; that is a layer following the object level analysis.

### Verbal Transformations

These are, in concept pretty straightforward: the big challenge here would be
decision making. We would also likely need additional features.

## Older comments

### things I might want the model to handle:

- individual data for aba task? (maybe already)
- fill out the van-noorden diagram: handle ambiguity in temporal rates (unknown, but probably not working)
- bulid-up & context (seems like the most obvious place to improve)
    - many other aba manipulations depend on build-up, so that would be a good
      place to start
- the pure-tones organized by speech stimuli from stephen david's lab (maybe already)
- handle shepard tones (not yet, could be pretty straightforward, might come along with build-up/context)
- something about the role of attention? (relate to top-down effects)
    incorporate auditory flashes

one basic approach: incorporate a number of these different stimuli,
and examine how the improvements below change behavior, and use this as a guide
for what stimuli to focus on

things I might do to help the model handle more:

in some ways it would be more interesting if I could handle a broader range
of paradigms than a broader range of stimuli at the moment; seems like there
is lower hanging fruit there; it would also probably be easier to 

### broad ways to change the model:

1. more intra-layer communication
2. more sophisticated adapt/inhibit/noise
3. tighter integration of adapt/inhibit/noise and computations
4. more realistic features

#### top-down stuff

this stuff might address the first three of these ways of changing:
- extend object level to include top-down effects
- perhaps follow brascamp's framework to create layer-to-layer communication,
    possibly replace current adapt/inhibit/noise implementation?

#### better features stuff

this stuff might address the fourth way of changing the model:

- more features:
    - pitch
    - better use of temporal rate
    - a second "Central" layer?

I'm inclined to first try a data driven approach, as it might handle many
more types of stimuli before adding a lot more features.

- more data driven object-level (e.g. some sort of NN trained to predict)
- data driven other layers - beta-VAE
    - incorporate adapt, inhib, noise into computations? of data-driven layers?
    - if I went this route might make sense to consider just a single trained layer to start with
    could replace earlier layers later on

