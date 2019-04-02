# Plans

## Research plans

What are the next obvious steps:

basic approach: justify subsequent innovations on the basis of new
data/stimuli to incorporate into the model

even more basic: examine new stimuli/conditions the model can already handle

another important question: what would improve the pilot data for my grant?

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
