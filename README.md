
# Bistability

This is a work-in-progress bistability modeling project. 

- [x] Step 1 - replicate the basic behavior described in [Rankin et al. 2015](https://doi.org/10.1371/journal.pcbi.1004555)
- [ ] Step 2 - replicate outcomes of Deb's (Chakrabarty) model
  - [x] asynchronous vs. synchronous
  - [x] buildup
  - [ ] [Micheyl et al 2013](https://doi.org/10.1121/1.4789866) data
  - [ ] generate and save pdfs of results
- [ ] Step 3 - introduce simple bi-stable components to Deb's model
  - [ ] figure out how to generate a continuous response 
        from Deb's model
  - [ ] refactor the model to ensure I can modify it with ease

# Setup

If you modify files that end in `.inputs.ipynb`, they should be stripped of
their output upon commit. Locally your output will remain unchanged unless you
use git checkout.

To set this up, you need to setup local configuration and enable execution privlegeds for the output filtering script:

```sh
$ git config --local filter.dropoutput_ipynb.clean ./util/ipynb_output_filter.py
$ git config --local filter.dropoutput_ipynb.smudge cat
$ chmod +x util/ipynb_output_filter.py
```

# File organization

* **data**: all files generated by model simulations
* **models**: all of the source code to run model simulations
  * **rankin2015**: Replication of Rankin et al. 2015
  * **chakrabarty2017**: Replication of Deb's model results.
    * **model**: Auditory-model related code
	* **stimulus**: Stimulus-generation related code
	* **original**: some of the original scripts used by Deb, for reference
* **plots**: all of the plots generated by model simulations
* **util**: utility files for managing the project (see above)
* **plans.md**: Notes about the plans for this project
