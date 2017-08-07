# Bistability

This is a work-in-progress bistability modeling project. Step 1 is to replicate the basic behavior described in [Rinzel et al. 2015](https://doi.org/10.1371/journal.pcbi.1004555)

# To avoid generated output being saved to git

You can have all files that end in inputs.ipynb stripped of their output, upon
commit. Locally your output will remain unchanged unless you use checkout.

For this to work, you must setup the local configuration and enable execution privlegeds for the filtering script:

```sh
git config --local filter.dropoutput_ipynb.clean ./ipynb_output_filter.py
git config --local filter.dropoutput_ipynb.smudge cat
chmod +x ipynb_output_filter.py
```
