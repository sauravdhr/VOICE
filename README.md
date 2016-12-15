# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###
* VOICE: Viral Outbreak InferenCE
* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

### Project Structure ###
* Normalization

* Simulation for all pairs
run_simulation_for_all_pairs.py
The script runs pair simulations for all pairs of all hosts specified in a given folder.

* Graph coloring
color.py
The script draws graphs for simulations results for one pair.

Input parameters:

simulation result folder with .out files;
first host graph template (.dot file);
second host graph template (.dot file).

Graphs are currently saved to simulation folder.

* Analyses Simulation Statistics
calculate_simulation_stats.py

Input parameters:

Simulations resulting folder (e.g. out/graphs).

Result text files are saved to statistics folder by default.

