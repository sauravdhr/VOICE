### VOICE: Viral Outbreak InferenCE ###
* This tool determines the relatedness between viral samples if they consist of haplotypes.

### Demo Run ###
Change absolute paths in the file test_data/test_outbreak_tasks
Run command:
python main.py -t test_data/test_outbreak_tasks

### Input Format for Different Tasks ###
* Running simulation
main.py has a mandatory argument -t which is responsible for a task file.
Every single line in the task file is a separate task. Task can be one of two types:
determining relatednes between a pair of samples, and determining relatednes within an outbreak.
The line in the task file should consists of two or three parameters separated with a space,
and it looks like following lines:

p path/to/a.fas path/to/b.fas
d path/to/outbreak/fasta/files

The first parameter is just a letter 'p' or 'd'. 'p' is responsible for the task for pair of samples,
while 'd' is responsible for determining relatedness within outbreak.
Samples should be represented as fasta files. Every fasta file should end with '.fas' extension and
it should contain all known haplotypes from the same viral sample.
'p' mode works with two fasta files which should be specified in a task as the second and the third parameter.
'd' mode works with multiple fasta files allocated at the same directory, it will run multiple pair simulations
for all possible pairs in the folder.

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

