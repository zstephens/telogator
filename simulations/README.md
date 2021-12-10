## Scripts for simulating data

The scripts in this directory were used to generate the simulated data results in our manuscript.

`python3 simulate_reads.py \ `  
`    -o sim_dir/ \ `  
`    -c default.cfg \ `  

If you're interested to run these scripts on your machine you'll need to point the script towards an edited config file, which specifies the paths to certain executables, and which parameters you want to sweep over.

Accuracies can be generated with this command:

`python3 parse_simulated_results.py \ `  
`    -i sim_dir/ \ `  
`    -o accuracies.tsv \ `  

The resultant values can be used to make plots and what-have-you.
