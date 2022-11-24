# Execute the files
Once you downloaded the EXECUTION.Beta folder, open the terminal and execute the file execution.py
Is not needed to change any path to execute and get the output of any script.

If when you are executing the scripts you find the following error: 
Permission denied: '/home/USERNAME/Downloads/EXECUTION.Beta/Files/soft/NACCESS/naccess'

Then type the following in your terminal:
$ cd Files/soft/NACCESS/
$ csh install.scr

If you now go back to the EXECUTION.Beta folder and execute again the execute.py file you should not have more errors.

Is recomended to execute the files in cluster because two of them take a lot of time if you execute it in your computer.
In the output you should get a pdb fixed protein file, a plot for the alanine scanning and a txt with the energy outputs.
