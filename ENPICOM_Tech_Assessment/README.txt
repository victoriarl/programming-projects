ENPICOM Tech assignment answers - Victoria Ramírez López
26/08/2022

Question 1: 

I used the Levenshtein distance, defined as: "The minimum number of single-character edits 
(insertions, deletions or substitutions) required to change one word into the other."

With this, I calculate how many edits separate a from each sequence in B and I select that
with the shortest distance (i.e. amount of edits) as the 'best matching sequence' b. 

The file must be run with the following command line:
python3 question1.py -a a.txt -B B.fa
	-a = .txt file containing the target sequence 'a' 
	-B = FASTA file containing the database 'B' 
	
NOTE: The script has not been optimised. Specially the calculation of the Levenshtein distance.
I am aware that there are oportunities to optimise the for loops that I am utilising in this
implementation, but in order to not extend the amount of time that I spent on the assignment 
more than it was expected I have decided to submit this first draft implementation and I'm 
open to discuss its optimisation when we go over it in detail. 

Question 2:

Based on the description of the challenge in Rosalind, I found the bug to be that the script
was not checking whether the incorrect sequences were at a Hamming distance of 1 with respect 
to EXACTLY ONE correct read in the dataset (or its reverse complement). I implemented this 
check in the 'find_corrections' function and added inline comments to explain my reasoning. 



 
