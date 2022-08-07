The Java code can be compiled as follows:

javac dtl_algo/DTLReconciliation.java
jar cvfe ARTra.jar dtl_algo.DTLReconciliation *


Brief description of source code:
The Java code implements DTL Reconciliation and all the rule-based heuristic 
algorithms to classify inferred transfer nodes, and the Python code implements 
the machine learning classifier. The Java code runs DTL reconciliation and 
the rule-based classification heuristics and generates the input file used by 
the machine learning classifier (Python code). The Python Code is executed on 
this input file and it computes the machine learning based classification of 
transfer events. This classification is then returned to the Java code, which 
then prints the final output file. 

