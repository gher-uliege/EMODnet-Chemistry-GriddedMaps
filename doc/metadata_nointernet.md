## Add the netCDF metadata when no internet connection is available

This is the typical situation when your code is run from a cluster where the nodes are not connected to internet.        
In this case, the metadata have to be prepared beforehand, from a desktop computer for example.

### Instructions
1. Define the metadata, in the form of an ordered dictonary.
2. Define the file name, the variable and the domain. 
3. Get the attributes from the Vocab server (using a machine that has an internet connection)
4. Save the information in a text file
5. Upload the text file to the cluster where the computations will be run.


### Example

A full example can be found here:      
https://github.com/gher-uliege/Diva-Workshops/blob/master/notebooks/4-Postprocessing/save_attributes_file.ipynb
