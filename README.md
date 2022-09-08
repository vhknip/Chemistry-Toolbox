# About the Project
This is a personal collection of useful functions for data analysis in chemistry. This includes mass spectra match factor calculations, mass spectra visualizations, and common sample grouping functions for consolidates result tables.

## Features
* Mass Spectra Comparisons Calculations
  * Direct matching: normal match factor found in spectral comparison software
  * Forward matching: match factor calculation evaluated only on m/z values detected in the sample
    * Particularly useful with lower resolution spectra or low concenstration analyte.
  * Reverse matching: match factor caluclation evluated only on m/z values detected in the refereence
    * Generally useful for high-noise spectra
 
* Spectra Comparision Visualizer
  * Generates interactive mass spectrum comparison plot with match factor information displayed
  * Option to export as html file to share with others for reports or compound rationalization
  * Plots are generated in plotly.graph_objects
  * Allows for spectra comparison generally only available on instrument software to be done anywhere
  
* Sample String Analysis
  * Searches through a sample name for common names / abbreviations of an attribute and assigns the correct valule of that attribute for thar sample.
    * Useful for analyses where operators need to group samples by attributes with information included in the sample name before doing and calculations.
  * Finds longest common prefix in a collection of strings
    * great for grouping different replicates of the same sample using some enumeration feature in instrumentation software
  * Converts string representation of a list to a list
    * Mass Spectra and other chromatographic data points are exported and formatted as lists.
      * When read in to a dataframe, these are taken in as lists, so conversion before analysis is required.
      
* To-Do
  * Add jupyter notebooks showing examples of how functions work
    

