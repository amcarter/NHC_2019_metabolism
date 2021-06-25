Readme for NHC 50 years later coding workflow                         |

# Calculate physical stream properties
  - Depth, discharge, velocity, K600

 1. gap fill missing NHC and UNHC level data from the other sensor, 
	snap filled data to match

 2. Build a discharge rating curve based on measurements
   
 3. Correct rating curve to deal with no flow days
 
 4. Interpolate discharge at all sites based on accumulated watershed 
  area. (see Leach et al 2017) Calculate watershed area based on NHD?

 5. Build a relationship between channel thalweg depth and avg depth
     Data: geomorphology survey data from NHC/UNHC done by E. Moore

 6. Fit average depth by discharge relationship based on 
  Leopold & Maddock 1953 (D = c*Q ^f)
 	D is depth, Q is flow, c is depth at unit discharge

 7. Use average depths and widths and discharge to calculate velocity
 
 8. Get slope from NHD? from Amanda's thing?

 9. Calculate K600 according to Churchill 1952 and Raymond 2012


# Calculate Metabolism
# using stream Metabolizer
# using Hall 1970 method

 