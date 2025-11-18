# The Use of Bayesian Block Representations with NuSTAR Observations

This package was developed to apply Bayesian Block statistical methods to NuSTAR X-ray data of Sgr A\* with the goal of determining "flaring" periods where the count rate of Sgr \* is above the baseline quiesent rate [Zhang 2017](10.3847/1538-4357/aa74e8).

The Bayesian Block statistical method [Scargle 2013](10.1088/0004-637X/764/2/167) detects statistically signifigant change points in time series data and was developed with astronomical applications in mind. The idea is to optimally segment data from a target, apply likelihood functions, and surpress observational errors.

Most Baysian Block analysis applies to three types of data measurements: time-tagged events, binned events, and measuremnts of data at various time intervals (e.g. gamma ray observations). These three types of data have different likelihood statitics, being, generally speaking, Bernoulli, Poisson, and Gaussian. 

NuSTAR data can be best described by time-tagged event data, as the telescope reports photon arrival times. Segmentation for time-tagged works by creating time bins of photon arrival times. The width of the bin is defined by the time since the last photon arrival and the time until the next photon arrival. To account for observational errors, which mostly include telescope effects such as detection effeciency, each time bin is divided by a correction factor that has a value between 0 and 1, with 1 being perfect detector effeciency. Luckily, the NuSTAR telescope already reports the detection effeciency per photon hit, so we don't have to do any calculations to determine this. 

A likelihood function is applied to track the probability of how many photons should arrive in a given time before it is considered to be a change point in counts. This is done iteratively until another change point is detected. These change points then mark the start and stop time of "blocks" that have similar count rates. Since the purpose of my project is looking for periods of flaring from Sgr A\*, blocks with higher count rates than Sgr A* quiescence is labeled a flare. 

A data chalenge for applying a Bayesian Block algorithm to NuSTAR data in particular is the observation gaps in the data and the fact that NuSTAR has two observation modules. Observational gaps in NuSTAR data is due to the telescope shutting down as it passes over the South Atlantic Anomaly to protect the instruments on board. Because of this, NuSTAR reports "Good Time Intervals" (GTIs) of good and usable data. Time between these GTIS must be removed before running the data through the algorithm, or it will detect every gap as a change point. For accurate time reporting, the observational gaps are inserted back into the data so that the times of change points reflect actual time. 

To account for the fact that NuSTAR has two observational modules that simultanously look at a target, the data must be combined. That combination happens at many stages. First common GTIs must be found between the two modules, and then photon arrival times and correction factors are combined into one list to go through the bayesian block algorithm.




