# Using the Bayesian Block Routine

!!! NOTE

    This package is designed to work exclusively with NuSTAR observations.

Here we are going to run though an example, showing how to run observation 30801024002 through the Bayesian Block routine, and how to get the files you need.

## Getting the data you need

To run this package, you will need four types of data files:

1. A barycenter corrected event file
2. Event files of the observation from FPMA, FPMB, or both (not barycenter corrected)
3. Light curve correction files for the event files used above.

This tutorial assumes that you have already obtained cleaned data products through NuPipeline and barycentered events and light curve correction files through NuProducts. These files are a standard part of the NuSTAR data pipeline procedure, and more information can be found in the [NuSTAR Manual](https://heasarc.gsfc.nasa.gov/docs/nustar/analysis/nustar_swguide.pdf).

### Event Files

Once you have cleaned data, we can extract event files. For accurate analysis, the event files must exclusively contain data from the target of analysis. To achieve this, we can use another HEASOFT tool called XSELECT.

The way we use XSELECT here is to essentially take a spatial cut out of data in the region we are interested in. For example, we can take a 30" source region centered on Sgr A \* by opening the XSELECT dialogue in the command line:

```
xselect 
```
XSELECT will then launch and you will be prompted to enter a session name. Lets call our session tutorial:

```
 
                         **  XSELECT V2.5b  **
 

>Enter session name >[xsel30916] tutorial
```

We will then enter the event data file that we want to extract information from:

```
read event nu30801024002A01_cl.evt
```

you will then be prompted to enter the event file directory. After you input the directory, it will ask to reset the mission to NUSTAR. Select yes. It will then output something like:

```
 
 Notes: XSELECT set up for      NUSTAR
 Time keyword is TIME       in units of s
 Default timing binsize =   5.0000
 
Setting...
 Image  keywords   = X          Y           with binning =    1
 WMAP   keywords   = X          Y           with binning =    1
 Energy keyword   = PI                     with binning =    1
 
Getting Min and Max for Energy Column...
Got min and max for PI:     0   4095
 
could not get minimum time resolution of the data read
MJDREF =  5.5197000766019E+04 with TIMESYS = TDB
 Number of files read in:            1
 
******************** Observation Catalogue ********************
 
Data Directory is: /data/directory/
HK Directory is: /data/directory/
 
 
        OBJECT      OBS_ID      DATE-OBS
      1 SgrAstar    30801024002 2023-04-13T04:49:49
 

```
We can then input our region size by giving the path and filename of the .reg region file that we want to select:

```
tutorial:NUSTAR-FPMA > filter region SgrA_30ac.reg
```

Then we extract the event:

```
tutorial:NUSTAR-FPMA > extract event
```

Then we can save the event:

```
save event
```

And give the output file name when prompted.

When asked to use fultered events as input data file, select no.
You must then quit the program before repeating the process for other data. You can quit the program by:

```
quit
```

and saving the session if you wish (not required).


## Running the data

After you have obtained the required files, you can run it through the bayesian block routine. Do do this, you must edit some lines in the `main.py` script. 

For example, lets run through an example using observation ID 30801024002. This observation was taken in 2023, so I would set:

```
year = "2023"
obsID = "30801024002"
path = ("path/to/prepared/files")

# Define output directory
output_dir = path + "/desired_output/"

[.....]

# Define input file paths
barycorr_event = path + "/nu" + obsID + "A01_cl_barycorr.evt"
event_file_a = path + "/nu" + obsID + "A01_xselected.evt"
event_file_b = path + "/nu" + obsID + "B01_xselected.evt"
lccorrfileA = path + "/SgrA_correct_50ac_fpmA_lcsrccorrfile.fits"
lccorrfileB = path + "/SgrA_correct_50ac_fpmB_lcsrccorrfile.fits"

# --- Define Parameters / General Used Ones ---
energy_min = 3.0
energy_max = 30.0

```
With the correct names for the event files, light curve correction files, and desired energy range for analysis. 

!!! tip "Using data from only one module"
    If you are only using data from one observation module, you can comment out the filenames of the ones that you are not using. For example, if you were only using data from module A, your script might look like:
    ```
    year = "2023"
    obsID = "30801024002"
    path = ("path/to/prepared/files")
    # Define output directory
    output_dir = path + "/desired_output/"
 
    [.....]
 
    # Define input file paths
    barycorr_event = path + "/nu" + obsID + "A01_cl_barycorr.evt"
    event_file_a = path + "/nu" + obsID + "A01_xselected.evt"
    #event_file_b = path + "/nu" + obsID + "B01_xselected.evt"
    lccorrfileA = path + "/SgrA_correct_50ac_fpmA_lcsrccorrfile.fits"
    #lccorrfileB = path + "/SgrA_correct_50ac_fpmB_lcsrccorrfile.fits"

    # --- Define Parameters / General Used Ones ---
    energy_min = 3.0
    energy_max = 30.0
    ```

Output csv files, plots, and txt files will then be put in the output directory specified above. 
