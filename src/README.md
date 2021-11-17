## Some notes on what each module here does

#### Configfile

For reading in text files specifying settings in `keyword = value` format. We use this for reading in some of the data.

#### clusters

Handles the cluster **modeling** aspects of the code, specifically gas profiles, the NFW mass model, and the actual fgas model whose parameters will be marginalized over. Constants: m_He/m_p.

#### fgas

Contains the main data organizing structure, `Cluster`. Also covers nuisance parameters specific to the measurements themselves. Includes the function that predicts fgas measurements based on parameters in `fgas` and `clusters`.

#### fwrapper

Fortran interface to `wrapper`; see below.

#### lensing

Data structure and likelihood calculation for weak lensing measurements. Constants: c^2/G.

#### util

pi and some inline integer power functions.

#### wrapper

Top-level interface for cosmology codes. Functions to initialize, load data, evaluate the overall likelihood, etc. C++ and C versions. Look here for vague intructions on how to call the code from elsewhere.

#### xray

Data handling for X-ray measurements, and related unit conversions. Constants: G, m_p, Mpc.
