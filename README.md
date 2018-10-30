# DIVAndNN

DIVAnd with neural network


```bash
git clone git@gitlab.com:gher-ulg/DIVAndNN.git
cd DIVAndNN/
```

See `REQUIRE` for necessary packages.

Run in `DIVAndNN/src`

```julia
include("emodnet_bio3.jl")
```

Data is assumed to be in the directory defined in the environment variable `DATADIR` (`$HOME/tmp/Emodnet-Bio` if this variable is not defined)

A copy of the example data is in `nic4:/home/ulg/gher/abarth/tmp/Emodnet-Bio`




## Zooplankton products

Analysis of 40 zooplankton species for the years 2007, 2008, 2010, 2011, 2012 and 2013 have been made. There is no data for the year 2009.

The method uses DIVAnd and a neural network with the following inputs:
* Distance from coast (from GSFC, NASA)
* Bathymetry (from EMODNET bathymetry)
* Salinity (from SeaDataNet)
* Temperature (from SeaDataNet)
* Dissolved oxygen (from EMODNET Chemistry)
* Chlorophyll concentration (MODIS-Aqua from NASA)

Additionally the position (latitude and longitude) and the year are provided as inputs to the neural network.

Abundance values (in the figures and NetCDF files) are expressed in number per m2 and transformed by the function log(x/a + 1) where a is 1 m-2. 


The full list of the species is:

* Acartia (Acanthacartia) bifilosa
* Acartia (Acanthacartia) tonsa
* Acartia (Acartiura) clausi
* Acartia (Acartiura) longiremis
* Amphibalanus improvisus
* Appendicularia
* Bivalvia
* Bosmina (Eubosmina) coregoni
* Bryozoa
* Calanus finmarchicus
* Centropages
* Centropages hamatus
* Cercopagis (Cercopagis) pengoi
* Cnidaria
* Cyclopoida
* Daphnia
* Daphnia cristata
* Echinodermata
* Eurytemora
* Evadne nordmanni
* Fritillaria borealis
* Gastropoda
* Harpacticoida
* Keratella cochlearis
* Keratella cruciformis
* Keratella eichwaldi
* Keratella quadrata
* Limnocalanus macrurus macrurus
* Mysidae
* Oikopleura
* Oithona
* Paracalanus
* Pleopis polyphemoides
* Podon intermedius
* Podon leuckartii
* Polychaeta
* Pseudocalanus
* Rotifera
* Synchaeta
* Temora longicornis
