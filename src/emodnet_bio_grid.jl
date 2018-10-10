
# bounds for the Baltic Sea
gridlon = 9.0 : 0.1 : 30.8
gridlat = 53.0 : 0.1 : 66.1

datadir = get(ENV,"DATADIR",expanduser("~/tmp/Emodnet-Bio"))

years = 2007:2013

# bounds for the Benthos product (Atlantic)
dlon, dlat = 1/10., 1/10.
gridlonBenthos = -10. : dlon : 35.
gridlatBenthos = 36. : dlat : 73.
