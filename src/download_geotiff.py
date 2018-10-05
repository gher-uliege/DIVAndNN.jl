# download SeaBed Habitat using WCS


from owslib.wcs import WebCoverageService


width=2000
height=2000

wcs = WebCoverageService('https://ows.emodnet-seabedhabitats.eu/wcs?',version='1.0.0')
print(wcs.contents)
print(wcs.contents.keys())

identifier = "emodnet:eusm_2016_baltic_substrate_raster"

tmax = wcs[identifier]
dir(tmax)
print(tmax.boundingBoxWGS84)
print(tmax.timepositions)
print(tmax.supportedFormats)


output=wcs.getCoverage(identifier=identifier,bbox=tmax.boundingBoxWGS84,format='GeoTIFF',crs="EPSG:4326", width=width, height=height);


f=open('foo3.tif','wb')
f.write(output.read())
f.close()
