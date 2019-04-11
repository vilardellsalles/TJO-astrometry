from glob import iglob
from astropy.io import fits

for imname in iglob('*imc.fits'):
    image = fits.open(imname)[0]
    image.scale('float32')
    image.writeto(imname, clobber=True)
