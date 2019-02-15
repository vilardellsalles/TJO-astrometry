Required commands to clean and reduce images:

```
icat gaia-imgqs <images>
icat gaia-reduce *
icat light-curve *_catalog/*.dat
```

Issues found:
* No flats available for V band filter at the proper focus position => use flats with other focus positions
* Observations were obtained with two different instruments (MEIA and MEIA2), use proper calibration images for each instrument