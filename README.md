The data for this project can be summarized in the following table:


| object     | instrume | filter | Number of images |
| ---------- | -------- | ------ | -------- |
| 12 Oph     | MEIA     | B      |       56 |
| 12 Oph Ref | MEIA     | B      |        3 |
| HD 135599  | MEIA     | B      |       78 |
| HD 180617  | MEIA     | B      |       19 |
| HD 98712   | MEIA     | B      |       78 |
| OT Ser     | MEIA2    | V      |      534 |
| V2306 Oph  | MEIA2    | V      |      150 |

The following commands are required to clean and reduce the images:

```
icat gaia-imgqs <images>
icat gaia-reduce *
icat light-curve *_catalog/*.dat
```

However, the following issues were found:
* No flats available for V band filter at the proper focus position => use flats with other focus positions
* Observations were obtained with two different instruments (MEIA and MEIA2), use proper calibration images for each instrument