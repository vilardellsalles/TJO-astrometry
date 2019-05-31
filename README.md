# Documentation

## Introduction

This repo contains the steps followed to reduce the data from the TJO proposal
93 entitled "First astrometric measurements of stellar activity with 
spectroscopy". The data for this project can be summarized in the following 
table:

<div align="center">

| object     | instrume | filter | focus range | images |
| ---------- | -------- | ------ | ----------- | ------ |
| 12 Oph     | MEIA     | B      | 4464 - 4474 |    56 |
| 12 Oph Ref | MEIA     | B      | 4466 - 4467 |     3 |
| HD 135599  | MEIA     | B      | 3935 - 4054 |    78 |
| HD 180617  | MEIA     | B      | 3217 - 3463 |    19 |
| HD 98712   | MEIA     | B      | 3440 - 4231 |    78 |
| OT Ser     | MEIA2    | V      | 7859 - 8884 |   534 |
| V2306 Oph  | MEIA2    | V      | 7865 - 8871 |   150 |

</div>

## Data issues

I have observed the following issues regarding data quality:

* **HD 98712.** The focus position for this object is not constant and several 
images were obtained with a different focus:

<div align="center">

| focus | images   |
| ----- | -------- |
|  3300 |        4 |
|  3450 |       62 |
|  3600 |        4 |
|  3900 |        4 |
|  4050 |        2 |
|  4200 |        2 |

</div>

* **OT Ser and V2306 Oph.** None of the V band images have flats with the
proper focus position. Therefore, flats of two posterior nights (2015-04-05 and
2015-04-07) were downloaded and used regardless of focus position.

* **12 Oph.** This target has no comparison stars in the field of view and,
therefore, it was supposed to have a reference exposure (12 Oph Ref) for each
target exposure. However, reference images were obtained for only one of the
two nights that 12 Oph was observed. In addition, the night with both target
and reference images only contains 6 images of the target (of the 56 obtained)
and 3 refence images.

Considering the lack of flats with the proper focus, I have finally decided to
use all the flats available, regardless of the focus position. 

## Reduction

The following commands are required to clean and reduce the images:

```
icat gaia-imgqs <images>
icat gaia-update --imgqs 0099 <images with invalid IMQS>
icat gaia-reduce *
icat light-curve *_catalog/*.dat
```

IMGQS is a keyword inserted in the FITS header after ```gaia-imgqs``` that 
contains several quality indicators. All the science images used are valid and, 
hence, the images with invalid IMGQS can be updated with 0099.

All the images could be cleaned from bias, darks and flats ([see issues
regarding flats](#data-issues)). However, no catalog could be
produced for 4 images:

<div align="center">

| image | object |
| ----- | ------ |
| 2457071/TJO2457071.69501_R_iml.fits.gz | 12 Oph |
| 2457071/TJO2457071.69892_R_iml.fits.gz | 12 Oph |
| 2457071/TJO2457071.70402_R_iml.fits.gz | 12 Oph |
| 2457116/TJO2457116.45622_U_iml.fits.gz | OT Ser |

</div>

Therefore, the resulting reduction provided:

<div align="center">

| object     | night   | images | stars |
| ---------- | ------- | ------ | ----- |
| 12 Oph     | 2457071 |     47 |    64 |
| 12 Oph     | 2457073 |      6 |    10 |
| 12 Oph Ref | 2457073 |      3 |    47 |
| HD 135599  | 2457032 |     14 |   125 |
| HD 135599  | 2457033 |     14 |   107 |
| HD 135599  | 2457034 |     14 |   131 |
| HD 135599  | 2457036 |     14 |   135 |
| HD 135599  | 2457037 |      4 |    24 |
| HD 135599  | 2457045 |     12 |    73 |
| HD 135599  | 2457046 |      6 |    33 |
| HD 180617  | 2457073 |      4 |   386 |
| HD 180617  | 2457081 |     15 |  1129 |
| HD 98712   | 2457061 |     10 |    84 |
| HD 98712   | 2457062 |      8 |    70 |
| HD 98712   | 2457063 |     30 |   425 |
| HD 98712   | 2457064 |     30 |   426 |
| OT Ser     | 2457109 |     75 |  1416 |
| OT Ser     | 2457110 |     75 |  1176 |
| OT Ser     | 2457111 |     75 |   949 |
| OT Ser     | 2457112 |     75 |  1308 |
| OT Ser     | 2457113 |     75 |  1483 |
| OT Ser     | 2457114 |     45 |   680 |
| OT Ser     | 2457115 |     25 |   364 |
| OT Ser     | 2457116 |     38 |   142 |
| OT Ser     | 2457117 |     50 |   541 |
| V2306 Oph  | 2457109 |     25 |  1101 |
| V2306 Oph  | 2457110 |     25 |  1228 |
| V2306 Oph  | 2457111 |     25 |  1036 |
| V2306 Oph  | 2457112 |     27 |  1760 |
| V2306 Oph  | 2457113 |     27 |  1399 |
| V2306 Oph  | 2457114 |      7 |   184 |
| V2306 Oph  | 2457115 |      7 |   188 |
| V2306 Oph  | 2457117 |      7 |   146 |

</div>

