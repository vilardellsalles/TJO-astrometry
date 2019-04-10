# Documentation

## Introduction

This repo contains the steps followed to reduce the data from the TJO proposal
93 entitled "First astrometric measurements of stellar activity with 
spectroscopy". The data for this project can be summarized in the following 
table:

| object     | instrume | filter | Number of images |
| ---------- | -------- | ------ | -------- |
| 12 Oph     | MEIA     | B      |       56 |
| 12 Oph Ref | MEIA     | B      |        3 |
| HD 135599  | MEIA     | B      |       78 |
| HD 180617  | MEIA     | B      |       19 |
| HD 98712   | MEIA     | B      |       78 |
| OT Ser     | MEIA2    | V      |      534 |
| V2306 Oph  | MEIA2    | V      |      150 |

## Reduction

The following commands are required to clean and reduce the images:

```
icat gaia-imgqs <images>
icat gaia-update --imgqs [0099|209] <images with invalid IMQS>
icat gaia-reduce *
icat light-curve *_catalog/*.dat
```

IMGQS is a keyword inserted in the FITS header after ```gaia-imgqs``` that 
contains several quality indicators. All the images used are valid and, hence, 
the IMGQS can be updated with either 0099 (science images) or 209 (calibration 
images).

## Reduction details

### HD 98712

I have observed that the focus position for this object is not constant and 
several images where obtained with a different focus:

| focus | images   |
| ----- | -------- |
|  3300 |        4 |
|  3450 |       62 |
|  3600 |        4 |
|  3900 |        4 |
|  4050 |        2 |
|  4200 |        2 |

Therefore, some of the images of this target will be rejected.

### 12 Oph

This target has no comparison stars in the field of view and therefore, it was 
supposed to have a reference exposure (12 Oph Ref) for each target exposure.
However, only reference images were obtained for one of the two nights that
12 Oph was observed. In addition, only 6 images were obtained on the same night 
(of the 56 obtained).