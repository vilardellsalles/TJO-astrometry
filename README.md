# Documentation

## Introduction

This repo contains the steps followed to reduce the data from the TJO proposal
93 entitled "First astrometric measurements of stellar activity with 
spectroscopy". The data for this project can be summarized in the following 
table:

<div align="center">

| object     | instrume | filter | focus range | images | nights |
| ---------- | -------- | ------ | ----------- | ------ | ------ |
| 12 Oph     | MEIA     | B      | 4464 - 4474 |     56 |      2 |
| 12 Oph Ref | MEIA     | B      | 4466 - 4467 |      3 |      1 |
| HD 135599  | MEIA     | B      | 3935 - 4054 |     78 |      7 |
| HD 180617  | MEIA     | B      | 3217 - 3463 |     19 |      2 |
| HD 98712   | MEIA     | B      | 3440 - 4231 |     78 |      4 |
| OT Ser     | MEIA2    | V      | 7859 - 8884 |    534 |      9 |
| V2306 Oph  | MEIA2    | V      | 7865 - 8871 |    150 |      8 |

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

Data is provided in compressed *tar* files. One file is provided for each
observation night. Compressed files have to be uncompressed to have FITS images with the following command:

```
cd obsdata
for night in `ls *.tar.gz`; do tar xzvf $night; done
cd ..
```

Extracted directories, can then be processed with specific
[ICAT](https://gitlab.ice.csic.es/rocs/icat) pipelines:

```
for image in `ls obsdata/2015*/*.fits*; do icat gaia-imgqs $image; done
icat gaia-reduce rawdata/*
```

With the commands above, images are stored into specific directories and in a
custom MySQL database. They are also cleaned from bias, darks and flats ([see
issues regarding flats](#data-issues)). Photometry and astrometry is also
extracted for all stars and stored in a SExtractor file (one file per reduced
image) and in the database. Unfortunately, a valid catalog could not be
produced for 12 of the 12 Oph images. Considering the [lack of reference
stars](#data-issues), no further attempts were done to reduce these images.

Therefore, the resulting reduction provided:

<div align="center">

| object     | night   | images | stars | stars per image |
| ---------- | ------- | -----  | ----- | --------------- |
| 12 Oph     | 2457071 |     38 |    38 |             1.0 |
| 12 Oph     | 2457073 |      6 |     6 |             1.0 |
| 12 Oph Ref | 2457073 |      3 |    29 |             9.7 |
| HD 135599  | 2457032 |     14 |   125 |             8.9 |
| HD 135599  | 2457033 |     14 |   107 |             7.6 |
| HD 135599  | 2457034 |     14 |   131 |             9.4 |
| HD 135599  | 2457036 |     14 |   135 |             9.6 |
| HD 135599  | 2457037 |      4 |    24 |             6.0 |
| HD 135599  | 2457045 |     12 |    73 |             6.1 |
| HD 135599  | 2457046 |      6 |    33 |             5.5 |
| HD 180617  | 2457073 |      4 |   705 |           176.3 |
| HD 180617  | 2457081 |     15 |  1498 |            99.9 |
| HD 98712   | 2457061 |     10 |   115 |            11.5 |
| HD 98712   | 2457062 |      8 |    83 |            10.4 |
| HD 98712   | 2457063 |     30 |   855 |            28.5 |
| HD 98712   | 2457064 |     30 |   923 |            30.8 |
| OT Ser     | 2457109 |     75 |  2080 |            27.7 |
| OT Ser     | 2457110 |     75 |  1943 |            25.9 |
| OT Ser     | 2457111 |     75 |  1584 |            21.1 |
| OT Ser     | 2457112 |     75 |  1754 |            23.4 |
| OT Ser     | 2457113 |     75 |  1991 |            26.5 |
| OT Ser     | 2457114 |     45 |  1060 |            23.6 |
| OT Ser     | 2457115 |     25 |   625 |            25.0 |
| OT Ser     | 2457116 |     39 |   272 |             7.0 |
| OT Ser     | 2457117 |     50 |   723 |            14.5 |
| V2306 Oph  | 2457109 |     25 |  1474 |            59.0 |
| V2306 Oph  | 2457110 |     25 |  1759 |            70.4 |
| V2306 Oph  | 2457111 |     25 |  1395 |            55.8 |
| V2306 Oph  | 2457112 |     27 |  2690 |            99.6 |
| V2306 Oph  | 2457113 |     27 |  1684 |            62.4 |
| V2306 Oph  | 2457114 |      7 |   227 |            32.4 |
| V2306 Oph  | 2457115 |      7 |   249 |            35.6 |
| V2306 Oph  | 2457117 |      7 |   211 |            30.1 |

</div>

The next step is to produce a light curve from the extracted photometry.
Currently, the pipeline responsible for that can only read the data from the
database and generate a table for each object in HTML form:

```
cd reddata
icat light-curve -o "12 Oph"
icat light-curve -o "12 Oph Ref"
icat light-curve -o "HD 135599"
icat light-curve -o "HD 180617"
icat light-curve -o "HD 98712"
icat light-curve -o "OT Ser"
icat light-curve -o "V2306 Oph"
cd ..
```
