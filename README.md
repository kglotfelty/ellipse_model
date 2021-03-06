![ds9 image](ds9.png)

<!--
ds9 -geometry 1135x670 img2.fits -scale log -tile column \
  -scale limits 0 175 \
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  img2.fits -region ellipses.fits -scale limits 0 175 \
  model.fits -scale limits 0 175 \
  -view colorbar no  -saveimage png ds9.png -quit
-->


# Elliptical model

This script takes the output from the [CIAO](http://cxc.cfa.harvard.edu/ciao/index.html)
[dmellipse](http://cxc.cfa.harvard.edu/ciao/ahelp/dmellipse.html) tool
and uses it to create a model image.

## Example

### Run `dmellipse`

We start by running `dmellipse` to create an ellipse profile of 
the image.  `dmellipse` takes the input image `img2.fits` and a 
[stack](http://cxc.cfa.harvard.edu/ciao/ahelp/stack.html) of 
integrated fractions.  

In this example I've filtered the input image with a region so that
the edge of the image is stored in the file's data-subspace.

```bash
dmimghull img.fits hull.fits clob+
dmcopy "img.fits[sky=region(hull.fits)][opt full]" img2.fits clob+
dmellipse img2.fits ellipses.fits fraction="lgrid(0.05:1.0:0.025)" step=20 mode=h clob+
```

This example generates the ellipse profile which enclose from 5% to 100%
of the total flux in the image, in steps of 2.5%.  The `step` parameter is
increased to match the scale of the image which allows the tool to run
faster.

```bash
ds9 img2.fits -region ellipses.fits -scale log \
  -scale limits 0 175 \
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  -view colorbar no  -saveimage png ellipses.png -quit
```


![Ellipses](ellipses.png)


### Reconstruct using `emodel`

Now we will reconstruct a model of the image based on the 
ellipse profile

```bash
emodel ellipses.fits img2.fits model.fits weight=flat clob+ norm=486119
```

`emodel` takes the information about the ellipses and the fraction 
of the total flux in each ellipse to reconstruct a 2D model.  The `weight`
parameter controls how the flux is distributed within the ellipse.  The
default `flat` equally distributes the flux into each pixel in the ellipse.



```bash
ds9 model.fits -scale log -scale limits 0 175 \
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  -view colorbar no  -saveimage png model.png -quit
```

![Model](model.png)


### Weights

**Experimental**

The `weight` option is still experimental.  The `flat` and `linear` 
are well behaved but the other options may produces some artifacts, 
especially at the boundary between ellipses.

```bash
emodel ellipses.fits img2.fits model_flat.fits weight=flat clob+ norm=486119
emodel ellipses.fits img2.fits model_linear.fits weight=linear clob+ norm=486119
emodel ellipses.fits img2.fits model_gaussian.fits weight=gaussian clob+ norm=486119
emodel ellipses.fits img2.fits model_square.fits weight=square clob+ norm=486119
emodel ellipses.fits img2.fits model_sqrt.fits weight=sqrt clob+ norm=486119 
emodel ellipses.fits img2.fits model_hemi.fits weight=hemi clob+ norm=486119
```

The effect of the weighting is most obvious in the center of the
image


```
ds9 model_flat.fits -scale log -scale limits 0 175 -zoom 4 -tile -view colorbar no\
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  model_linear.fits model_gaussian.fits model_square.fits model_sqrt.fits model_hemi.fits \
  -pan to 4270 3960 physical -match frames physical \
  -saveimage png weights.png -quit
```

![weights](weights.png)


## Installation

This script requires [ciao-4.10](http://cxc.cfa.harvard.edu/ciao/download/index.html)
or later.   It will not work with earlier versions of CIAO.

Users must have sourced the CIAO setup script before installing this
script.

```bash
$ source /soft/ciao/bin/ciao.bash
CIAO configuration is complete... 
CIAO 4.10 Thursday, April 12, 2018
  bindir      : /soft/ciao-4.10/bin
  CALDB       : 4.7.8
```

replace `/soft/ciao` with your location where you installed CIAO, then

```bash
$ git clone https://github.com/kglotfelty/ellipse_model 
Cloning into 'ellipse_model'...
remote: Counting objects: 110, done.
remote: Compressing objects: 100% (76/76), done.
remote: Total 110 (delta 53), reused 86 (delta 29), pack-reused 0
Receiving objects: 100% (110/110), 960.65 KiB | 0 bytes/s, done.
Resolving deltas: 100% (53/53), done.
Checking connectivity... done.

$ cd ellipse_model/
$ python setup.py install
```

