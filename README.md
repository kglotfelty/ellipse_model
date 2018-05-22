
![ds9 image](ds9.png)


# Elliptical model

This script takes the output from the [CIAO](http://cxc.cfa.harvard.edu/ciao/index.html)
[dmellipse](http://cxc.cfa.harvard.edu/ciao/ahelp/dmellipse.html) tool
and uses it to create a model image.

## Example

### Run `dmellipse`

We start by running `dmellipse` to create an ellipse profile of 
the image.  `dmellipse` takes the input image `img.fits` and a 
[stack](http://cxc.cfa.harvard.edu/ciao/ahelp/stack.html) of 
integrated fractions.  

```bash
dmellipse img.fits ellipses.fits fraction="lgrid(0.05:1.0:0.025)" step=20 mode=h clob+
```

This example generates the ellipse profile which enclose from 5% to 100%
of the total flux in the image, in steps of 2.5%.  The `step` parameter is
increased to match the scale of the image which allows the tool to run
faster.

```bash
ds9 img.fits -region ellipses.fits -scale log \
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  -view colorbar no  -saveimage png ellipses.png -quit
```

![Ellipses](ellipses.png)


### Reconstruct using `emodel`

Now we will reconstruct a model of the image based on the 
ellipse profile

```bash
emodel ellipses.fits img.fits model.fits weight=flat clob+ 
```

`emodel` takes the information about the ellipses and the fraction 
of the total flux in each ellipse to reconstruct a 2D model.  The `weight`
parameter controls how the flux is distributed within the ellipse.  The
default `flat` equally distributes the flux into each pixel in the ellipse.



```bash
ds9 model.fits -scale log \
  -cmap load ~/ds9_hacks/LUT/Neota/neota_sunset-in-atlantis.lut \
  -view colorbar no  -saveimage png model.png -quit
```

![Model](model.png)


