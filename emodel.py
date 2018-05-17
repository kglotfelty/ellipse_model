

from region import *
import numpy as np
from pycrates import read_file




infile = "ellipses.fits[#row=1:3]"
infile="k.fits"
image_file = "img.fits"
outfile = "model.fits"
normalization = 1.0
clobber = True


def load_ellipses( infile ):

    def scalar(x):
        if len(x.shape) > 1:
            x = x[:,0]
        return x
    
    tab = read_file(infile,"r")
    
    s = tab.get_column("shape").values
    if not all([_s.lower() == 'ellipse' for _s in s]):
        raise ValueError("Only ellipses, inclusive ellipses, are supported")
    
    x = scalar(tab.get_column("x").values)
    y = scalar(tab.get_column("y").values)
    r = tab.get_column("r").values
    r0 = r[:,0]
    r1 = r[:,1]
    a = scalar(tab.get_column("rotang").values)
    f = tab.get_column("fraction").values
    
    ee = [ ellipse(_x,_y,_r0,_r1,_a) for _x,_y,_r0,_r1,_a in zip(x,y,r0,r1,a)]
    retval = [ (_f,_e) for _f,_e in zip(f,ee)]
    retval.sort()
    return retval


def load_image(infile):
    
    img = read_file(infile,"r")
    pixels = img.get_image().values
    img.get_image().values = np.zeros_like(pixels)*1.0

    _a = [x.lower() for x in img.get_axisnames()]

    if 'sky' in _a:
        sky = 'sky'
    elif 'pos' in _a:
        sky = 'pos'
    else:
        raise IOError("Unknown sky column in image")
    
    return(img,sky)


def region_limits( e, xmax, ymax ):
    b = e.extent()
    retval = img.get_transform(sky).invert([[b["x0"],b["y0"]],[b["x1"],b["y1"]]])
    retval = retval-1 # Convert to 0 based array indexing
    retval = np.floor(retval)
    retval = retval.astype("i4") + np.array( [[0,0],[1,1]]) # +1 to upper limit
    retval = np.clip(retval, 1, (xmax,ymax))
    return retval


def convert_coordinates( img ):
    ylen,xlen = img.get_image().values.shape
    
    ivals = [x+1.0 for x in range(xlen)]
    jvals = [y+1.0 for y in range(ylen)]
    ivals = np.array( ivals*ylen)
    jvals = np.repeat( jvals,xlen)
    ijvals = list(zip(ivals,jvals))
    xyvals = img.get_transform(sky).apply( ijvals )
    ijvals = [(int(i),int(j)) for i,j in ijvals ]
    coords = dict(zip(ijvals, xyvals))
    return coords


def union_previous_shapes( edata, N ):

    shapes_already_counted = CXCRegion()
    for e in edata[0:N]:
        shapes_already_counted = shapes_already_counted+e[1]
    
    return(shapes_already_counted)
    


    
edata = load_ellipses( infile )
img,sky = load_image(image_file)

ylen,xlen = img.get_image().values.shape
coords = convert_coordinates(img)

ovals = img.get_image().values


for N in range(0,len(edata)):

    fN = edata[N][0] # fraction in 1st region
    eN = edata[N][1] # the region itself


    shapes_already_counted = union_previous_shapes(edata,N)

    # Intersect current shape with existing ones.
    region_already_counted = eN * shapes_already_counted

    # Now we compute the elliptical "annulus" to fill in.
    # We have to treat the 1st shape special.
    if 0 == len(shapes_already_counted.shapes):
        region_not_already_counted = eN
    else:
        region_not_already_counted = eN - shapes_already_counted


    # Determine sum of fraction already included 
    # in the intersection of the current ellipse and 
    # the previous ellipses.  
    erange = region_limits(region_already_counted,xlen, ylen)
    current_frac = 0.0
    for _y in range( erange[0][1], erange[1][1]+1):
        for _x in range(erange[0][0], erange[1][0]+1):
            pos = coords[(_x,_y)]
            if region_already_counted.is_inside( pos[0],pos[1]):
                current_frac = current_frac + ovals[_y-1,_x-1]

    # The delta is the fraction 
    delta_fN = fN - current_frac
    if delta_fN < 0:
        raise ValueError("Somethings wrong, no negative flux please")


    
    # Now we count how many pixels in the current image are 
    # in the annulus.  We could use the .area() method but 
    # this gets us what we need.
    #
    erange = region_limits(region_not_already_counted,xlen,ylen)
    inside = []
    weight = []
    for _y in range( erange[0][1], erange[1][1]+1):
        for _x in range(erange[0][0], erange[1][0]+1):
            pos = coords[(_x,_y)]
            if region_not_already_counted.is_inside( pos[0],pos[1]):
                inside.append (  (_x-1,_y-1), )
                weight.append( 1.0 )

    #
    # Due to pixel size, there may not be any pixels so we check
    # otherwise, we evenly distibuted the delta-flux into that many
    # pixels.
    #
    if len(inside)>0:
        val = delta_fN / sum(weight)

        # Save the values.
        for _ij,_w in zip(inside,weight):
            _i,_j = _ij
            ovals[_j,_i] = ovals[_j,_i]+(val*_w)
        
    # DEBUG: img.write("{:03d}.fits".format(N),clobber=True)


img.get_image().values = ovals * normalization

img.write(outfile, clobber=clobber)


    



