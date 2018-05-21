#! /usr/bin/env python

#
# Copyright (C) 2018 Smithsonian Astrophysical Observatory
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import print_function

toolname = "emodel"
__revision__ = "21 May 2018"

import os

import ciao_contrib.logger_wrapper as lw
lw.initialize_logger(toolname)
lgr = lw.get_logger(toolname)
verb0 = lgr.verbose0
verb1 = lgr.verbose1
verb2 = lgr.verbose2
verb3 = lgr.verbose3
verb4 = lgr.verbose4
verb5 = lgr.verbose5


import sys


from region import *
import numpy as np
from pycrates import read_file





class ModelImage( object ):
    """
    
    """
    def __init__( self, infile ):
        self.infile = infile
        self.ovals = None
        self.xlen = None
        self.ylen = None
        self.sky = None
        
        self.load_image()
        self.get_sky_axis()
        self.convert_coordinates()

        
    def load_image(self):
        """Load image
        
        We don't care about the image pixel values.
        
        TBD OK -- we do.  I will need them at some point to check for NaN/NULLs.
        
        The same input crate is used for the output.  
        
        TBD: Units/Content/Creator/etc.
        
        """

        self._img_crate = read_file(self.infile,"r")
        ### TBD: remove pixels, replace with self.ovals?
        pixels = self._img_crate.get_image().values
        self._img_crate.get_image().values = np.zeros_like(pixels)*1.0
        self.ovals = self._img_crate.get_image().values
        self.ylen,self.xlen = self.ovals.shape


    def get_sky_axis(self):
        """Determine sky coord axis name"""
        
        _a = [x.lower() for x in self._img_crate.get_axisnames()]
        if 'sky' in _a:
            self.sky = 'sky'
        elif 'pos' in _a:
            self.sky = 'pos'
        else:
            raise IOError("Unknown sky column in image")
        

    def convert_coordinates(self):
        """
        Convert pixel coordinates to physical coordinates

        Create (i,j) -> (x,y) mapping 
        """
        ivals = [x+1.0 for x in range(self.xlen)]
        jvals = [y+1.0 for y in range(self.ylen)]
        ivals = np.array( ivals*self.ylen) # this is why using lists instead of np.arrays
        jvals = np.repeat( jvals,self.xlen)
        ijvals = list(zip(ivals,jvals))
        
        # Due to memory leaks/etc, best to send all i,j in at once.
        xyvals = self._img_crate.get_transform(self.sky).apply( ijvals )

        # Cast to ints.  We don't want to before now, or crates/transform
        # will use integers for x,y which we don't want.        
        ijvals = [(int(i),int(j)) for i,j in ijvals ]
        
        self.coords = dict(zip(ijvals, xyvals))

    def region_logical_limits( self, extent):
        """
        Convert extent from physical to image array indexes, clipping at
        image boundaries.
        """

        limits = np.array( [[extent["x0"],extent["y0"]],
                            [extent["x1"],extent["y1"]] ])                
        retval = self._img_crate.get_transform(self.sky).invert(limits)
        retval = retval-1 # Convert to 0 based array indexing
        retval = np.floor(retval)
        retval = retval.astype("i4") + np.array( [[0,0],[1,1]]) # +1 to upper limit
        retval = np.clip(retval, 1, (self.xlen,self.ylen))

        logical_limits = { 'xlo' : retval[0][0],
                           'ylo' : retval[0][1],
                           'xhi' : retval[1][0],
                           'yhi' : retval[1][1] }

        return logical_limits


    def write( self, ovals, outfile, normalization, clobber):
        """
        
        """
        self._img_crate.get_image().values = ovals * normalization
        self._img_crate.write(outfile, clobber=clobber)



class Weights(object):
    """
    
    """
    def __init__(self, bias=None):
        self.bias = bias
    
    def set_region(self, reg):
        self.reg = reg
        self.decompose_region()

    def decompose_region( self ):        
        reg = self.reg.shapes[0]
        self.x0 = reg.xpoints[0]
        self.y0 = reg.ypoints[0]
        self.r0 = reg.radii[0]
        self.r1 = reg.radii[1]
        self.ang = reg.angles[0]
        self.ang *= 1.0        
        self.cos_ang = np.cos(np.deg2rad(self.ang))
        self.sin_ang = np.sin(np.deg2rad(self.ang))
    
    def calc_weight( self, x,y):
        raise NotImplementedError("Implement in derived class")

    def calc_distance( self, x,y):
        
        xh = x-self.x0
        yh = y-self.y0
        xr =      xh*self.cos_ang + yh*self.sin_ang
        yr = -1.0*xh*self.sin_ang + yh*self.cos_ang
        xs = xr/self.r0
        ys = yr/self.r1

        dist = np.hypot( ys, xs)
        assert dist <= 1.0, "Whoops, wrong scale"

        dist_from_edge = 1.0 - dist        
        assert dist_from_edge >= 0, "Whoops"
        
        return dist_from_edge

        
class FlatWeight( Weights ):
    
    def __init__(self, bias=None):
        super( FlatWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        return 1.0


class LinearWeight(Weights):
    def __init__(self,bias=2):
        super( LinearWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        d = self.calc_distance( x, y )
        
        # Humm, we need a bias for the weights
        d=d+self.bias
        return d

    
class SquareWeight(Weights):
    def __init__(self,bias=1):
        super( SquareWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        d = self.calc_distance( x, y )
        
        d = (d*d)+self.bias
        return d

class ISquareWeight(Weights):
    def __init__(self,bias=1):
        super( ISquareWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        d = self.calc_distance( x, y )
        d=1-d
        d=1-d*d
        d=d+self.bias
        
        return d




class SquareRootWeight(Weights):
    def __init__(self,bias=1):
        super( SquareRootWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        d = self.calc_distance( x, y )
        d= np.sqrt(d)
        d = d+self.bias
        return d


    
class ExpWeight(Weights):
    def __init__(self,bias=0.5):
        super( ExpWeight, self).__init__(bias)
        
    def calc_weight( self, x, y ):
        d = self.calc_distance( x, y )
        d=1-d
        d=d*d
        d=d*0.5
        #d=d*0.5
        #d=d*5.0
        retval = np.exp( -(d) )
        retval = retval+self.bias
        return retval




class EllipseModel(object):
    
    def __init__(self, infile, weight_class):
        self.infile = infile
        self.ellipse_list = None
        self.img = None
        
        self.load_ellipses()

        if not isinstance( weight_class, Weights):
            raise ValueError("Expecting a weight class")
        self.weight_class = weight_class


    def load_ellipses( self ):
        """Load ellipses from region file.
        
        Convert to CXCRegion ellipses so we get all the wonderful
        math logic
        """

        def _scalar(x):
            if len(x.shape) > 1:
                x = x[:,0]
            return x
        
        tab = read_file(self.infile,"r")
        s = tab.get_column("shape").values
        x = _scalar(tab.get_column("x").values)
        y = _scalar(tab.get_column("y").values)
        r = tab.get_column("r").values
        r0 = r[:,0]
        r1 = r[:,1]
        a = _scalar(tab.get_column("rotang").values)
        f = tab.get_column("fraction").values
        
        if not all([_s.lower() == 'ellipse' for _s in s]):
            raise ValueError("Only ellipses, inclusive ellipses, are supported")
        ee = [ ellipse(_x,_y,_r0,_r1,_a) for _x,_y,_r0,_r1,_a in zip(x,y,r0,r1,a)]
        retval = [ (_f,_e) for _f,_e in zip(f,ee)]
        retval.sort() 
        self.ellipse_list = retval


    def load_image( self, infile ):
        """
        Load image and setup coordinates
        """        
        self.img = ModelImage(infile)
        
        
    def union_previous_shapes( self, up_to_N ):
        """
        Combine ellipses with lower fractions.
        """
        shapes_already_counted = CXCRegion()
        for e in self.ellipse_list[0:up_to_N]:
            shapes_already_counted = shapes_already_counted+e[1]
        
        return(shapes_already_counted)
        

    def _calc_covered_frac(self, reg ):
        # Determine sum of fraction already included 
        # in the intersection of the current ellipse and 
        # the previous ellipses.          
        extent = reg.extent()
        erange = self.img.region_logical_limits( extent )
        current_frac = 0.0
        for _y in range( erange["ylo"], erange["yhi"]+1):
            for _x in range(erange["xlo"], erange["xhi"]+1):
                pos = self.img.coords[(_x,_y)]
                if reg.is_inside( pos[0],pos[1]):
                    current_frac = current_frac + self.img.ovals[_y-1,_x-1]
        return current_frac


    def calc_area_and_weights( self, reg, curshape ):
        # Now we count how many pixels in the current image are 
        # in the annulus.  We could use the .area() method but 
        # this gets us what we need.
        #
        self.weight_class.set_region( curshape )

        extent = reg.extent()
        erange = self.img.region_logical_limits( extent )
        inside = []
        weight = []
        for _y in range( erange["ylo"], erange["yhi"]+1):
            for _x in range(erange["xlo"], erange["xhi"]+1):
                pos = self.img.coords[(_x,_y)]
                if reg.is_inside( pos[0],pos[1]):
                    inside.append (  (_x-1,_y-1), )
                    weight.append( self.weight_class.calc_weight(pos[0],pos[1] ))

        return inside,weight


    def _add_to_output(self, delta_fN, inside, weight ):
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
                self.img.ovals[_j,_i] = self.img.ovals[_j,_i]+(val*_w)


    def process_region( self, N):
        """
        
        """
        fN = self.ellipse_list[N][0] # fraction in 1st region
        eN = self.ellipse_list[N][1] # the region itself

        shapes_already_counted = self.union_previous_shapes(N)

        # Intersect current shape with existing ones.
        region_already_counted = eN * shapes_already_counted
        current_frac = self._calc_covered_frac( region_already_counted)

        # The delta is the fraction 
        delta_fN = fN - current_frac
        if delta_fN < 0:
            raise ValueError("Somethings wrong, no negative flux please")

        # Now we compute the elliptical "annulus" to fill in.
        # We have to treat the 1st shape special.
        if 0 == len(shapes_already_counted.shapes):
            region_not_already_counted = eN
        else:
            region_not_already_counted = eN - shapes_already_counted

        inside, weight = self.calc_area_and_weights(region_not_already_counted, eN)

        self._add_to_output( delta_fN, inside, weight )


    def doit(self):
        """
        """        
        for N in range(len(self.ellipse_list)):
            self.process_region(N)
        

    def write( self, pars ):
        """
        TODO clobber check
        """
        self.img.write( self.img.ovals, pars["outfile"], 
            pars["normalization"], pars["clobber"])
        

def map_weights( weight_name, bias ):    
    weight_class = { 'flat' :     (FlatWeight,1.0),
                     'linear' :   (LinearWeight,2.0),
                     'square' :   (SquareWeight,1.0),
                     'gaussian' : (ExpWeight,0.5),
                     'sqrt' :     (SquareRootWeight,1.0),
                     'hemi' :     (ISquareWeight,1.0)
                     }
    wn = weight_name.lower()
    if wn not in weight_class:
        raise ValueError("Unknown weight name")

    my_weight = weight_class[wn][0]

    if bias in [None, "INDEF"]:
        bias = weight_class[wn][1]
    else:
        bias = float(bias)

    return my_weight(bias)


def process_parameters( pars ):
    """
    """
    
    



@lw.handle_ciao_errors( toolname, __revision__)
def main():

    # Load parameters
    from ciao_contrib.param_soaker import get_params
    pars = get_params(toolname, "rw", sys.argv, 
        verbose={"set":lw.set_verbosity, "cmd":verb1} )
    
    my_weight = map_weights( pars["weight"], pars["bias"])
    my_ellipse = EllipseModel( pars["infile"], my_weight)
    my_ellipse.load_image( pars["image_file"] )
    my_ellipse.doit()
    my_ellipse.write( pars ) # outfile, normalization, clobber )



if __name__ == "__main__":
    try:
        main()
    except Exception as E:
        print ("\n# "+toolname+" ("+__revision__+"): ERROR "+str(E)+"\n", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)
  

    



