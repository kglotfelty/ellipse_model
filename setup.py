

import os
import sys

assert "ASCDS_INSTALL" in os.environ, "Please setup for CIAO before installing"

# Set PYVER env variable so that it'll get replaced in setup.cfg
ver = sys.version_info
os.environ["PYVER"] = "python{}.{}".format(ver[0],ver[1]) 


from distutils.core import setup

setup( name='emodel',
       version='0.1.0',
       description='Elliptical model builder',
       author='Anonymous',
       author_email='WhoDat@cfa.harvard.edu',
       url='https://github.com/kglotfelty/ellipse_model/',
       py_modules=["masked_image_crate"],
       scripts=["emodel"],
       data_files=[('param',['emodel.par']),
                    ('share/doc/xml',['emodel.xml'])
                    ]
                    
    )
