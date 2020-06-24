

import os
import sys

assert "ASCDS_INSTALL" in os.environ, "Please setup for CIAO before installing"

# Set PYVER env variable so that it'll get replaced in setup.cfg
ver = sys.version_info
os.environ["PYVER"] = "python{}.{}".format(ver[0],ver[1]) 


from setuptools import setup
from setuptools.command.install import install


class InstallAhelpWrapper(install):
    'A simple wrapper to run ahelp -r after install to update ahelp index'

    @staticmethod
    def update_ahelp_database():
        print("Update ahelp database ...")
        from subprocess import check_output
        sout = check_output(["ahelp","-r"])
        for line in sout.decode().split("\n"):
            for summary in ["Processed", "Succeeded", "Failed", "Purged"]:
                if line.startswith(summary):
                    print("    "+line)

    
    def run(self):
        install.run(self)
        self.update_ahelp_database()


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
                    ],
        cmdclass={'install': InstallAhelpWrapper},
                    
    )

