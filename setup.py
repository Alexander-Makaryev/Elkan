from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

if __name__ == '__main__':
	module1 = Extension('elkan', ['elkan_wrapper.pyx'], language="c++"
						)

	setup (name = 'elkan',
		   version = '1.0',
		   ext_modules = cythonize([module1]),
           include_dirs=[numpy.get_include()])