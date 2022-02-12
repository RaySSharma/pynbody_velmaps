from __future__ import absolute_import
install_requires = [
    'setuptools',
    'numpy >= 1.19.2',
    'scipy >= 1.5.2',
    'matplotlib >= 3.0.0',
    'astropy >= 4.0.0',
    'pynbody >= 1.0.2'
    ]

tests_require = []

from setuptools import setup, find_packages

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(name='pynbody_velmaps',
      version='1.0.0',
      description='Generate stellar/gas velocity maps and position angles using pynbody.',
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          "Programming Language :: Python",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "License :: OSI Approved :: BSD License",
      ],
      author="Ray Sharma",
      author_email="raysharma@physics.rutgers.edu",
      license="BSD",
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6',
      install_requires=install_requires,
      tests_require=tests_require,
      test_suite="nose.collector",
      long_description=long_description,
      long_description_content_type='text/markdown'
      )
