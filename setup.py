__author__ = 'Stephen G. Gaffney'

from setuptools import setup

setup(name='maf_matrix',
      version='0.1',
      description='Matrix plotting for TCGA Maf files.',
      url='http://github.com/sggaffney',
      author='Stephen G. Gaffney',
      author_email='stephen.gaffney@yale.edu',
      license='GPLv3',
      packages=['maf_matrix'],
      install_requires=[
          'pandas', 'requests', 'numpy', 'matplotlib', 'future'
      ],
      test_suite='nose.collector',
      tests_require=['nose'],
      scripts=['bin/plot_maf_matrix'],
      zip_safe=False)
