#!/usr/bin/env python

from distutils.core import setup

setup(name='basilisk',
      version='0.1',
      description='A generative, probabilistic model of side chains in proteins.',
      author='Tim Harder and Jes Frellsen',
      author_email='mail@tim-harder.de',
      url='http://www.binf.ku.dk/',
      license='GPLv3',
      packages=['basilisk_lib',
                'basilisk_lib.Mocapy',
                ],
     )
