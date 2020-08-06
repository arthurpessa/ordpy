from setuptools import setup

setup(name='ordpy',
      version='1.0',
      description='A Python package for data analysis with permutation entropy and ordinal networks methods.',
      url='https://github.com/hvribeiro/pyhvr',
      author='Haroldo V. Ribeiro, Arthur A. B. Pessa',
      author_email='hvr@dfi.uem.br, arthur_pessa@hotmail.com',
      license='GPL',
      packages=['ordpy'],
      install_requires=[
          'numpy',
      ],
      zip_safe=False)