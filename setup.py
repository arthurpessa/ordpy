#old version
# from setuptools import setup

# setup(name='ordpy',
#       version='1.0.0',
#       description='A Python package for data analysis with permutation entropy and ordinal networks methods.',
#       url='https://github.com/arthurpessa/ordpy',
#       author='Arthur A. B. Pessa and Haroldo V. Ribeiro',
#       author_email='arthur_pessa@hotmail.com, hvr@dfi.uem.br',
#       license='MIT',
#       packages=['ordpy'],
#       install_requires=['numpy'],
#       python_requires=">=3.6",
#       zip_safe=False
# )

import setuptools

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ordpy",
    version="1.0.6",
    author="Arthur A. B. Pessa and Haroldo V. Ribeiro",
    author_email="arthur_pessa@hotmail.com, hvr@dfi.uem.br",
    description="A Python package for data analysis with permutation entropy and ordinal networks methods.",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    url="https://github.com/arthurpessa/ordpy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
