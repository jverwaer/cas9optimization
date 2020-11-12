# -*- coding: utf-8 -*-

import setuptools

with open("./README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cas9_optimization_JVERWAER", # Replace with your own username
    version="0.0.1",
    author="Jan Verwaeren",
    author_email="jan.verwaeren@ugent.be",
    description="A package for codon-optimization of cas9",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)