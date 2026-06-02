from setuptools import setup

__author__ = "Yasas Wijesekara"
__copyright__ = "Copyright 2024, Yasas Wijesekara"
__email__ = "yasas.Wijesekara@uni-greifswald.de"
__license__ = "BSD-3"

# read the contents of your README file
from os import path, listdir

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

scripts = [path.join("scripts", i) for i in listdir(path.join(this_directory, "scripts")) if i.endswith(".py")]

setup(
    name="skadi",
    version="0.0.4",
    url="https://github.com/Yasas1994/skadi",
    license=__license__,
    author=__author__,
    author_email=__email__,
    zip_safe=False,
    description="SKADI: Sequence-based Knowledgebase for Annotation, Detection, and Identification",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=["skadi"],
    include_package_data=True,
    package_data={
        "skadi": ["pipeline/*", "pipeline/rules/*"],
    },
    python_requires=">=3.10",
    install_requires=[
        "click>=8.0",
        "polars>=0.20",
        "pandas>=1.5",
        "pyyaml>=6.0",
        "tqdm>=4.60",
        "taxopy>=0.15",
        "numpy>=1.24",
        "biopython>=1.80",
        "requests>=2.28",
        "rich>=13.0",
        "pyfastx>=2.0",
    ],
    entry_points={"console_scripts": ["skadi = skadi.cli:cli"]},
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    scripts=scripts,
)
