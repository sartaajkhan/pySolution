from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

#with codecs.open(os.path.join(here, "README.md"), encoding="utf-8") as fh:
#    long_description = "\n" + fh.read()

VERSION = '0.0.2.1.4.3'
DESCRIPTION = 'Computes details of solution for modeling purposes'
LONG_DESCRIPTION = """
A package that was created for the purpose of process modeling in different membranes. This is a solution chemistry toolkit
that can compute details of "saline" solutions such as the activity coefficient of each ion (using the Pitzer model),
osmotic coefficient, ionic strength, density of the solution, etc.
"""

# Setting up
setup(
    name="pySolution",
    version=VERSION,
    author="Sartaaj Khan",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    url = "https://github.com/sartaajkhan/pySolution",
    packages=find_packages(),
    #package_data={'pySolution': ['core/*.xlsx', 'core/*.pickle', 'core/*.txt']},
    package_data = {'pySolution' : ['core/*']},
    install_requires=['numpy', 'pandas', 'scipy', 'joblib', 'proplot', 'matplotlib', 'scikit-learn==0.24.2'],
    keywords=['python', 'solution chemistry', 'activity coefficient', 'osmosis', 'reverse osmosis', 'chemistry'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)