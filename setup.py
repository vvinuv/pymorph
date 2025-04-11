import os
import sys
from setuptools import setup, find_packages
from setuptools.command.easy_install import easy_install

# Check if pymouse could run on the given system
if os.name != 'posix':
    raise ValueError(
        'Detected unsupported operating system: {}.'.format(sys.platform)
    )

if sys.version_info < (3, 5):
    raise ValueError(
        'Unsupported Python version {}.{}.{} found. pymouse requires Python 3.5 or higher.'.format(sys.version_info.major, sys.version_info.minor, sys.version_info.micro)
    )


current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'requirements.txt')) as f:
    install_reqs = [r.rstrip() for r in f.readlines() if not r.startswith('#')]

with open(os.path.join(current_dir, 'requirements.txt')) as f:
    extras_reqs = [r[1:].rstrip() for r in f.readlines() if r.startswith('#')]

extras_reqs = dict(zip(extras_reqs, extras_reqs))

import _version
#with open("pymorph/_version.py") as f:
#    version = f.readlines()[-1].split()[-1].strip("\"'")
version = _version.__version__

long_description = 'This pipeline was written as a part of some projects on galaxy morphology by Vinu Vikraman and Alan Meert. Yogesh Wadadekar, Ajit Kembhavi, G V Vijayagovidan and Mariangela Bernardi were working on those projects. It REQUIRES SExtractor and GALFIT.'


setup(
    name='pymorph',
    author='Vinu Vikraman, Alan Meert',
    author_email='vvinuv@gmail.com',
    description='Galaxy Models using SExtractor and GALFIT',
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=install_reqs,
    version=version,
    packages=find_packages(),
    include_package_data=True,
    package_data={'pymorph': ['SEx/*', 'html/*']},
    #extras_require=extras_reqs,
    scripts=['scripts/pymorph_bin'],
    license='BSD',
    platforms=['Linux'],
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Galaxy Analysis",
        'Programming Language :: Python :: 3.6',
    ],
    python_requires='==3.12.2',
)

