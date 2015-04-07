from setuptools import setup
from distutils.core import Extension
from glob import glob
import os
import shutil

libraspa2_so = Extension("RASPA2.simulations.lib.libraspa2",
                         sources=glob("src/*.c"),
                         include_dirs=["src"])

raspa_data = (glob("forcefield/*/*.def") +
              glob("framework/*/*.def") +
              glob("molecules/*/*.def"))

# The structures have a weird install pattern that has to be manually managed
# This copies files around to match the proper directory format.
structure_types = ["cif", "block", "ions"]
for structure_type in structure_types:
    structures = glob("structures/*/{t}/*.{t}".format(t=structure_type))
    new_path = "structures/{t}".format(t=structure_type)
    if not os.path.exists(new_path):
        os.mkdir(new_path)
    for structure in structures:
        shutil.copy(structure, new_path)
    raspa_data += glob("{p}/*.{t}".format(p=new_path, t=structure_type))

setup(
    name="RASPA2",
    version="2.0.2",
    description="A general purpose classical simulation package.",
    url="http://github.com/numat/RASPA2/",
    author="David Dubbeldam",
    author_email="D.Dubbeldam@uva.nl",
    package_dir={"RASPA2": "python",
                 "RASPA2.share.raspa": "."},
    package_data={"RASPA2.share.raspa": raspa_data},
    include_package_data=True,
    packages=["RASPA2", "RASPA2.share.raspa"],
    entry_points={
        "console_scripts": [
            "simulate = RASPA2:run_command_line",
            "raspa-dir = RASPA2:get_raspa_dir"
        ]
    },
    ext_modules=[libraspa2_so],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Programming Language :: C",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Framework :: IPython",
        "Topic :: Scientific/Engineering :: Chemistry"
    ]
)

# This deletes the intermediate structure data files.
for structure_type in structure_types:
    shutil.rmtree("structures/{t}".format(t=structure_type))
