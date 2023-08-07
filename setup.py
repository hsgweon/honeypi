from setuptools import setup, find_packages
import os

__version__ = os.environ.get("VERSION", "1.2")

setup(
    name = "honeypi",
    version = __version__,
    packages = ["honeypi"],
    scripts = ["bin/honeypi_createreadpairslist", "bin/honeypi", "bin/honeypi_dada2", "bin/honeypi_reformatAssignedTaxonomy", "bin/honeypi_filterASVtable", "bin/honeypi_mergeDuplicateASV", "bin/honeypi_joinTwoResults.py"],
    description = "honeypi (HONEY Pollen and Plants ITS Pipeline)",
    long_description = "An open source tools for processing NGS amplicon data.",
    author = "Hyun Soon Gweon",
    author_email = "h.s.gweon@reading.ac.uk",
	url = "https://github.com/hsgweon/honeypi",
	license = "GPA"
)
