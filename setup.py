import pathlib

from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="lsq-ellipse",
    version="2.2.1",
    description="Lease squares fitting of Ellipses",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/bdhammel/least-squares-ellipse-fitting",
    author="Ben Hammel",
    author_email="bdhammel@gmail.com",
    license="MIT",
    py_modules=['ellipse'],
    install_requires=["numpy"],
)
