from setuptools import setup

setup(
    name="pyvela",
    version="0.0.3",
    description="Interface between PINT and Vela.jl",
    author="Abhimanyu Susobhanan",
    author_email="abhimanyu.susobhanan@nanograv.org",
    url="https://github.com/abhisrkckl/Vela.jl/",
    install_requires=["numpy", "astropy", "pint-pulsar", "juliacall"],
    package_dir={"": "pyvela"},
    packages=["pyvela"],
    scripts=["pyvela/scripts/par_tim-to-jlso.py"],
)
