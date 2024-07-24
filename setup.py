from distutils.core import setup

setup(
    name="pint2vela",
    version="0.0.1",
    description="Interface between PINT and Vela.jl",
    author="Abhimanyu Susobhanan",
    author_email="abhimanyu.susobhanan@nanograv.org",
    url="https://github.com/abhisrkckl/pint2vela/",
    install_requires=["numpy", "astropy", "pint-pulsar", "juliacall"],
    packages=["pint2vela"],
    scripts=["scripts/par_tim-to-jlso.py"]
)
