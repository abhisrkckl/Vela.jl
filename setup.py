from distutils.core import setup

setup(
    name="pint2vela",
    version="0.0.3",
    description="Interface between PINT and Vela.jl",
    author="Abhimanyu Susobhanan",
    author_email="abhimanyu.susobhanan@nanograv.org",
    url="https://github.com/abhisrkckl/Vela.jl/",
    install_requires=["numpy", "astropy", "pint-pulsar", "juliacall"],
    package_dir={"": "pint2vela"},
    packages=["pint2vela"],
    scripts=["pint2vela/scripts/par_tim-to-jlso.py"],
)
