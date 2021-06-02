from setuptools import setup
setup(
    name="clear_sky_models",
    version="0.1",
    description="Python converter of clear sky modules coded in R acquired from Jamie Brights GitHub: https://github.com/JamieMBright/clear-sky-models.",
    url="https://github.com/jonas-witthuhn/clear-sky-models",
    license="CC BY-NC",
    author="Jonas Witthuhn",
    author_email="witthuhn@tropos.de",
    packages=["clear_sky_models"],
    package_dir={"":"src"},
    package_data={"clear_sky_models":["R/*.R","R/*.txt"]},
    install_requires=["numpy",
                      "rpy2",
                      ]
        )
