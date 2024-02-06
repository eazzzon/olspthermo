import setuptools

setuptools.setup(
    name="OSAT",
    version="v0.1",
    author="Yishen Zhang",
    author_email="yishen.zhang@kuleuven.be",
    packages= setuptools.find_packages(
        exclude= ['example', 'benchmark']
        ),


    install_requires= [
    'pandas',
    'numpy',
    'scipy',
    'periodictable'
    ]
)
