import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='eDNA-particle-modeling',
    version='0.0.1',
    author='Thiago Sanches',
    author_email='none',
    description='Testing installation of Package',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/sanchestm/eDNA-particle-modeling/',
    license='MIT',
    packages=['eDNA-particle-modeling'],
    install_requires=['requests','numpy','pymc3>=3.8','arviz','scipy','ipywidgets','scikit-learn','seaborn','matplotlib','pandas','statsmodels'],
)
