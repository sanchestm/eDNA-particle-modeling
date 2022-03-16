import setuptools


setuptools.setup(
    name='eDNA-particle-modeling',
    version='0.0.1',
    author='Thiago Sanches',
    author_email='none',
    description='Testing installation of Package',
    url='https://github.com/sanchestm/eDNA-particle-modeling/',
    license='MIT',
    packages=['eDNA-particle-modeling'],
    install_requires=['numpy','pymc3>=3.8','arviz','scipy','ipywidgets','scikit-learn','seaborn','matplotlib','pandas','statsmodels'],
)
