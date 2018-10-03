from setuptools import setup, find_packages

setup(name='jacks',
      version='0.2',
      description='JACKS package for processing CRISPR/Cas9 Screens',
      url='http://github.com/felicityallen/JACKS',
      author='Felicity Allen and Leopold Parts',
      license='MIT',
      include_package_data=True,
      packages=find_packages(),
      install_requires=['scipy','numpy>=1.9.0','matplotlib'],
      zip_safe=False)
      
      