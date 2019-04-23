from setuptools import setup

setup(name='filtools',
      version='2019.04.1',
      description='Tools for working on filterbank data',
      url='https://bitbucket.org/mkeith/filtools/',
      author='Michael Keith',
      author_email='mkeith@pulsarastrononmy.net',
      license='GPLv3',
      packages=['filtools','filtools.inject_pulsar'],
      package_data={'filtools':['site_data/*.yaml']},
      scripts=['bin/ft_inject_pulsar'],
      zip_safe=False)
