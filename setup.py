from setuptools import setup
from setuptools import find_packages


setup(
    name='vegas-lattice',
    version='0.1.0',
    author='Oscar D. Arbeláez <@odarbelaeze>, Juan D. Alzate <@jdalzatec>',
    author_email='odarbelaeze@gmail.com',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.9.2',
        'click>=6.2',
    ],
    entry_points='''
    [console_scripts]
    vegas-lattice=vlattice.scripts:cli
    ''',
    url='https://github.com/odarbelaeze/vegas-lattice',
    download_url='https://github.com/odarbelaeze/vegas-lattice/tarball/0.1.0',
    keywords=['lattice', 'graph', 'magnetite', 'construction', ],
    description='A magnetite sample generator',
)
