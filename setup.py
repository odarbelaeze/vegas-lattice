import os
import setuptools


with open(os.path.join(os.path.dirname(__file__), 'README.md')) as readme:
    README = readme.read()


setuptools.setup(
    name='vegas-lattice',
    version='0.2.1',
    author='Oscar D. Arbeláez <@odarbelaeze>, Juan D. Alzate <@jdalzatec>',
    author_email='odarbelaeze@gmail.com',
    license='MIT license',
    packages=setuptools.find_packages(),
    install_requires=[
        'click>=6.2',
        'numpy>=1.9.2',
    ],
    entry_points='''
    [console_scripts]
    vegas-lattice=vlattice.scripts:cli
    ''',
    url='https://github.com/odarbelaeze/vegas-lattice',
    download_url='https://github.com/odarbelaeze/vegas-lattice/tarball/0.2.1',
    keywords=['lattice', 'graph', 'magnetite', 'construction', 'grid' ],
    description='A grid lattice sample generator',
    long_description=README,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)
