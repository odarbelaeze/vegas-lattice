import os
import setuptools


with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as readme:
    README = readme.read()


setuptools.setup(
    name='vegas-lattice',
    version='0.4.0',
    description='A grid lattice sample generator',
    long_description=README,
    author='Oscar D. Arbeláez <@odarbelaeze>, Juan D. Alzate <@jdalzatec>',
    author_email='odarbelaeze@gmail.com',
    packages=setuptools.find_packages(),
    entry_points='''
    [console_scripts]
    vegas-lattice=vlattice.scripts:cli
    ''',
    install_requires=[
        'click>=6.2',
        'numpy>=1.9.2',
    ],
    include_package_data=True,
    license='MIT license',
    url='https://github.com/odarbelaeze/vegas-lattice',
    download_url='https://github.com/odarbelaeze/vegas-lattice/tarball/0.4.0',
    keywords=['lattice', 'graph', 'magnetite', 'construction', 'grid', ],
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
