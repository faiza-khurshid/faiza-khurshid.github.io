#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = [
    'Click',
    'networkx',
    'matplotlib',
    'pandas',
    'requests',
    'tqdm',
    'scipy',
    'sqlalchemy',
    'xmltodict',
]

test_requirements = ['pytest>=3', ]

setup(
    author="Bruce Schultz",
    author_email='bruce.schultz@scai.fraunhofer.de',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Exercise06 Package",
    entry_points={
        'console_scripts': [
            'plab2=plab2.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='plab2',
    name='plab2',
    packages=find_packages(include=['plab2', 'plab2.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
