"""
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from setuptools import setup

setup(
    name = "fitDNM",
    version = "0.2.1",
    author = "Yu Jiang",
    author_email = "yu.jiang@duke.edu",
    description = ("Enrichment of de novo mutations within genes."),
    license = "MIT",
    packages=["fitDNM"],
    install_requires=['pandas >= 0.15.0',
                      'scipy >= 0.15.0',
                      'numpy >= 1.9.0',
                      'denovonear >= 0.9.9',
                      'pysam >= 0.9.0'
    ],
    entry_points={'console_scripts': ['fitdnm = fitDNM.__main__:main']},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
    ],
    test_suite="tests",
    tests_require=[
        'coverage',
    ]
)
