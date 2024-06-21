# coding utf8
import setuptools
from seqreffilter.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="seqreffilter",
    version=get_versions(),
    author="Yuxing Xu",
    author_email="xuyuxing@mail.kib.ac.cn",
    description="Here's an example of a repository",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/SouthernCD/seqreffilter",
    include_package_data = True,

    entry_points={
        "console_scripts": ["SeqRefFilter = seqreffilter.cli:main"]
    },    

    packages=setuptools.find_packages(),

    install_requires=[
        "numpy>=1.18.1",
        "yxseq>=0.0.1",
        "yxmath>=0.0.3",
        "pysam>=0.21.0",
    ],

    python_requires='>=3.5',
)