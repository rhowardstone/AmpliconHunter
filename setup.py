from setuptools import setup, find_packages

setup(
    name="ampliconhunter",
    version="1.0.0",
    description="A scalable tool for PCR amplicon prediction from microbiome samples",
    author="Rye Howard-Stone",
    author_email="rye.howard-stone@uconn.edu",
    url="https://github.com/rhowardstone/AmpliconHunter",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.78",
        "numpy>=1.19.0",
        "pandas>=1.0.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "hyperscan>=0.7.0",
    ],
    entry_points={
        'console_scripts': [
            'ampliconhunter=ampliconhunter.AmpliconHunter:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)