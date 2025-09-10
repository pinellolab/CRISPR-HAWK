from setuptools import setup, find_packages
import os

# Read README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [
        line.strip() for line in fh if line.strip() and not line.startswith("#")
    ]

setup(
    name="crisprhawk",
    version="0.1.0",
    author="Manuel Tognon",
    author_email="manu.tognon@gmail.com",
    description="CRISPR-HAWK: Haplotype and vAriant-aWare guide design toolKit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ManuelTgn/CRISPR-HAWK",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": ["pytest", "black", "flake8", "mypy"],
        "docs": ["sphinx", "sphinx-rtd-theme"],
    },
    entry_points={
        "console_scripts": [
            "crisprhawk=crisprhawk.__main__:main",
        ],
    },
    include_package_data=True,
    package_data={
        "crisprhawk": [
            "scripts/*",
        ],
        "crisprhawk.scores.azimuth": [
            "saved_models/*.pickle",
            "azure_models/*.pickle",
            "data/*",
            "tests/*",
            "LICENSE.txt",
            "*.csv",
            "*.zip"
        ],
        "crisprhawk.scores.cfdscore": [
            "models/*.pkl",
            "*.zip"
        ],
        "crisprhawk.scores.deepCpf1": [
            "weights/*.h5",
            "*.zip"
        ],
        "crisprhawk.scores.elevation": [
            "*.zip",
        ],
        "crisprhawk.config": [
            "*.json",
        ],
    },
    zip_safe=False,
)
