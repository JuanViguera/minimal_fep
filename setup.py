from setuptools import setup, find_packages

setup(
    name="minimal_fep",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[],
    author="Juan Viguera",
    description="Minimal FEP package for relative binding free energy calculations.",
    include_package_data=True,
    python_requires='>=3.7',
)
