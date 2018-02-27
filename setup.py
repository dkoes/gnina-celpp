import setuptools

setuptools.setup(
    name="gnina-celpp",
    version="0.1.0",
    url="https://github.com/dkoes/gnina-celpp",

    author="David Koes",
    author_email="dkoes@pitt.edu",

    description="gnina contestant for D3R CELPP competition",
    long_description=open('README.rst').read(),

    packages=setuptools.find_packages(),

    install_requires=["d3r"],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
     scripts = ['gnina-celpp/gnina-celpp_dock.py',
                'gnina-celpp/gnina-celpp_ligand_prep.py', 
                'gnina-celpp/gnina-celpp_protein_prep.py']
)
