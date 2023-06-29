from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()


setup(
   name='matrex',
   version='0.0.1',
   description='A package grouping matrices row reordering algorithms to minimize their front size.',
   license="MIT License",
   long_description=long_description,
   long_description_content_type="text/markdown",
   author=["benjaminlanthier", "cotejer"],
   author_email=['benjamin.lanthier2@usherbrooke.ca', 'jeremy.cote6@usherbrooke.ca'],
   url="https://github.com/benjaminlanthier/MatRexAlgs/",
   packages=['MatRexAlgs'],
   classifiers=[
       "Programming Language :: Python :: 3",
       "License :: OSI Approved :: MIT License",
       "Operating System :: OS Independent",
   ],

   install_requires=[
       'numpy>=1.23.3',
       'networkx>=3.1',
   ],

   python_requires='>=3.9'
)