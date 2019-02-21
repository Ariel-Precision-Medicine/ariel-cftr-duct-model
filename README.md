# ariel-cftr-duct-model
Mechanistic Modeling CFTR Duct Functionality

This program is a updated translation of the model presented in "A Mathematical Model of the Pancreatic Duct Cell Generating High Bicarbonate Concentrations in Pancreatic Juice" written by David C. Whitcomb, MD, PhD,* and G. Bard Ermentrout, PhD.

Paper can be found online here: https://pdfs.semanticscholar.org/b312/3f29dbb27090bac662e0fc6bb452ed450a79.pdf

Model has been translated from XPPAUT to Python for flexibility and future use.

Variant functionality has not been rigorously verified, but source is Cutting Paper for residual CFTR fxn found here: https://www.ncbi.nlm.nih.gov/pubmed/29805046

Dependencies:
* Python 3
* NumPy
* SciPy
* Matplotlib
* Pandas
* PyQt5

Run using:
```
python3 duct_model_events.py
```