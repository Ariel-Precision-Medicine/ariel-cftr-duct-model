# ariel-cftr-duct-model
Mechanistic Modeling CFTR Duct Functionality

This program is a updated translation of the model presented in "A Mathematical Model of the Pancreatic Duct Cell Generating High Bicarbonate Concentrations in Pancreatic Juice" written by David C. Whitcomb, MD, PhD,* and G. Bard Ermentrout, PhD.

Paper can be found online here: https://pdfs.semanticscholar.org/b312/3f29dbb27090bac662e0fc6bb452ed450a79.pdf

Model has been translated from XPPAUT to Python for flexibility and future use.

Variant functionality has not been rigorously verified, but source is Cutting Paper for residual CFTR fxn found here: https://www.ncbi.nlm.nih.gov/pubmed/29805046

Currently, two display options are available:
(1) a user-interactive GUI application to load in desired variants and see their impacts
(2) HTML-based plotting via Bokeh that outputs cleaner, visually appealing plots

Future work will add interactivity to (2) to allow the user to change criteria for the model (similar to (1)).

All scripts are found in the "scripts" directory.
All image and HTML outputs are found in the "outputs" subdirectory.
PDFs of all supporting papers can be found in the "supporting_papers" directory.

For best performance when using (2) to generate HTML plots, please install the Ariel-specific 'Gilroy' font packages found in the "fonts" directory.

Dependencies:
* Python 3
* NumPy
* SciPy
* Matplotlib
* Pandas
* PyQt5
* Bokeh

Run (1) using:
```
python3 duct_model_events.py
```

Run (2) using:
```
python3 oop_duct_model.py
```