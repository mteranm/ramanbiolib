# Custom script to build the html docs using sphinx
# To have the updated changes in code is necesary to install the package in dev mode pip install -e . 

sphinx-build -M html . _build/ -W -a -j auto -n --keep-going