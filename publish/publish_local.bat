python setup.py sdist bdist_wheel
pip install --find-links=dist\ fragannot
python clean.py