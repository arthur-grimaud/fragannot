python setup.py sdist bdist_wheel
pip install -U --find-links=dist\ fragannot
python clean.py