from setuptools import find_packages, setup
setup(
    name='cosmolenslib',
    packages=find_packages(include=['cosmolenslib']),
    version='0.1.0',
    description='my local python lib',
    author='Carlo Giocoli',
    license='INAF - OAS',
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)
