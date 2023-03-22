from setuptools import setup, find_packages

setup(
    name='sciso',
    version='1.1',
    description='Analysis of single-cell isotope tracing data',
    url='https://github.com/Buglakova/13C-SpaceM',
    packages=find_packages(include=['sciso']),
    python_requires='>=3.6',
    install_requires=[],
    author='Elena Buglakova',
    author_email='elena.buglakova@embl.de',
    license='MIT'
)
