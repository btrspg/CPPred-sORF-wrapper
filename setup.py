from setuptools import setup

setup(
    name='CPPredsorf-wrapper',
    version='0.1.0',
    packages=['cppredsorf'],
    scripts=['CPPred-sORF.py'],
    include_package_data=True,
    package_data={'': ['Model/*.*']},
    url='',
    license='',
    author='CHEN Yuelong',
    author_email='yuelong.chen.btr@gmail.com',
    description='''
This is a wrapper of CPPred-sORF algorithm. Download from original source, [http://www.rnabinding.com/CPPred-sORF/CPPred-sORF/CPPred-sORF.tar.gz](http://www.rnabinding.com/CPPred-sORF/CPPred-sORF/CPPred-sORF.tar.gz), would be perfect.

This wrapper is just for myself easy use.

Almost whole codes were from [http://www.rnabinding.com/CPPred-sORF/CPPred-sORF/CPPred-sORF.tar.gz](http://www.rnabinding.com/CPPred-sORF/CPPred-sORF/CPPred-sORF.tar.gz). I just modified the structure and install setup.py for myself easier use.
    '''
)
