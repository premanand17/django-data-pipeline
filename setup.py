import os
from setuptools import setup, find_packages

# with open(os.path.join(os.path.dirname(__file__), 'README.rst')) as readme:
#     README = readme.read()

# allow setup.py to be run from any path
os.chdir(os.path.normpath(os.path.join(os.path.abspath(__file__), os.pardir)))
ROOT = os.path.abspath(os.path.dirname(__file__))

setup(
    name='data-pipeline',
    version='0.0.1',
    packages=find_packages(),
    package_data={'data_pipeline': ['tests/data/DOWNLOAD/INTACT/*zip',
                                    'tests/data/DOWNLOAD/GENE_HISTORY/*gz',
                                    'tests/data/DOWNLOAD/ENSEMBL_GENE_GTF/*gz',
                                    'tests/data/DOWNLOAD/GENE2ENSEMBL/*gz',
                                    'tests/data/DOWNLOAD/GENE_INFO/*gz',
                                    'tests/data/DOWNLOAD/GENE_PUBS/*gz',
                                    'tests/data/DOWNLOAD/DBSNP/*gz',
                                    'tests/data/DOWNLOAD/RSMERGEARCH/*gz',
                                    'tests/data/DOWNLOAD/DISEASE/*txt',
                                    'tests/*.ini'], },
    include_package_data=True,
    zip_safe=False,
    url='http://github.com/D-I-L/django-data-pipeline',
    description='A data pipeline app.',
    long_description=open(os.path.join(ROOT, 'README.rst')).read(),
    install_requires=["requests>=2.7.0", "Django>=1.8.2", "ftputil>=3.2"],
    classifiers=[
        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.4',
        'Topic :: Internet :: WWW/HTTP',
        'Topic :: Internet :: WWW/HTTP :: Dynamic Content',
    ],
)
