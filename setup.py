from distutils.core import setup
setup(name='crabs',
      description='CRABS: Creating Reference databases for Amplicon-Based Sequencing',
      author='Gert-Jan Jeunen',
	  author_email='gjeunen@gmail.com',
	  url='https://github.com/gjeunen/reference_database_creator',
	  version='1.12.1',
      packages=['function'],
      scripts=['crabs']
)
