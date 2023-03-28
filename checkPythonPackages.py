print ("")

try:
    import MySQLdb
    print ('MySQLdb exist')
except ImportError:
    print ('please install MySQLdb: yum install MySQL-python')

try:
    import sqlalchemy
    print ('sqlalchemy exist')
except ImportError:
    print ('please install sqlalchemy: pip install sqlalchemy==1.1')

try:
    import pympler
    print ('pympler exist')
except ImportError:
    print ('please install pympler: pip install pympler==0.5')

try:
    import numpy
    print ('numpy exist')
except ImportError:
    print ('please install numpy: yum install numpy==1.7.1')

try:
    import networkx
    print ('networkx exist')
except ImportError:
    print ('please install networkx: pip install networkx==1.11')
		
try:
    import pymzml
    print ('pymzml exist')
except ImportError:
    print ('please install pymzml: pip install pymzml==0.7.8')

try:
    import pyteomics
    print ('pyteomics exist')
except ImportError:
    print ('please install pyteomics: pip install pyteomics==3.5.1')

print ("")
