'''
Created on Oct 10, 2010

@author: areich
'''

from time import time
from ddfgraph import DDFGraph

g = DDFGraph() # Create a new, empty DDFGraph object

# NOTE: All of the of the values returned from the calls to the concept & predicate
# methods below are 'keys'.

# Add Concept definitions to the DDFGraph, along with their datatype information.
# The default datatype is "String".  In each case, the returned result is the ID
# of the Concept in the DDFGraph.  These IDs are then used, below, to specify
# domain and range Concepts of Predicates.
dollars = g.concept("Dollars", "Float")
contribution = g.concept("Contribution")
person = g.concept("Person")
tradingSymbol = g.concept("TradingSymbol")
location = g.concept("Location")
industry = g.concept("Industry")
company = g.concept("Company")

# Add Predicate definitions to the DDFGraph, along with the domain and range
# concepts of each predicate.  Each call to the predicate method returns the
# new predicate's ID in the DDFGraph, but they are not saved in the expressions,
# below, because there is no use for them (yet) outside of the DDF Graph object
# -- unlike the IDs captured above in the concept expressions.
# 
# predicate( PredicateName, DomainConcept, RangeConcept )
#---------------------------------------------------------
g.predicate("amount", contribution, dollars)
g.predicate("contributor", contribution, tradingSymbol)
g.predicate("director", person, tradingSymbol)
g.predicate("headquarters", tradingSymbol, location)
g.predicate("industry", tradingSymbol, industry)
g.predicate("name", tradingSymbol, company)
g.predicate("recipient", contribution, person)
g.predicate("revenue", tradingSymbol, dollars)

# Now that the Concepts and Predicates are defined, we can load a set of
# instance data.  There are two in this example: The full dataset consists
# of over 36,000 DDF statements and a test dataset with about 100 DDF
# statements.  Both are stored as .csv files.

# Full dataset file and threshold to use for $ amounts reported
#g.load("data/bizDataFULL.csv")
#threshold = 300000 # Dollars

# Test dataset file and threshold to use for $ amounts reported
g.load("data/bizDataTEST.csv")
threshold = 0 # Dollars

# Here's an example query.  It's wrapped with some start/end time stuff so we
# can measure how long it takes.  The query consists of several clauses that
# are implicitely AND'ed.  Variables to be bound in the query are strings that
# start with '?'.  The answers returned consists of a list of dictionaries
startTime = time() # Record the start time of the query
answers = g.query([(contribution,'?contrib','recipient',person,'?person'),
                   (contribution,'?contrib','amount',dollars,'?dollars')])
queryRunTime = time() - startTime
print "\nQuery runtime: %s seconds\n" % (queryRunTime,)

# Uncomment the code below to see what the answers returned look like; but
# don't do this for the full dataset.  Or just look at what is produced
# by looping through the answers in the next section, below.
print "Answers:"
print answers

# This code will loop through the answers returned by the query and print
# results for contributions over the threshold amount.
print "\nThe data set contains %r contributions.  Here are the ones over $%r:" % (len(answers),threshold)
for answer in answers:
    if answer['dollars'] > threshold:
        print '%r received a contribution of $%r' % (answer['person'],answer['dollars'])

g.printGraph() # Only uncomment this when using the test dataset; the full dataset is too large

#-------------------------------------------------------------------------
# Resulting Output from Test Dataset.  See the file, 'example1_output_with_printgraph.py'
# for an example of the test dataset output that has the printGraph() statement, above,
# uncommented.
# 
#Loading from 'testData1.csv'
#Loading completed.
#
#Query runtime: 0.00300002098083 seconds
#
#Answers:
#[{'contrib': u'contrib252', 'person': u'Dave Camp', 'dollars': 5000.0},
#{'contrib': u'contrib250', 'person': u'John Campbell', 'dollars': 5000.0},
#{'contrib': u'contrib251', 'person': u'Jerry Moran', 'dollars': 3250.0},
#{'contrib': u'contrib1556', 'person': u'Pete King', 'dollars': 5000.0}]
#
#The data set contains 4 contributions.  Here are the ones over $0:
#u'Dave Camp' received a contribution of $5000.0
#u'John Campbell' received a contribution of $5000.0
#u'Jerry Moran' received a contribution of $3250.0
#u'Pete King' received a contribution of $5000.0
#-------------------------------------------------------------------------
