'''
Created on Oct 1, 2010
@author: areich
'''

import csv

#from uuid import uuid4
#def uniqueID(prefix, suffix=None):
#    return uuid4().__str__()

__keyCounts = {}
def uniqueID(prefix, suffix=None):
    result = __keyCounts.setdefault(prefix, 0)
    __keyCounts[prefix] = result + 1
    result = prefix + "-Key-" + str(result)
    if suffix:
        result = result + suffix
    return result


class PredicateException(Exception):
    '''Raise for errors wrt Predicate domains and ranges.'''
    pass


def identity(arg):
    return arg


class DDFGraph():

    def __init__(self):
        self.__typeConverters = {"Float":float, "Int":int, "String":identity} # type converters
        # TODO: Add ability to 'register' new type converters
        self.__conceptSignTypes = {} # Sign types, indexed by Concept key
        self.__concepts = {} # Concept names indexed by Concept key
        self.__conceptKeys = {} # Concept keys indexed by Concept name
        self.__terms = {} # Term tuples indexed by Term key
        self.__termKeys = {} # Term keys indexed by Term tuple
        self.__predicates = {} # Predicate names indexed by Predicate key
        self.__predicateKeys = {} # Predicate key indexed by Predicate name
        self.__predicateDomains = {} # Predicate domain names indexed by Predicate key
        self.__predicateRanges = {} # Predicate range names indexed by Predicate key
        # TODO: Fix how I'm counting terms and statements. 
        self.__numberOfStatements = 0 # Tracked for convenience
        self.__numberOfTerms = 0 # Ditto
        # Term indexes:
        #                FIRST INDEX    SET OF ITEMS
        self.__cs = {} #  concept         sign
        self.__sc = {} #  sign            concept
        # Statement indexes:
        #                 FIRST INDEX   2ND INDEX    SET OF ITEMS
        self.__spo = {} # subjectTerm   predicate    objectTerm
        self.__pos = {} # predicate     objectTerm   subjectTerm
        self.__osp = {} # objectTerm    subjectTerm  predicate

    def concept(self, conc, signType="String"):
        '''Create a Concept in the graph. Returns the Concept's key.
        If the Concept already exists, then the existing key is returned.'''
        conceptKey = None
        if conc not in self.__conceptKeys:
            conceptKey = uniqueID("Concept")
            self.__concepts[conceptKey] = conc
            self.__conceptKeys[conc] = conceptKey
            self.__conceptSignTypes[conceptKey] = signType
        else:
            conceptKey = self.__conceptKeys[conc]
        return conceptKey

    def term(self, conceptKey, sign):
        '''Create a Term, as a tuple, in the graph.  Returns the Term's key.
        If the Term already exists, then the existing key is returned.'''
        self.__addToTermIndex(self.__cs, conceptKey, sign)
        self.__addToTermIndex(self.__sc, sign, conceptKey)
        self.__numberOfTerms += 1
        #term = (sign, conceptKey) # Make Term a tuple
        term = (conceptKey, sign) # Make Term a tuple
        termKey = None
        if term not in self.__termKeys:
            termKey = uniqueID("Term")
            self.__terms[termKey] = term
            self.__termKeys[term] = termKey
        else:
            termKey = self.__termKeys[term]
        return termKey

    def predicate(self, pred, domainConceptKey, rangeConceptKey):
        '''Create a Predicate in the graph.'''
        predKey = None
        # If new predicate, index it and its domain and range concepts
        if pred not in self.__predicateKeys:
            predKey = uniqueID("Predicate")
            self.__predicates[predKey] = pred
            self.__predicateKeys[pred] = predKey
            self.__predicateDomains[predKey] = domainConceptKey
            self.__predicateRanges[predKey] = rangeConceptKey
        else:
            predKey = self.__predicateKeys[pred]
            dcKey = self.__concepts(self.__predicateDomains[predKey])
            rcKey = self.__concepts(self.__predicateRanges[predKey])
            if not (dcKey == domainConceptKey):
                pass # TODO: Throw exception
            if not (rcKey == rangeConceptKey):
                pass # TODO: Throw exception
        return predKey

    def getTerm(self, termKey):
        '''Get the Term indexed by the given key.'''
        return self.__terms[termKey]

    def getConcept(self, conceptKey):
        '''Get the Concept indexed by the given key.'''
        return self.__concepts[conceptKey]

    def getPredicate(self, predKey):
        '''Get the Predicate indexed by the given key.'''
        return self.__predicates[predKey]

    def getPredicateKey(self, predName):
        '''Get the Predicate key corresponding to the given name'''
        return self.__predicateKeys[predName]

    def statement(self, subjTermKey, predKey, objTermKey):
        '''Create a Statement in the graph.'''
        self.__addToStatementIndex(self.__spo, subjTermKey, predKey, objTermKey)
        self.__addToStatementIndex(self.__pos, predKey, objTermKey, subjTermKey)
        self.__addToStatementIndex(self.__osp, objTermKey, subjTermKey, predKey)
        self.__numberOfStatements += 1

    def __addToStatementIndex(self, index, a, b, c):
        '''Index the elements of a statement.'''
        if a not in index:
            index[a] = {b:set([c])}
        else:
            if b not in index[a]:
                index[a][b] = set([c])
            else:
                index[a][b].add(c)

    def __addToTermIndex(self, index, key, val):
        '''Index the elements of a term.'''
        if key not in index:
            index[key] = set([val])
        else:
            index[key].add(val)

    def removeStatement(self, (sub, pred, obj)):
        '''Remove a Statement from the graph.'''
        stmts = list(self.triples((sub, pred, obj)))
        for (delSub, delPred, delObj) in stmts:
            self.__removeFromStatementIndex(self.__spo, delSub, delPred, delObj)
            self.__removeFromStatementIndex(self.__pos, delPred, delObj, delSub)
            self.__removeFromStatementIndex(self.__osp, delObj, delSub, delPred)
        self.__numberOfStatements -= 1

    def __removeFromStatementIndex(self, index, a, b, c):
        '''Remove a the indexing of a statement.'''
        try:
            bs = index[a]
            cset = bs[b]
            cset.remove(c)
            if len(cset) == 0: del bs[b]
            if len(bs) == 0: del index[a]
        # KeyErrors occur if a term was missing, which means that it wasn't a valid delete:
        except KeyError:
            pass

    def statements(self, (subjConcept, subjSign, pred, objConcept, objSign)):
        '''Generate Statements according to the input pattern.'''
        try:
            for subjConcKey, subjSgn in self.terms((subjConcept, subjSign)):
                for objConcKey, objSgn in self.terms((objConcept, objSign)):
                    predKey = self.__predicateKeys[pred]
                    subjTermKey = self.__termKeys[(subjConcKey, subjSgn)]
                    objTermKey = self.__termKeys[(objConcKey, objSgn)]
                    for su, pr, ob in self.triples((subjTermKey, predKey, objTermKey)):
                        suTerm = self.__terms[su]
                        obTerm = self.__terms[ob]
                        prName = self.__predicates[pr]
                        suConc = self.__concepts[suTerm[0]]
                        obConc = self.__concepts[obTerm[0]]
                        suSign = suTerm[1]
                        obSign = obTerm[1]
                        yield (suConc, suSign, prName, obConc, obSign)
        except KeyError:
            pass

    def terms(self, (conc, sgn)):
        '''Generate Terms according to the input pattern.'''
        try:
            if conc is not None:
                if sgn is not None:
                    if sgn in self.__cs[conc]:
                        yield (conc, sgn)
                else:
                    for s in self.__cs[conc]:
                        yield (conc, s)
            else:
                if sgn is not None:
                    for c in self.__sc[sgn]:
                        yield (c, sgn)
                else:
                    for c, sgnSet in self.__cs.items():
                        for s in sgnSet:
                            yield (c, s)
        except KeyError:
            pass

    def triples(self, (sub, pred, obj)):
        '''Generate triples according to the input pattern.'''
        # check which terms are present in order to use the correct index:
        try:
            if sub is not None:
                if pred is not None:
                    # sub pred obj
                    if obj is not None:
                        if obj in self.__spo[sub][pred]:
                            yield (sub, pred, obj)
                    # sub pred None
                    else:
                        for retObj in self.__spo[sub][pred]:
                            yield (sub, pred, retObj)
                else:
                    # sub None obj
                    if obj is not None:
                        for retPred in self.__osp[obj][sub]:
                            yield (sub, retPred, obj)
                    # sub None None
                    else:
                        for retPred, objSet in self.__spo[sub].items():
                            for retObj in objSet:
                                yield (sub, retPred, retObj)
            else:
                if pred is not None:
                    # None pred obj
                    if obj is not None:
                        for retSub in self.__pos[pred][obj]:
                            yield (retSub, pred, obj)
                    # None pred None
                    else:
                        for retObj, subSet in self.__pos[pred].items():
                            for retSub in subSet:
                                yield (retSub, pred, retObj)
                else:
                    # None None obj
                    if obj is not None:
                        for retSub, predSet in self.__osp[obj].items():
                            for retPred in predSet:
                                yield (retSub, retPred, obj)
                    # None None None
                    else:
                        for retSub, predSet in self.__spo.items():
                            for retPred, objSet in predSet.items():
                                for retObj in objSet:
                                    yield (retSub, retPred, retObj)
        # KeyErrors occur if a query term wasn't in the index, so we yield nothing:
        except KeyError:
            pass

#    def value(self, sub=None, pred=None, obj=None):
#        for retSub, retPred, retObj in self.statements((sub, pred, obj)):
#            if sub is None: return retSub
#            if pred is None: return retPred
#            if obj is None: return retObj
#            break
#        return None

    def load(self, filename):
        print "Loading from %r" % (filename,)
        f = open(filename, "rb")
        reader = csv.reader(f)
        for sub, pred, obj in reader:
            sub = unicode(sub, "UTF-8")
            pred = unicode(pred, "UTF-8")
            obj = unicode(obj, "UTF-8")
            if pred in self.__predicateKeys:
                predKey = self.__predicateKeys[pred]
                domainConceptKey = self.__predicateDomains[predKey]
                rangeConceptKey = self.__predicateRanges[predKey]
                domainConceptSignType = self.__conceptSignTypes[domainConceptKey]
                rangeConceptSignType = self.__conceptSignTypes[rangeConceptKey]
                subTypeConverter = self.__typeConverters[domainConceptSignType]
                objTypeConverter = self.__typeConverters[rangeConceptSignType]
                subjTermKey = self.term(domainConceptKey, subTypeConverter(sub))
                objTermKey = self.term(rangeConceptKey, objTypeConverter(obj))
            else:
                raise PredicateException("No definition for predicate: %s" % (pred,))
            self.statement(subjTermKey, predKey, objTermKey)
        f.close()
        print "Loading completed."

    def save(self, filename):
        f = open(filename, "wb")
        writer = csv.writer(f)
        for subj in self.__spo:
            for pred in self.__spo[subj]:
                for obj in self.__spo[subj][pred]:
                    predName = self.__predicates[pred]
                    subjTerm = self.__terms[subj]
                    subjSign = subjTerm[0]
                    #subjConcept = self.__concepts[subjTerm[1]]
                    objTerm = self.__terms[obj]
                    objSign = objTerm[0]
                    #objConcept = self.__concepts[objTerm[1]]
                    writer.writerow([subjSign.encode("UTF-8"), predName.encode("UTF-8"),
                                     objSign.encode("UTF-8")])
        f.close()

    def query(self, clauses):
        bindings = None
        for clause in clauses:
            bpos = {}
            qc = []
            for x, pos in zip(clause, range(5)):
                if x.startswith('?'):
                    qc.append(None)
                    bpos[x[1:]] = pos
                else:
                    qc.append(x)
            rows = list(self.statements((qc[0], qc[1], qc[2], qc[3], qc[4])))
            if bindings is None:
                bindings = []
                for row in rows:
                    binding = {}
                    for var, pos in bpos.items():
                        binding[var] = row[pos]
                    bindings.append(binding)
            else:
                newb = []
                for binding in bindings:
                    for row in rows:
                        validmatch = True
                        tempbinding = binding.copy()
                        for var, pos in bpos.items():
                            if var in tempbinding:
                                if tempbinding[var] != row[pos]:
                                    validmatch = False
                            else:
                                tempbinding[var] = row[pos]
                        if validmatch: newb.append(tempbinding)
                bindings = newb
        return bindings

    def applyinference(self, rule):
        queries = rule.getqueries()
        bindings = []
        for query in queries:
            bindings += self.query(query)
        for b in bindings:
            new_statements = rule.makeStatements(b)
            for triple in new_statements:
                self.addStatement(triple)

    def printGraph(self):
        print "\n=================================================="
        print "             D D F   G r a p h \n"
        self.printPredicates()
        self.printConcepts()
        self.printTerms()
        self.printStatements()
        self.prettyPrint()
        self.summary()
        print "=================================================="

    def printPredicates(self):
        # TODO: Also print out the number of instances of each predicate
        print "                PREDICATES"
        print "PredKey, PredName, DomainConcept, RangeConcept"
        print "--------------------------------------------------"
        for p in self.__predicates:
            domain = self.__concepts[self.__predicateDomains[p]]
            range = self.__concepts[self.__predicateRanges[p]]
            print "%s, %s, %s, %s" % (p, self.__predicates[p], domain, range)
        print # newline

    def printConcepts(self):
        # TODO: Also print out the number of instances of each concept
        print "                 CONCEPTS"
        print "ConceptKey, ConceptName, ConceptSignType"
        print "--------------------------------------------------"
        for c in self.__concepts:
            print "%s, %s, %s" % (c, self.__concepts[c], self.__conceptSignTypes[c])
        print # newline

    def printTerms(self):
        print "                  TERMS"
        print "TermKey, ConceptKey, Sign"
        print "--------------------------------------------------"
        for t in self.__terms:
            (conceptKey, sign) = self.__terms[t]
            print "%s, %s, %r" % (t, conceptKey, sign)
        print # newline

    def printStatements(self):
        print "               STATEMENTS (RAW)"
        print "SubjectTermKey, PredicateKey, ObjectTermKey"
        print "--------------------------------------------------"
        for subj in self.__spo:
            for pred in self.__spo[subj]:
                for obj in self.__spo[subj][pred]:
                    print "%s, %s, %s" % (subj, pred, obj)
        print # newline

    def prettyPrint(self):
        print "            STATEMENTS (READABLE)"
        print "--------------------------------------------------"
        for subj in self.__spo:
            for pred in self.__spo[subj]:
                for obj in self.__spo[subj][pred]:
                    print self.statementString(subj, pred, obj)
        print # newline

    def statementString(self, subj, pred, obj):
        predName = self.__predicates[pred]
        subjTerm = self.__terms[subj]
        subjSign = subjTerm[1]
        subjConcept = self.__concepts[subjTerm[0]]
        objTerm = self.__terms[obj]
        objSign = objTerm[1]
        objConcept = self.__concepts[objTerm[0]]
        return "(%s:%r) <%s> (%s:%r)" % (subjConcept, subjSign, predName, objConcept, objSign)

    def summary(self):
        print "                GRAPH SUMMARY"
        print "--------------------------------------------------"
        print "# of Concepts: %s" % len(self.__concepts)
        print "# of Predicates: %s" % len(self.__predicates)
        #print "# of Terms: %s (graph vertices)" % len(self.__terms)
        print "# of Terms: %s (graph vertices)" % (self.__numberOfTerms,)
        print "# of Statements: %s (graph edges)" % (self.__numberOfStatements,)

    def printTermIndex(self, index): # DEBUG PRINT
        print "------------------------------------------------"
        print "DEBUG PRINTOUT -- TERM INDEX\n"
        for a in index:
            print "%s" % (a,)
            for b in index[a]:
                print "    %s" % (b,)

    def printTermIndices(self): # DEBUG PRINT
        self.printTermIndex(self.__cs)
        self.printTermIndex(self.__sc)

    def printBindings(self, bindings):
        for binding in bindings:
            for queryVar, termKey in binding.iteritems():
                (sign, conceptKey) = self.getTerm(termKey)
                print "%s = (%s, %s)" % (queryVar, sign, (self.getConcept(conceptKey)))


class InferenceRule:
    def getqueries(self):
        return None
    def makeStatements(self, binding):
        return self.__makeStatements(**binding)


if __name__ == "__main__":
    pass
