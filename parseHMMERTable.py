'''
Created on Aug 2, 2013

@author: ranko
'''

''' 
parses HMMER results TABLE
NOTE: 
 -> so far doesn't parse domains
    or alignments!    
 -> hmmer MUST be called with --notextw
      
Other Notes:
 -> hmmer can have multiple queries, and multiple iterations per query, 
 therefore parsed results are in format:
  qry<1+> -> iteration<1+> -> results<1+> 

MODEL
   f: parseHMMResults(filePath):    --> parses HMMER results (input: filename)
                                     --> returns: multiQueryResults (see below)
multiQueryResults
   -> qry <1 or more>                     {HMMER results for xth query sequence}
       -> hmmerType <1>                   {type of HMMER used (jack or normal) }
       -> queryFile <1>                   {file that was used to query HMM db}
       -> resultsFile <1>                 {results file for HMMER run}
       -> targetDB <1>                    {HMM db that was queried}
       -> qryNr    <1>                    {number of query, as entered}
       -> iterationResults <1 or more>    {all results for iterations}
            -> iteration                  {simply number of iteration}
                f: getLastIteration --> returns last iteration
                f: getIterationsNr  --> returns # of iterations 
            -> hmmTable                   {tabular results}
                    -> seqID
                    -> fullSeqEV
                    -> fullSeqBS
                    -> fullSeqBSBias
                    -> bestDomEV
                    -> bestDomBS
                    -> bestDomBSBias
                    -> nrDom
                    -> nrDomExp
                    -> seqID
                    -> description
#NYI            -> domainResults    #NYI

#TODO: write rest of parser!
'''

class HmmResultsMQuery:
    def __init__(self):
        self.qry = []

class HmmTableResult:
    def __init__(self):
        self.fullSeqEV = 0.0
        self.fullSeqBS = 0.0
        self.fullSeqBSBias = 0.0
        self.bestDomEV = 0.0
        self.bestDomBS = 0.0
        self.bestDomBSBias = 0.0
        self.nrDom = 0.0
        self.nrDomExp = 0.0
        self.seqID = ""
        self.description = ""
        
class HmmResultsOneQuery:
    def __init__(self):
        self.iterationResults = []
        self.queryFile = ""
        self.resultsFile = ""
        self.targetDB = ""
        self.hmmerType = ""
        self.qryNr = 0  
    def getIterationsNr(self):
        return len(self.iterationResults)
    def getLastIteration(self):
        return self.iterationResults[self.getIterationsNr()-1]

class HmmResultsOneIteration:
    def __init__(self):
        self.hmmTable = []
        self.iteration = 0
    
def parseHMMResults (hmmRes):
#    print "parsing..."
#    hmmRes = '/home/ranko/workspace_e4/PyParsers/Results/HMM_RESULTS/run3_23.08/bp_nrf2_hs_bzip_nr_ref_e0.01_c_0.5_BAC.jhmr'
    with open(hmmRes) as hmmFile:
        # file lines
        hmmFL = hmmFile.readlines()
        # all results
        hmmAllResults = HmmResultsMQuery()     
        # tmp results (one query)
        tmpHmmResultsOneQuery = HmmResultsOneQuery()    
        # tmp table (one table for one set of results)    
        tmpHmmTable = []
        # tmp one iteration
        tmpHmmerOneIteration = HmmResultsOneIteration()    
        # for decisions on hmmfullSeqEer type parsing
        hmmerType = ""    
        # PARSE HEADER    
        doParseHeader = True
        # PARSE TABLE
        doParseTable = False        
        startParseTable = False
        # PARSE QUERY    
        queryNR = 0
        newQuery = True
        newIteration = True
        # for managing jackHmmer
        jackHmmerIteration = 0
    
        for line in hmmFL:
            # catch new query
            if '#' not in line and '//' in line:
                newIteration = True
                newQuery = True
            
            # catch jackhmmer iteration
            if 'jackhmmer' in hmmerType and '@@ Round:' in line:                        
                newIteration = True
        
            if newIteration:
                doParseTable = False        
                startParseTable = False
                tmpHmmTable = []
                tmpHmmerOneIteration.iteration = jackHmmerIteration
                if jackHmmerIteration > 0:
                    tmpHmmResultsOneQuery.iterationResults.append(tmpHmmerOneIteration)            
                tmpHmmerOneIteration = HmmResultsOneIteration()
                jackHmmerIteration +=1            
                doParseHeader = True
                # PARSE TABLE
                doParseTable = False                    
                startParseTable = False        
                newIteration = False
          
            if newQuery:
                # initialize            
                doParseHeader = True
                # PARSE TABLE
                doParseTable = False        
                startParseTable = False
                if queryNR >= 1:
                    hmmAllResults.qry.append(tmpHmmResultsOneQuery)
                queryNR += 1
                newQuery = False
                tmpHmmResultsOneQuery = HmmResultsOneQuery()
                tmpHmmerOneIteration = HmmResultsOneIteration()                                        
                
            if doParseHeader:
                if "# hmmsearch :: " in line:
                    hmmerType = "hmmsearch"
                    tmpHmmResultsOneQuery.hmmType = "hmmsearch"
                if "# jackhmmer :: " in line:
                    hmmerType = "jackhmmer"         
                    tmpHmmResultsOneQuery.hmmType = "jackhmmer"                    
                if "# query HMM file:" in line or '# query sequence file:' in line: 
                    tmpHmmResultsOneQuery.queryFile = line[35:].strip()
                if "# target sequence database:" in line: 
                    tmpHmmResultsOneQuery.targetDB = line[35:].strip()
                if "# output directed to file:" in line: 
                    tmpHmmResultsOneQuery.resultsFile = line[35:].strip()
                if "Scores for complete sequences" in line:
                    doParseHeader = False
                    doParseTable = True
                    startParseTable = False
                
            elif doParseTable:
                line = line.strip()
                if len(line) > 0:
                    if line[0] == '-':
                        line = ' '+line[1:]
                    if line[0] == '+':
                        line = ' '+line[1:]
            # end of table
                if line.strip() == '':
                    startParseTable = False
                    doParseTable = False 
                    tmpHmmerOneIteration.hmmTable = tmpHmmTable  
                                                     
                if startParseTable:
                    line = line.strip()
                
                    if 'inclusion threshold' not in line:
                        while '  ' in line:                    
                            line = line.replace('  ',' ')
                        
                        sl = line.split(' ')
                        # RECORD RESULT
                        tabResult = HmmTableResult()
                        # FULL SEQ EV                            
                        tabResult.fullSeqEV = float(sl[0].strip())
                        # FULL SEQ SCORE
                        tabResult.fullSeqBS = float(sl[1].strip())
                        # FULL SEQ BIT SCORE BIAS
                        tabResult.fullSeqBSBias = float(sl[2].strip())
                    
                        # BEST DOMAIN EV
                        tabResult.bestDomEV = float(sl[3].strip())
                        # BEST DOMAIN BITSCORE
                        tabResult.bestDomBS = float(sl[4].strip())
                        # BEST DOMAIN BITSCORE BIAS
                        tabResult.bestDomBSBias = float(sl[5].strip())
                    
                        # NR DOMAINS
                        tabResult.nrDom = float(sl[6].strip())
                        # EXPECTED NR DOMAINS
                        tabResult.nrDomExp = float(sl[7].strip())
                
                        # SEQ ID
                        tabResult.seqID = sl[8].strip()
                
                        # DESCRIPTION (OPTIONAL)
                        tabResult.description = ""
                        if len(sl) > 9:
                            d = ''
                            for desL in range (9,len(sl)):
                                d = d + ' '+sl[desL] 
                            tabResult.description = d
#                            print tabResult.description
                
                        # SAVE IT TO TABLE
                        #print tabResult.seqID,tabResult.fullSeqEV,"added"                
                        tmpHmmTable.append(tabResult)            
                                
                if "    ------- ------ -----" in line: 
                    startParseTable = True            
            
    # now parse results table                                
    return hmmAllResults
# ------------------------ END OF PARSE MODULE ------------------
            
# MAIN: TEST
'''
target = '/home/ranko/Development/workspace/TaxMapper/src/testJACK.jhmr'
hmmAllResults = parseHMMResults(target)              
print " --- ALL QUERIES --- "
cntQ = 0
for hmmQR in hmmAllResults.qry:
    cntQ +=1
    print " **** QUERY ",cntQ," ****"
    print "res file: ", hmmQR.resultsFile                
    print "target DB: ", hmmQR.targetDB                
    print "query file:", hmmQR.queryFile
    print "iterations: ", len(hmmQR.iterationResults)    
    for hmmQI in hmmQR.iterationResults:        
        print " #### ITERATION: ", hmmQI.iteration," ####"        
        print " @@@@ RESULTS TABLE @@@@ "        
        cnt = 0
        for h in hmmQI.hmmTable:
            cnt +=1
            if cnt < 50:
                print cnt, h.seqID, "eV: ",h.fullSeqEV        
                
print ' ----------------------------------------------------- '            
print ' ----------------------------------------------------- '
print ' ----------------------------------------------------- '
cntQ = 0
for hmmQR in hmmAllResults.qry:
    cntQ +=1
    print " **** QUERY ",cntQ," ****"
    print "res file: ", hmmQR.resultsFile                
    print "target DB: ", hmmQR.targetDB                
    print "query file:", hmmQR.queryFile
    print "iterations: ", len(hmmQR.iterationResults)    
    
    hmmQI = hmmQR.getLastIteration()        
    print " #### ITERATION: ", hmmQI.iteration," ####"        
    print " @@@@ RESULTS TABLE @@@@ "        
    cnt = 0
    for h in hmmQI.hmmTable:
        cnt +=1
        if cnt < 50:
            print cnt, h.seqID, "eV: ",h.fullSeqEV                                        
'''            