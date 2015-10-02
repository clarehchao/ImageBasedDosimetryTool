import json
import pandas as pd

def LoadInputParameter(fname):
    """
    load the input file and read the parameters into a dictionary via JSON
    """
    thedict = {}
    jdec = json.JSONDecoder()
    ff = open(fname)
    for line in ff:
        thedict = jdec.decode(line)
    return thedict 
    
    
def Dict2DF(thedict,thecolname,fname):
    """
       convert a dictionary to a data frame for file saving
    """
    df = pd.DataFrame(thedict.values(),index=thedict.keys(),columns=thecolname)
    df.to_csv(fname)
    print 'Save the dataframe to {}'.format(fname)
    return df
