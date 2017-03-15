

"""
libGeneral
Various useful routines
"""


import time


class EmptyObject: 
    """
    Define an empty class for passing arguments - if you need to pass multiple 
    arguments to a callback function, just create one of these and add 
    the properties needed.
    """
    pass


def now():
    "returns a string with today's date and time in a nice format"
    return time.strftime('%a %Y-%m-%d, %I:%M:%S %p')


def pause():
    "waits for user to hit enter before continuing"
    print
    raw_input('hit enter:')
    print


def wait(secs=4):
    time.sleep(secs)


def showvars(*vars):
    import sys
    try:
        1/0
    except ZeroDivisionError:
        callers_globals = sys.exc_info()[2].tb_frame.f_back.f_globals
        callers_locals = sys.exc_info()[2].tb_frame.f_back.f_locals
    varids=[]
    for eachvar in vars:
        varids.append(id(eachvar))
    for varname, varvalue in callers_locals.items():
        if id(varvalue) in varids:
            print varname, "=", varvalue
    for varname, varvalue in callers_globals.items():
        if id(varvalue) in varids:
            print varname, "=", varvalue



if __name__=='__main__':
    x=5
    y=6
    showvars(x,y)
    
    
## # save data to file for plotting
## filename = 'hw5epsilon.txt'
## print 'saving data to ', filename
## file = open(filename,"w")
## for i in range(npoints+1):
    ## print >> file, aRho0inv[i], aEpsilon[i]
## file.close()

## filename = 'hw5psi.txt'
## file = open(filename,"w")
## for i in range(npoints):
    ## print >> file, aRho[i], aR[i]
## file.close()
