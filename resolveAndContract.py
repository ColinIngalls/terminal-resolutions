# %load Resolve.py
# Resolve.py
"""
We consider a normal order ramified on snc but no secondary ramification. We
blow up repeatedly at nodes recording the log discrepancy and
b-discrepancy until we have all possible negative discrepancy curves.

The exceptional curves should probably be indexed by dyadic fractions, written
as binary numbers, but we will instead clear denominators and work with
integers. 
"""

from fractions import Fraction 
from pprint import pprint



# Enter initial ramification data.
# The ramification cover of the curve i is gi copies of zi modulo n
#initial_ram_data = input("Please enter the ramification data tuple n,z1,z2,g1,g2\n")
#n_str, z1_str, z2_str, g1_str, g2_str = initial_ram_data.split(',')
#n = int(n_str.strip())
#z1 = int(z1_str.strip())
#z2 = int(z2_str.strip())
#g1 = int(g1_str.strip())
#g2 = int(g2_str.strip())


#This function computes the order of an element in a cyclic group.
def order(num,modulus):
    ord = 1
    while ord*num %modulus !=0:
        ord +=1
    return ord


#print("Corresponding maximal order has ram data", end=':')
#print(z1 % n, z2 % n, ' modulo ', n, end='')
#f1 = order(z1,n) # these are ramification indices of max order
#f2 = order(z2,n)
#print(' giving ramification indices ', f1, f2,'.')
# Check working by showing ramification of containing maximal order



"""
This is the function that computes log and b-discrepancies of
all exceptional curves in a good log resolution. Here the good log
resolution is one chosen so all negative discrepancy curves are
guaranteed to be found there. The data will be recorded in a dictionary.
"""
def resolve(n,z1,z2,g1,g2):
    f1 = order(z1,n) # these are ramification indices of max order
    f2 = order(z2,n)
    no_blowups = max(g1,g2)*n -1
    # crude bound on number of blowups to ensure good log resolution
    # print(f'We need to blowup set of nodes  at most {no_blowups} times.')
    no_curves = 2**no_blowups # number of exc curves actually 2**no_blowups -1
    curves = {0:[Fraction(1,f1*g1)-1,z1] , no_curves:[Fraction(1,f2*g2)-1,z2]}
    # Create dictionary of curve data. Keys are 0 up to no_curves.
    # Current values are lists of [log discrepancy, ram of max order] 
    step = no_curves # Final indices for strict transforms of initial ram 
                     # curves are `step' apart

    
    for i in range(1,no_blowups+1):

 #       print(f'This is the {i}-th blowup. New curves have  log and b- discrepancy and ramification:')
        l = len(curves)-1 # the number of new exceptional curves
       
        oldcurves = set(curves.keys()).difference({0,no_curves})  # Gives indices for exceptionals in previous iteration
        for index in oldcurves:
            curves[index] = [curves[index][0], curves[index][1], curves[index][2],curves[index][3]-2]                                                         
        step = int(step/2) # indices of i-th blowup curves have indices 'step' apart
        
        is_good_res = True
        # We'll use this to test if the i-th blowup gives a good resolution

        for j in range(l):
            log_disc = 1 + curves[2*j*step][0] + curves[(2*j+2)*step][0]
            #  Give log discrepancy of j-th new exceptional in i-th blowup
            is_good_res =  is_good_res and (log_disc >=0)
            ram = (curves[2*j*step][1] + curves[(2*j+2)*step][1])%n
            # Gives ramification along j-th new exceptional in i-th blowup
            b_disc = log_disc + 1-Fraction(1,order(ram,n))
            # Compute the b-discrepancy
            curves[(2*j +1)*step] = [log_disc, ram, b_disc,-1]
            # Update curves dictionary with entry for j-th curve in i-th blowup
#            print(log_disc,b_disc, ram, end=';')

        if is_good_res:
            break
            #print(' DONE', end='\n\n')
    indices = sorted(curves)  
    curvesList = [ curves[i] for i in indices] # Create list from dictionary so not dyadic indices.
    curvesList[0].extend([0,0])
    curvesList[len(curvesList)-1].extend([0,0])
    #print(curvesList)    
    return curvesList
#    print('\n\n Good log resolution has been achieved.')

def toNegativeFractions(sequence):
    return [-Fraction(x,1) for x in sequence]

#This functions inputs a sequence of integers appearing 
#in a continued fraction and spits out the actual fraction. 
def HJcontinuedFraction(sequence):
#    print(sequence)
    if len(sequence) == 1:
#        print(sequence[0])
        return sequence[0]
    else:
        first = sequence.pop(0)
        return first-1/HJcontinuedFraction(sequence)
    
# Input sequence of negatives of self-intersections of 
# string of exceptionals, returns cyclic group action
# (1,b)/r
def cyclicGroup(seq):
    frac = HJcontinuedFraction(toNegativeFractions(seq))
    r = frac.numerator
    b = frac.denominator
    return([r,1,b])

# Input list of curves with ramification & 
# self-intersection data. We contract (-1)-curves
# with negative b-disc. I'm worried that on contracting 
# the i-th exceptional, the (i-1)-st exceptional becomes 
# a (-1)-curve which should be contracted but misses out
# as i only increases. 
def contract2smooth(curves):
    noCurves=len(curves)
    i = 1
    while i < noCurves:
        if curves[i][2] >= 0 and curves[i][3] == -1:
#            print(i)
            noCurves=len(curves)
            curves[min(i+1,noCurves-2)][3] += 1
            # Confused about use of min above and max below.
            curves[max(i-1,1)][3] += 1
            curves.pop(i)
     #       print([x[3] for x in curves])        
            i = 1
            noCurves -= 1
        else:
            i += 1
    return curves


def contract2min(curves):
    # scan until you find a b des >= 0
    # then scan until you have b des < 0
    noCurves = len(curves)
    smoothPt = [0]
    for b in range (0,noCurves):
        curves.insert(b*2,smoothPt)
    curves.pop(0) 
    # this is a sequence for data for curves and points
    # print(curves)
    i = 1
    while i < noCurves-1:        
        # print(curves)
        # print(noCurves)
        # print(i)
        if curves[2*i][2] >= 0:
            start = i
            j = i+1
            while j < noCurves and curves[2*j][2] >= 0:
                j += 1
            stop = j
            # print('start stop are:',start,stop)
            # print('to contract are:',[curves[2*x] for x in range(start,stop)])
            cg = cyclicGroup([curves[2*x][3] for x in range(start,stop)])
            # print(cg)
            # print(curves)
            del curves[2*start-1:2*stop]
            # print(curves)
            curves.insert(2*start-1,cg)
            # print(curves)
        i += 1    
        noCurves = (len(curves)-1)/2
    return(curves)

def frPrint(list):
    [print(*x) for x in list]
    return

def resolveAndContract(n,z1,z2,g1,g2):
    resolved = resolve(n,z1,z2,g1,g2)
    #print('resolution')
    #frPrint(resolved)
    smth = contract2smooth(resolved)
    #print('contracted to smooth')
    #frPrint(smth)
    minTerminRes = contract2min(smth)
    print('minimal Terminal resoluion')
    frPrint(minTerminRes)
    return(minTerminRes)


# running some sample output
primes  = [2,3,5,7,11,13,17,19,23]
primes = [3,5,7]
for n in primes:
    print('input data',n,0,0,n,n)
    resolveAndContract(n,0,0,n,n)
    for z1 in range(n):
        print('input data',n,z1,0,1,n)
        resolveAndContract(n,z1,0,1,n)
    for z1 in range(n):
        for z2 in range(z1+1):
            print('input data',n,z1,z2,1,1)  
            resolveAndContract(n,z1,z2,1,1)

            
