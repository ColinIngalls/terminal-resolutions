# Resolve.py
"""
We consider a normal order ramified on snc but no secondary ramification. We
blow up repeatedly at nodes recording the log discrepancy and
b-discrepancy until we have all possible negative discrepancy curves.

The exceptional curves should probably be indexed by dyadic fractions, written
as binary numbers, but we will instead clear denominators and work with
integers. 
"""
# Enter initial ramification data.
# The ramification cover of the curve i is gi copies of zi modulo n
initial_ram_data = input("Please enter the ramification data tuple n,z1,z2,g1,g2\n")
n_str, z1_str, z2_str, g1_str, g2_str = initial_ram_data.split(',')
n = int(n_str.strip())
z1 = int(z1_str.strip())
z2 = int(z2_str.strip())
g1 = int(g1_str.strip())
g2 = int(g2_str.strip())


#This function computes the order of an element in a cyclic group.
def order(num,modulus):
    ord = 1
    while ord*num %modulus !=0:
        ord +=1
    return ord


print("Corresponding maximal order has ram data", end=':')
print(z1 % n, z2 % n, ' modulo ', n, end='')
f1 = order(z1,n) # these are ramification indices of max order
f2 = order(z2,n)
print(' giving ramification indices ', f1, f2,'.')
# Check working by showing ramification of containing maximal order



"""
This is the function that computes log and b-discrepancies of
all exceptional curves in a good log resolution. Here the good log
resolution is one chosen so all negative discrepancy curves are
guaranteed to be found there. The data will be recorded in a dictionary.
"""
def resolve(n,z1,z2,g1,g2):
    
    no_blowups = max(g1,g2)*n -1
    # crude bound on number of blowups to ensure good log resolution
    print(f'We need to blowup set of nodes  at most {no_blowups} times.')
    no_curves = 2**no_blowups # number of exc curves actuall 2**no_blowups -1
    curves = {0:[1/(f1*g1)-1,z1] , no_curves:[1/(f2*g2)-1,z2]}
    # Create dictionary of curve data. Keys are 0 up to no_curves.
    # Current values are lists of [log discrepancy, ram of max order] 
    step = no_curves # Final indices for strict transforms of initial ram 
                     # curves are `step' apart

    
    for i in range(1,no_blowups+1):

        print(f'This is the {i}-th blowup. New curves have  log and b- discrepancy')
        step = int(step/2) # indices of i-th blowup curves have indices 'step' apart
        l = len(curves)-1 # the number of new exceptional curves
        is_good_res = True
        # We'll use this to test if the i-th blowup gives a good resolution

        for j in range(l):
            log_disc = 1 + curves[2*j*step][0] + curves[(2*j+2)*step][0]
            #  Give log discrepancy of j-th new exceptional in i-th blowup
            is_good_res =  is_good_res and (log_disc >=0)
            ram = (curves[2*j*step][1] + curves[(2*j+2)*step][1])%n
            # Gives ramification along j-th new exceptional in i-th blowup
            b_disc = log_disc + 1-1/order(ram,n)
            # Compute the b-discrepancy
            curves[(2*j +1)*step] = [log_disc, ram, b_disc]
            # Update curves dictionary with entry for j-th curve in i-th blowup
            print(log_disc,b_disc, end=';')

        if is_good_res:
            break
        else:
            print(' DONE', end='\n\n')

    print('\n\n Good log resolution has been achieved.')

resolve(n,z1,z2,g1,g2)
