print('starting')
# intial data
n=3
# f is in Z/n
f1=2
g1=1
f2=1
g2=1
# example intial data
curves=[[0,f1,g1],[0,g2,g2]] # our list of curves with [C^2,f,g]
# our list of curves as a list of lists
# we will insert new curves as we move to the log resolution
# we may expand the numbers recorded for each curve to include
# things like ramification index, discrepency, etc.
print(curves)
# first blowup
curves.insert(1,[-1,f1+f2,1])
print(curves)

i=1
# blowuping up at the intersection of curve[i] and curve[i+1]
leftcurve=curves[i]
rightcurve=curves[i+1]
curves[i][0] = curves[i][0]-1
curves[i+1][0] = curves[i+1][0]-1
curves.insert(i,[-1,curves[i][2]+curves[i+1][2],1])
print(curves)
