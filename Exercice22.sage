Pr.<x,y,u> = PolynomialRing(ZZ)

def changevar2(f):
    #On change u en 1+x*y 
    f = f.substitute(u = 1+x*y)
    return f

# On teste une petite fonction 
f = x*u +y*u**2 +2
print("f =",f)
print("changevar2(f) =", changevar2(f)) 