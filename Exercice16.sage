# On crée le Polynomial Ring de Z
Pr.<x,y,u> = PolynomialRing(ZZ)

def gen(k,gamma):
    # On prend au hasard p le premier nombre premier.
    p = random_prime(2**k,lbound = 2**(k-1), proof=false)
    
    # On prend au hasard q le deuxième nombre premier.  
    q = random_prime(2**k,lbound = 2**(k-1), proof=false)
    
    # Si P et Q sont égal, on recherche pour une autre valeur qui est
    # premier.
    if q==p:
        q=next_prime(q)
        
    N = p*q
    
    # On calcule phi(N)= (p-1)*(q-1)
    phi =(p-1)*(q-1)
    
    # On calcule d'abord d
    d = floor(N**gamma)
    
    # On cherche le plus grand d possible inferieur à N^gamma
    # et qui soit premier avec phi
    while gcd(phi,d) != 1:
        d = d-1
        
    # On calcule e=d^(-1)mod(phi)
    e = d.inverse_mod(phi)
    return N,e,d

print("gen(2048,0.5) =", gen(2048,0.5))  
