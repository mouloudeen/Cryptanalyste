load('Exercice16.sage')

def changevar1e(f):
    
    # On recupère le degrée du polynômes pour créer les listes de monômes.
    d = f.degree()
    
    # On recherche dans la fonction tous les coefficient des mônomes de la forme x**i*y**j avec i et j entre 1 et d
    for i in range(1,d):
        for j in range(1,d):
        # On calcule tous les coefficients de chaque monôme de la forme x**i*y**j
            k = f.coefficient(x**i*y**j)
            if i < j:
                # Comme la puissance de x est supérieur à la puissance de y.
                f = f  - k*((x**i*y**j) - ((u-1)**i*y**(j-i)))
             
            if i ==j :
                #Comme la puissance de x est égale à la puissance de y.
                f = f   - k*((x**i*y**j) - (u-1)**i)
            
            if j<i :
                # Comme la puissance de y est supérieur à la puissance de x.
                f = f - k*((x**i*y**j) -  ((u-1)**j*x**(i-j)))
            
    return f

def lattice(m,t,N,e,gamma):
    e1 = ZZ(e)
    # D'après l'énoncé:
    X = integer_ceil(e1**gamma)
    Y = integer_ceil(e1**(0.5))
    U = 1+X*Y
    
    # On crée la fonction ft
    
    ft = u + (N+1)*x
    
    # On crée la liste de famille de polynomes ht
    H = [y**j*ft**k*e**(m-k) for j in range(1,t+1) for k in [integer_floor(m/t)*j..m]]
    
    # Pour chaque polynôme ht on fait le changement de x*y par (u-1)
    H1 = [changevar1e(i) for i in H]
   
    # On crée la liste des deux familles de pôlynomes.
    L =[x**i*ft**k*e**(m-k) for k in range(m+1) for i in range(m-k+1)]+H1 
    
    # On rajoute x*X,y*Y,u*U
    L1 = [i.subs({x:x*X, y:y*Y, u:u*U}) for i in L]
    
    # On crée la liste des monômes.
    M =[x**(k-i)*u**i for k in range(m+1) for i in range(k+1)]+[y**j*u**k for j in range(1,t+1) for k in [integer_floor(m/t)*j..m]]
    
    # On crée la matrice avec les coefficients de chaques monômes par rapport aux polynômes.
    Mat = matrix(ZZ,len(L1),len(M))
    for i in range(len(L1)):
        for j in range(len(M)):
            Mat[i,j] = L1[i].monomial_coefficient(M[j])
    return Mat

# On essaie avec un petit exemple
N,e,d =gen(5,0.7)
print("N = %d e =%d" %(N, e))
gamma = 0.7
X = integer_ceil(e**gamma)
Y = integer_ceil(e**(0.5))
U = 1+X*Y
print("X = %d  Y = %d  U = %d" %(X, Y, U))
print("\n")
print(" lattice(2, 1, %d, %d, %.2f) = " %(N, e, gamma))
print(lattice(2,1,N,e,gamma))
