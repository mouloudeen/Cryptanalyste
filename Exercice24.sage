load('Exercice21.sage')
load('Exercice22.sage')

def factorisation(N,e,gamma,m):
    e1 = ZZ(e)
   # D'après l'énoncé:
    X = integer_ceil(e1**gamma)
    Y = integer_ceil(e1**(0.5))
    U = 1+X*Y
    t = integer_floor(m*(1-2*gamma))
    Mat =lattice(m,t,N,e,gamma)
    Matl = Mat.LLL()
    
    Pr.<x,y,u> = PolynomialRing(QQ)
    # On crée une liste de monômes
    M =[(x/X)**(k-i)*(u/U)**i for k in range(m+1) for i in range(k+1)]+[(y/Y)**j*(u/U)**k for j in range(1,t+1) for k in [integer_floor(m/t)*j..m]]
    
    # On crée les polynomes de chaque ligne de la matrice.
    P=[sum(Matl[i][j]*M[j] for j in range(len(M))) for i in range(len(M))]
    
    
    
    # On remplace u par 1+x*y
    Pchange = [changevar2(i) for i in P]
    
    
   
    # on initialise i=0
    i = 0
    P1 = Pchange[i]
    P2 = Pchange[i+1]
    
    # On met Pxy un nouvelle anneau polynomial en fonction de x et y
    Pxy.<x, y> = PolynomialRing(QQ)
    # On caste P1 et P2 en Pxy
    
    PP1 = Pxy(P1)
    PP2 = Pxy(P2)
   
    
    res = PP1.resultant(PP2,y)
    
     # On met Px un nouvelle anneau polynomial en fonction de x
    Px.<x> = PolynomialRing(QQ)
    
    resx = Px(res)
    
    
    while res.degree() < 0  and i!=len(M)-2:
        i = i+1
        P1 = Pchange[i]
        P2 = Pchange[i+1] 
        PP1 = Pxy(P1)
        PP2 = Pxy(P2)
        
    
        # On met Px un nouvelle anneau polynomial en fonction de x
   
        res = PP1.resultant(PP2,y)
        
        Px.<x> = PolynomialRing(QQ)
    
        resx = Px(res)
        
        
    
        
    
    # On regarde sur on est sorti de la liste de fonctions
    if i == len(M)-2:
        return (-1)
    # On cherche la premiere valeur   
    x0 = resx.roots()[0][0]
    
    
    #On reinjecte cette valeur sur le premier polynôme
    PPP1 =P2.subs(x=x0)
    
    # On met Py un nouvelle anneau polynomial en fonction de y
    Py.<y> = PolynomialRing(QQ)
    
    # On le caste en Py
    equy = Py(PPP1)
    
    # On cherche la racine de cette équation en fonction de y
    y0 =equy.roots()[0][0]
    
    # D'aprés les résultats de la question 16 et 17
    d = (x0*(N+1+y0)+1)/e1
    
    return (d )


#premier exemple demandé
N,e,d =gen(1024,0.2)
gamma = 0.2
m = 2

# On teste si on a trouvé le bon d
d1 = factorisation(N, e, gamma, m)

print("d = factorisation(%d, %d,%.2f, %d) est %s" %(N, e, gamma, m, d==d1))

# Pour la deuxième exemple, je n'ai pas réussit à debugger mon programme.