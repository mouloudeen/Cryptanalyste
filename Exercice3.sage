load('Exercice1.sage')
# On crée les matrices de rétroaction avec cette fonction:
def matf(f):# matrice de transition 
    n = f.degree()
    A = Matrix(GF(2),n,n)
    for i in range(n-1):
        A[i,i+1] = 1
    tmp = f.list()[1:] # c_1,c_2,...,c_L
    tmp.reverse()   # c_L,c_{L-1}, ... , c_1
    A[n-1] = vector(tmp)
    return A

def A5_2_stepret(R1, R2, R3, R4, f1, f2, f3, f4):
    #nos matrices de retroactions:
    A1 = matf(f1)
    A2 = matf(f2)
    A3 = matf(f3)
    A4 = matf(f4)
    
    
    # On calcule la fonction majorité 
    # sur les indices 6, 13 et 9 de R4.
    m = Maj(R4[6], R4[13], R4[9])
    
    
    # On a 3 choix possibles.
    # Si à l'indice 6 de R4 est égale à m
    # on met à jour R1 en ignorant sa sortie.
    
    if R4[6] == m :
        x1 =vector(R1)
        x1 = A1*x1
        R1 = Sequence(x1)
       
        
        
    # Si à l'indice 13 de R4 est égale à m
    # on met à jour R2 en ignorant sa sortie.
    if  R4[13] == m :
        x2 = vector(R2)
        x2 = A2*x2
        R2 = Sequence(x2)
          
            
            
    # Si à l'indice 9 de R4 est égale à m 
    # on met à jour R3 en ignorant sa sortie.
    if  R4[9] == m :       
        x3 = vector(R3)
        x3 = A3*x3
        R3 = Sequence(x3)
        
        
   
        
        
    # on met à jour R4 en ignorant la sortie.
    x4 = vector(R4)
    x4 = A4*x4
    R4 = Sequence(x4)   
    
    # On calcule y1, y2 et y3 avec les formules données.
    y1 = R1[0] + Maj(R1[3], R1[4] + 1, R1[6])
    y2 = R2[0] + Maj(R2[8], R2[5] + 1, R2[12])
    y3 = R3[0] + Maj(R3[4], R3[9] + 1, R3[6])
    
    # On ressort le résultat de l'additions des 3 variables
    return (y1 + y2 + y3),(R1, R2, R3, R4)



# On suppose qu'on connait R4 et donc on va utiliser les K et IV données pour l'exemple du (1)

# La clef K
K = Sequence([GF(2)(0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

# IV
IV = Sequence([GF(2)(1), 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# On utilise les mêmes polynômes de rétroactions que l'algorithme "A5_2_init"
PR.<X> = PolynomialRing(GF(2))
f1 = X**19 + X**18 + X**17 + X**14 + 1
f2 = X**22 + X**21 + 1
f3 = X**23 + X**22 + X**21 + X**8 + 1
f4 = X**17 + X**12 + 1
    
# On crée le R4 que l'on connait.
_, _, _, R4 = A5_2_init(K, IV, f1, f2, f3, f4)  


# Pour déclarer les 64 inconnues que l'on utilisera dans les registres R1, R2 et R3.
BPR = BooleanPolynomialRing(64, 'x')
v = BPR.gens()

# On crée les R1, R2 et R3 avec les inconnues associés.
R1 = Sequence([v[i] for i in range(19)])
R2 = Sequence([v[i] for i in range(19,41)])
R3 = Sequence([v[i] for i in range(41,64)])

# On exécute 99 fois la fonction "A5_2_step" en ignorant 
# son bit de sortie 
for i in range(99):
     _, (R1 ,R2 ,R3 ,R4) = A5_2_stepret(R1, R2, R3, R4, f1, f2, f3, f4)
    
# On exécute 228 fois la fonction "A5_2_step" en utilisant en ignorant
# bit de sortie.
for i in range(228):
    _, (R1 ,R2 ,R3 ,R4) = A5_2_stepret(R1, R2, R3, R4, f1, f2, f3, f4)

# On crée une liste d'équations linéaires.

Eql = [R1[i] for i in range(19)] + [R2[i] for i in range(22)] + [R3[i] for i in range(23)]
print("Taille de la liste d'équation =", len(Eql))

# On voit les 64 équations.
print("Les 64 équations sont : ")
for i in range(64):
    print("xprime[%d] = %s" %(i, Eql[i]))