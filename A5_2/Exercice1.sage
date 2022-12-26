def LFSR_step(state,f):
    L = len(state)
    out = state[0]
    state = state[1:]+[sum(f[j+1]*state[L-1-j] for j in range(L))]
    return out,state


def A5_2_init(K, IV, f1, f2, f3, f4):
    
    # On initialise les registres à 0.
    R1 = Sequence([GF(2)(0) for _ in range(19)]) 
    R2 = Sequence([GF(2)(0) for _ in range(22)])  
    R3 = Sequence([GF(2)(0) for _ in range(23)])
    R4 = Sequence([GF(2)(0) for _ in range(17)])
    
    
    for i in range(64):
        # On met à jour chaque registre avec son polynome
        # de rétroaction respectif.
        _, R1 = LFSR_step(R1, f1)
        _, R2 = LFSR_step(R2, f2)
        _, R3 = LFSR_step(R3, f3)
        _, R4 = LFSR_step(R4, f4)
        
       
        
        # On Xor pour chaque derniers bits de chaque registre 
        # avec l'indice i (qui est le nombre tours faits) 
        # de la clef K.
        R1[18] =  R1[18] + K[i]
        R2[21] =  R2[21] + K[i]
        R3[22] =  R3[22] + K[i]
        R4[16] =  R4[16] + K[i]
        
        
        
        
    for i in range(22):
        
        # On met à jour chaque registre avec son polynome
        # de rétroaction respectif.
        _, R1 = LFSR_step(R1, f1)
        _, R2 = LFSR_step(R2, f2)
        _, R3 = LFSR_step(R3, f3)
        _, R4 = LFSR_step(R4, f4)
        
        
        
        # On Xor pour chaque derniers bits de chaque registre
        # avec l'indice i du vecteur initial IV. 
        R1[18] =  R1[18] + IV[i]
        R2[21] =  R2[21] + IV[i]
        R3[22] =  R3[22] + IV[i]
        R4[16] =  R4[16] + IV[i]
        
        
    
    # On fixe certains bits des registres.
    R1[3] = 1
    R2[5] = 1
    R3[4] = 1
    R4[6] = 1
    
    
    return R1, R2, R3, R4


def Maj(a, b, c):
    # 3 cas possibles pour la fonction majorité.
    m = 0 
    if a == b:
        m = a
        
    if a == c:
        m = a
        
    if b == c :
        m = b
        
    return m


def A5_2_step(R1, R2, R3, R4, f1, f2, f3, f4):
    
    # On calcule la fonction majorité 
    # sur les indices 6, 13 et 9 de R4.
    m = Maj(R4[6], R4[13], R4[9])
    
    
    # On a 3 choix possibles.
    # Si à l'indice 6 de R4 est égale à m
    # on met à jour R1 en ignorant sa sortie.
    
    if R4[6] == m :
        _,  R1 = LFSR_step(R1, f1)
       
        
        
    # Si à l'indice 13 de R4 est égale à m
    # on met à jour R2 en ignorant sa sortie.
    if  R4[13] == m :
        _,  R2 = LFSR_step(R2, f2)
          
            
            
    # Si à l'indice 9 de R4 est égale à m 
    # on met à jour R3 en ignorant sa sortie.
    if  R4[9] == m :       
        _,  R3 = LFSR_step(R3, f3)
        
        
   
        
        
    # on met à jour R4 en ignorant la sortie.
    _, R4 = LFSR_step(R4, f4)   
    
    # On calcule y1, y2 et y3 avec les formules données.
    y1 = R1[0] + Maj(R1[3], R1[4] + 1, R1[6])
    y2 = R2[0] + Maj(R2[8], R2[5] + 1, R2[12])
    y3 = R3[0] + Maj(R3[4], R3[9] + 1, R3[6])
    
    # On ressort le résultat de l'additions des 3 variables
    return (y1 + y2 + y3),(R1, R2, R3, R4)


def A5_2_production(R1, R2, R3, R4, f1, f2, f3, f4):
    
    # On exécute 99 fois la fonction "A5_2_step" en ignorant 
    # son bit de sortie 
    for i in range(99):
        _, (R1 ,R2 ,R3 ,R4) = A5_2_step(R1, R2, R3, R4, f1, f2, f3, f4)
    
    # On prend z une suite iniatilisé à une liste vide.
    z = []
    
    # On exécute 228 fois la fonction "A5_2_step" en utilisant son
    # bit de sortie.
    for i in range(228):
        zi, (R1 ,R2 ,R3 ,R4) = A5_2_step(R1, R2, R3, R4, f1, f2, f3, f4)
        z.append(zi)
        
    
    return z


def A5_2(K, IV):
    
    # On crée les polynomes de rétroaction
    PR.<X> = PolynomialRing(GF(2))
    f1 = X**19 + X**18 + X**17 + X**14 + 1
    f2 = X**22 + X**21 + 1
    f3 = X**23 + X**22 + X**21 + X**8 + 1
    f4 = X**17 + X**12 + 1
    
    # Phase d'initialisation.
    R1, R2, R3, R4 = A5_2_init(K, IV, f1, f2, f3, f4)
    
    
    # Phase de productions.
    return A5_2_production(R1, R2, R3, R4, f1, f2, f3, f4)


load('A5_2-test-vector.sage')

# On va tester si notre fonction nous donne les bonnes résultats
# prévues .
zprime = A5_2(K, IV)

print("z = A5_2(K, IV) est ", zprime == z)

