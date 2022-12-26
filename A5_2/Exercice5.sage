load('Exercice1.sage')

def Maj_al_n(a,b,c): # fonction majorité en forme algébrique normale.
    return a*b +b*c +a*c


def A5_2_stepv_2(R1, R2, R3, R4, f1, f2, f3, f4): 
    # On utilise la fonction majorité en forme algébrique normale.
    
    # On calcule la fonction majorité 
    # sur les indices 6, 13 et 9 de R4.
    m = Maj_al_n(R4[6], R4[13], R4[9])
    
    
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
    y1 = R1[0] +  Maj_al_n(R1[3], R1[4] + 1, R1[6])
    y2 = R2[0] +  Maj_al_n(R2[8], R2[5] + 1, R2[12])
    y3 = R3[0] +  Maj_al_n(R3[4], R3[9] + 1, R3[6])
    
    # On ressort le résultat de l'additions des 3 variables
    return (y1 + y2 + y3),(R1, R2, R3, R4)


def A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4): # On utilise la fonction majorité en forme algébrique normale
    # Dans cette version on rajoute N pour pourvoir utiliser avec n'importe quel N.
    # On exécute 99 fois la fonction "A5_2_step" en ignorant 
    # son bit de sortie 
    for i in range(99):
        _, (R1 ,R2 ,R3 ,R4) = A5_2_stepv_2(R1, R2, R3, R4, f1, f2, f3, f4)
    
    # On prend z une suite iniatilisé à une liste vide.
    z = []
    
    # On exécute N fois la fonction "A5_2_step" en utilisant son
    # bit de sortie.
    for i in range(N):
        zi, (R1 ,R2 ,R3 ,R4) = A5_2_stepv_2(R1, R2, R3, R4, f1, f2, f3, f4)
        z.append(zi)
        
    
    return z


def A5_2v_2(N,K, IV):
    
    # On crée les polynomes de rétroaction
    PR.<X> = PolynomialRing(GF(2))
    f1 = X**19 + X**18 + X**17 + X**14 + 1
    f2 = X**22 + X**21 + 1
    f3 = X**23 + X**22 + X**21 + X**8 + 1
    f4 = X**17 + X**12 + 1
    
    # Phase d'initialisation.
    R1, R2, R3, R4 = A5_2_init(K, IV, f1, f2, f3, f4)
    
    
    # Phase de productions.
    return A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)

# On suppose qu'on connait R4 et donc on va utiliser les K et IV données pour l'exemple du (1)

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

# On veut N =228
N = 228

# On calcule la suite chiffré z avec ses 64 inconnues.
zprime = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)

# On a crée une liste d'équations quadratiques et on l'affiche les 228 équations.
print("Taille de la liste d'équation =", len(zprime))

# On voit les 228 équations.
print("Les %d équations sont : " %(N))
for i in range(228):
    print("zprime[%d] = %s" %(i, zprime[i]))


