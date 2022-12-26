load('Exercice5.sage')
load('Exercice6.sage')

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

# On crée les R1, R2 et R3 avec les inconnues associés et leurs bits connus.
R1 = Sequence([v[i] for i in range(3)]+[1]+[v[i] for i in range(4,19)])
R2 = Sequence([v[i] for i in range(19,24)]+[1]+[v[i] for i in range(25,41)])
R3 = Sequence([v[i] for i in range(41,45)]+[1]+[v[i] for i in range(46,64)])

# On calcule la suite chiffré z avec ses 64 inconnues.
zprimev2 = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)

# On a crée une liste d'équations quadratiques et on l'affiche les 228 équations.
print("Taille de la liste d'équation =", len(zprimev2))

# On voit les 228 équations.
print("Les %d équations sont : " %(N))
for i in range(228):
    print("zprime%d = %s" %(i, zprimev2[i]))

print("\n")
# On crée la matrice Mat de taille N*L dans F2.
Mat = matrix(GF(2),N,L)

# On crée le vecteur de longueur N.
Vec = vector(GF(2),N)

# On utilise le zprime calculé dans le (5).
# On utilise les indications de l'énoncé.
for i in range(N):
    Vec[i] = zprimev2[i].monomial_coefficient(1)
    for j in range(L):
        Mat[i,j] = zprimev2[i].monomial_coefficient(M[j])
        
print("Vec =", Vec)
print("\n")
print("Mat =", Mat)