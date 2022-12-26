load('Exercice5.sage')
load('A5_2-700.sage')

N = 700
# On considère une exécution de A5/2 donnant N = 700 bits de suite chiffrante avec
# On a la valeur de R4 après la phase d'initialisation
# et z la suite chiffrante de 700 bits

# Pour déclarer les 64 inconnues que l'on utilisera dans les registres R1, R2 et R3.
BPR = BooleanPolynomialRing(64, 'x')
v = BPR.gens()

# On utilise les mêmes polynômes de rétroactions que l'algorithme "A5_2_init"
PR.<X> = PolynomialRing(GF(2))
f1 = X**19 + X**18 + X**17 + X**14 + 1
f2 = X**22 + X**21 + 1
f3 = X**23 + X**22 + X**21 + X**8 + 1
f4 = X**17 + X**12 + 1

# On crée les R1, R2 et R3 avec les inconnues associés et leurs bits connus.
R1 = Sequence([v[i] for i in range(3)]+[1]+[v[i] for i in range(4,19)])
R2 = Sequence([v[i] for i in range(19,24)]+[1]+[v[i] for i in range(25,41)])
R3 = Sequence([v[i] for i in range(41,45)]+[1]+[v[i] for i in range(46,64)])

# On calcule la suite chiffré z avec ses 64 inconnues.
zprime = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)

# On a crée une liste d'équations quadratiques et on l'affiche les 700 équations.
print("Taille de la liste d'équation =", len(zprime))

# On voit les 700 équations.
print("Les %d équations sont : " %(N))
for i in range(N):
    print("zprime[%d] = %s" %(i, zprime[i]))

print("\n")

# On crée la matrice Mat de taille N*L dans F2.
Mat = matrix(GF(2),N,L)

# On crée le vecteur de longueur N.
Vec = vector(GF(2),N)

# On utilise le zprime .
# On utilise les indications de l'énoncé.
for i in range(N):
    Vec[i] = zprime[i].monomial_coefficient(1)
    for j in range(L):
        Mat[i,j] = zprime[i].monomial_coefficient(M[j])
        
print("Vec =", Vec)
print("\n")
print("Mat =", Mat)
print("\n")

# Pour résoudre l'équation linéaire, on va utiliser la fonction solve_right() donnée dans l'énoncé.
# On additionne vec à z.
Right = [Vec[i] + z[i] for i in range(N)]
#on résoud l'équation grâce à la fonction, en mettant tous les inconnues à gauche de l'égalité
# et le reste à droite, on pose sol le résultat de la fonction:
Sol = Mat.solve_right(Right)
print("Sol =",Sol)
print("\n")

#On a besoin que les 61 premiers bits du sol pour récupérer 
# R1, R2 et R3 :
R1 = Sequence([Sol[i] for i in range(3)]+[1]+[Sol[i] for i in range(3,18)])
R2 = Sequence([Sol[i] for i in range(18,23)]+[1]+[Sol[i] for i in range(23,39)])
R3 = Sequence([Sol[i] for i in range(39,43)]+[1]+[Sol[i] for i in range(43,61)])
print("R1 =", R1)
print("R2 =", R2)
print("R3 =", R3)
print("R4 =", R4)
print("\n")

# On teste si ça nous donne le bon résultat:
ztest = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)
print(" z = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4) est", ztest == z)

