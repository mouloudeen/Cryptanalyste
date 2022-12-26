load('Exercice5.sage')
load('A5_2-700.sage')
# On recupère les Ri de l'exercice 8:
R1 = [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0]
R2 = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
R3 = [1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0]
R4 = [1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1]

# Son IV est nulle
IV = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


def matf(f):# matrice de transition 
    n = f.degree()
    A = Matrix(GF(2),n,n)
    for i in range(n-1):
        A[i,i+1] = 1
    tmp = f.list()[1:] # c_1,c_2,...,c_L
    tmp.reverse()   # c_L,c_{L-1}, ... , c_1
    A[n-1] = vector(tmp)
    return A


# On utilise les mêmes polynômes de rétroactions que l'algorithme "A5_2_init"
PR.<X> = PolynomialRing(GF(2))
f1 = X**19 + X**18 + X**17 + X**14 + 1
f2 = X**22 + X**21 + 1
f3 = X**23 + X**22 + X**21 + X**8 + 1
f4 = X**17 + X**12 + 1

# On construit les matrices de retroactions pour chaque fonction de rétroaction
A1 = matf(f1)
A2 = matf(f2)
A3 = matf(f3)
A4 = matf(f4)


# Pour déclarer les 64 inconnues que l'on utilisera dans la clef K.
BPR = BooleanPolynomialRing(64, 'K')
v = BPR.gens()

# Comme dans l'algorithme de A5/2-init les registres sont vides des 4 LFSR sont nulles.
X1 = vector([0 for _ in range(19)]) # registre du LFSR1
X2 = vector([0 for _ in range(22)]) # registre du LFSR2
X3 = vector([0 for _ in range(23)]) # registre du LFSR3
X4 = vector([0 for _ in range(17)]) # registre du LFSR4

# On calcule les registres comme dans le résultat du (9).
# Dans la première boucle for:
for i in range(64):
    X1 = A1*X1+vector([0 for _ in range(18)]+[v[i]])
    X2 = A2*X2+vector([0 for _ in range(21)]+[v[i]])
    X3 = A3*X3+vector([0 for _ in range(22)]+[v[i]])
    X4 = A4*X4+vector([0 for _ in range(16)]+[v[i]])
    
# Comme IV est nulle, dans la deuxième boucle for on multiplie pour chaque registre par leurs matrices 
# de rétroaction qui sont mis à la puissance 22.
X1 = A1**22*X1
X2 = A2**22*X2
X3 = A3**22*X3
X4 = A4**22*X4

# A la fin de cet algorithme, on met un seul bit à 1 dans chaque registre :
X1[3] = 1
X2[5] = 1
X3[4] = 1
X4[6] = 1

# On fait une liste équation
L = list(X1)+list(X2)+list(X3)+list(X4)

# On concatène les registres trouvés dans le (8)
R = list(R1)+list(R2)+list(R3)+list(R4)

# On va calculer l'idéal :
I  = ideal([L[i] + R[i] for i in range(81)])
print("I lidéal :",I)
print("\n")
# On utilise la base groebner.
Iq = I.groebner_basis();
print("list(II.groebner_basis() :", list(Iq))
print("\n")

# On cherche les coefficients qui additionnent les variables pour avoir une clef Kprime.
Kprime = [el.constant_coefficient() for el in Iq]
print("Kprime =", Kprime)
print("\n")

#On calcule avec Kprime pour chercher les Riprime avec 1<=i<=4
R1prime, R2prime, R3prime,R4prime =A5_2_init(Kprime, IV, f1, f2, f3, f4)

# On teste si ça nous donne les bons résultats:
print("R1prime = R1 est ", R1prime == R1)
print("R2prime = R2 est ", R2prime == R2)
print("R3prime = R3 est ", R3prime == R3)
print("R4prime = R4 est ", R4prime == R4)
print("On a trouvé le bon K")
