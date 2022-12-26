load('Exercice5.sage')
load('A5_2-3frames.sage')

IV= Sequence([GF(2)(0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) # pour z0
IVprime = Sequence([GF(2)(0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]) # pour z1
IVtilde = Sequence([GF(2)(0), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]) # pour z2

# On pose les registres des LFSRi Ri pour z0, Riprime pour z1 et Ritilde pour z2 avec 1 <= i <= 4
# Pour déclarer les 64 inconnues que l'on utilisera dans les registres R1, R2, R3, R1prime, R2prime, R3prime
# R1tilde, R2tilde et R3tilde
BPR = BooleanPolynomialRing(64, 'x')
v = BPR.gens()

# On crée les R1, R2 et R3 avec les inconnues associés et leurs bits connus.
R1 = Sequence([v[i] for i in range(3)]+[1]+[v[i] for i in range(4,19)])
R2 = Sequence([v[i] for i in range(19,24)]+[1]+[v[i] for i in range(25,41)])
R3 = Sequence([v[i] for i in range(41,45)]+[1]+[v[i] for i in range(46,64)])

print("R1 =",R1)
print("R2 =",R2)
print("R3 =",R3)
print("R4 =",R4)
print("\n")

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

# Pour faire les calculs je transforme les registres en vecteur.
X1 = vector(R1)
X2 = vector(R2)
X3 = vector(R3)
X4 = vector(R4)

# On calcul les registres des autres suites chiffrantes.
# On commence par celle de z1:
X1prime = X1 + vector([0 for _ in range(18)]+[IVprime[21]])
X2prime = X2 + vector([0 for _ in range(21)]+[IVprime[21]])
X3prime = X3 + vector([0 for _ in range(22)]+[IVprime[21]])
X4prime = X4 + vector([0 for _ in range(16)]+[IVprime[21]])

# on met un seul bit à 1 dans chaque registre :
X1prime[3] = 1
X2prime[5] = 1
X3prime[4] = 1
X4prime[6] = 1 

# On le remet en sequence
R1prime = Sequence(X1prime)
R2prime = Sequence(X2prime)
R3prime = Sequence(X3prime)
R4prime = Sequence(X4prime)

print("R1prime =",R1prime)
print("R2prime =",R2prime)
print("R3prime =",R3prime)
print("R4prime =",R4prime)
print("\n")

# et celle de z2:
X1tilde = X1 + A1 * vector([0 for _ in range(18)]+[IVtilde[20]])
X2tilde = X2 + A2 * vector([0 for _ in range(21)]+[IVtilde[20]])
X3tilde = X3 + A3 * vector([0 for _ in range(22)]+[IVtilde[20]])
X4tilde = X4 + A4 * vector([0 for _ in range(16)]+[IVtilde[20]])

# on met un seul bit à 1 dans chaque registre :
X1tilde[3] = 1
X2tilde[5] = 1
X3tilde[4] = 1
X4tilde[6] = 1 

# On le remet en sequence
R1tilde = Sequence(X1tilde)
R2tilde = Sequence(X2tilde)
R3tilde = Sequence(X3tilde)
R4tilde = Sequence(X4tilde)

print("R1tilde =",R1tilde)
print("R2tilde =",R2tilde)
print("R3tilde =",R3tilde)
print("R4tilde =",R4tilde)
print("\n")

# On pose N = 228
N = 228
# On calcule les z théoriques 
z = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)
zprime = A5_2_productionv_2(N, R1prime, R2prime, R3prime, R4prime, f1, f2, f3, f4)
ztilde = A5_2_productionv_2(N, R1tilde, R2tilde, R3tilde, R4tilde, f1, f2, f3, f4)


# On crée une liste d'équations
Eq =z+zprime+ztilde

# On crée une liste de résultats voulues
Res = z0+z1+z2

# On sait x3 = x24 = x45 = 1, sont les variables connues, on crée M1 la liste des monômes
# de degré 1 sans ses 3 variables.
M1 = [v[i] for i in range(3)]+[v[i] for i in range(4,24)]+[v[i] for i in range(25,45)]+[v[i] for i in range(46,64)]

# On crée M la liste de tous les monômes possibles.
M = [M1[i] for i in range(61)]+[M1[i]*M1[j] for i in range(18) for j in range(i+1,18)]+[M1[i]*M1[j] for i in range(18,39) for j in range(i+1,39)]+[M1[i]*M1[j] for i in range(39,61) for j in range(i+1,61)]

L = len(M)

# On crée la matrice Mat de taille N*L dans F2.
Mat = matrix(GF(2),3*N,L)

# On crée le vecteur de longueur N.
Vec = vector(GF(2),3*N)


# On utilise les indications de l'énoncé.
for i in range(3*N):
    Vec[i] = Eq[i].monomial_coefficient(1)
    for j in range(L):
        Mat[i,j] = Eq[i].monomial_coefficient(M[j])


# On rajoute le Vecteur aux résultats prévues
Right = [Vec[i]+Res[i] for i in range(3*N)]

Sol = Mat.solve_right(Right)

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

# Maintenant on calcul les autres Registres qui produisent les deux autres suites chiffrantes.
# On transforme les 4 registres trouvés en vecteur pour pouvoir calculer les autres registres.
X1 = vector(R1)
X2 = vector(R2)
X3 = vector(R3)
X4 = vector(R4)

# on calcul les registres qui produisent z1
X1prime = X1 + vector([0 for _ in range(18)]+[IVprime[21]])
X2prime = X2 + vector([0 for _ in range(21)]+[IVprime[21]])
X3prime = X3 + vector([0 for _ in range(22)]+[IVprime[21]])
X4prime = X4 + vector([0 for _ in range(16)]+[IVprime[21]])

# on met un seul bit à 1 dans chaque registre :
X1prime[3] = 1
X2prime[5] = 1
X3prime[4] = 1
X4prime[6] = 1 

# On le remet en sequence
R1prime = Sequence(X1prime)
R2prime = Sequence(X2prime)
R3prime = Sequence(X3prime)
R4prime = Sequence(X4prime)

print("R1prime =",R1prime)
print("R2prime =",R2prime)
print("R3prime =",R3prime)
print("R4prime =",R4prime)
print("\n")


# et celles qui produisentz2:
X1tilde = X1 + A1 * vector([0 for _ in range(18)]+[IVtilde[20]])
X2tilde = X2 + A2 * vector([0 for _ in range(21)]+[IVtilde[20]])
X3tilde = X3 + A3 * vector([0 for _ in range(22)]+[IVtilde[20]])
X4tilde = X4 + A4 * vector([0 for _ in range(16)]+[IVtilde[20]])

# on met un seul bit à 1 dans chaque registre :
X1tilde[3] = 1
X2tilde[5] = 1
X3tilde[4] = 1
X4tilde[6] = 1 

# On le remet en sequence
R1tilde = Sequence(X1tilde)
R2tilde = Sequence(X2tilde)
R3tilde = Sequence(X3tilde)
R4tilde = Sequence(X4tilde)

print("R1tilde =",R1tilde)
print("R2tilde =",R2tilde)
print("R3tilde =",R3tilde)
print("R4tilde =",R4tilde)
print("\n")

# on produit ztest et on teste si c'est égale au résultat voulu
ztest = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4)
print("z0 = A5_2_productionv_2(N, R1, R2, R3, R4, f1, f2, f3, f4) est ", ztest ==z0)

# On teste les deux suites chiffrantes avec leurs valeurs voulues
zprimetest = A5_2_productionv_2(N, R1prime, R2prime, R3prime, R4prime, f1, f2, f3, f4)
print("z1 A5_2_productionv_2(N, R1prime, R2prime, R3prime, R4prime, f1, f2, f3, f4) est ",  zprimetest == z1)

ztildetest = A5_2_productionv_2(N, R1tilde, R2tilde, R3tilde, R4tilde, f1, f2, f3, f4)
print("z2 = A5_2_productionv_2(N, R1tilde, R2tilde, R3tilde, R4tilde, f1, f2, f3, f4) est ",ztildetest == z2)
print("\n")

#On voit bien qu on a les bons registres.
# Maintenant on cherche la clef K 
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

# On va calculer l'idéal:
I  = ideal([L[i] + R[i] for i in range(81)])

print("L'idéal I :", I)
print("\n")

# On utilise la base groebner.
Iq = I.groebner_basis()
print("List(I.groebner_basis()) :" ,list(Iq))
print("\n")

# On cherche les coefficients qui additionnent les variables pour avoir une clef Kprime.
Kprime = [el.constant_coefficient() for el in Iq]
print("Kprime =", Kprime)
print("\n")

# On produit avec la clef et leurs IV respectifs
ztest1 = A5_2v_2(N,Kprime, IV)
zprimetest1 = A5_2v_2(N,Kprime, IVprime)
ztildetest = A5_2v_2(N,Kprime, IVtilde)

# on teste avec leurs valeurs voulus
print("z0 = A5_2v_2(N,Kprime, IV) est ", ztest1 == z0)
print("z1 = A5_2v_2(N,Kprime, IVprime) est ", zprimetest1 == z1)
print("z2 = A5_2v_2(N,Kprime, IVtilde) est ",  ztildetest == z2)

