# Pour déclarer les 64 inconnues que l'on utilisera dans les registres R1, R2 et R3.
BPR = BooleanPolynomialRing(64, 'x')
v = BPR.gens()

# On sait x3 = x24 = x45 = 1, sont les variables connues, on crée M1 la liste des monômes
# de degré 1 sans ses 3 variables.
M1 = [v[i] for i in range(3)]+[v[i] for i in range(4,24)]+[v[i] for i in range(25,45)]+[v[i] for i in range(46,64)]
len(M1)

# On crée M la liste de tous les monômes possibles.
M = [M1[i] for i in range(61)]+[M1[i]*M1[j] for i in range(18) for j in range(i+1,18)]+[M1[i]*M1[j] for i in range(18,39) for j in range(i+1,39)]+[M1[i]*M1[j] for i in range(39,61) for j in range(i+1,61)]
print("M =", M)
print("La longueur de M est :", len(M))
# On voit que sa longueur est bien 655
N = 228
L = 655