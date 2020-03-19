# generateur d'instances de EVRPTW
# generateur de graphe grille
# le sommet source s est le sommet=1 (0 en python), le sommet puit t est le sommet=nbre_sommet (nbre_sommet-1 en python)
import random
m=7 # nombre de rangÃ©es
n=12 # nombre de colonnes
nbre_sommet=m*n+2
nbre_arc=(n-2)*(5*m-4)+5*m-2
nbre_arc_penalise=10 #nombre arcs penalises
test_penalite=0
val_penalite=10000

#fichier_nom="grille_3x5.txt"
#fichier_nom="grille_3x5_penalite.txt"
#chemin=r"C:\Users\math-\OneDrive\Documents\Courses\MPRO\RORT\RORT\data\grille_test.txt"

cost_max=10 # cout max d'un arc = poids = distance max
penalite_max=10 # penalite max = temps max de parcours d'un arc
# ce qui concerne la batterie
taux_consom=0.5 #0.5
taux_recharge=1.0 # vitesse de recharge
capacite=penalite_max # capacite de la batterie
ratio_station=4 # une station de recharge tous les 2 (ratio_recharge) sommets

Table_sommet=nbre_sommet*[0] # donnees sur sommets
for i in range(0,nbre_sommet):
    Table_sommet[i]=3*[0] # 3 donnees mise a zero

A=nbre_sommet*[0] # matrice d'adjacence
# on la met Ã  zero
for i in range(0,nbre_sommet):
    A[i]=nbre_sommet*[0]
C=nbre_sommet*[0] # matrice des couts
# on la met Ã  zero
for i in range(0,nbre_sommet):
    C[i]=nbre_sommet*[0]
D=nbre_sommet*[0] # matrice des penalites
# on la met Ã  zero
for i in range(0,nbre_sommet):
    D[i]=nbre_sommet*[0]

# init donnees sur le sommet 0 et le dernier
# fenetre temps debut , fenetre temps fin , type client ou station
Table_sommet[0][0]=0;Table_sommet[0][1]=penalite_max;Table_sommet[0][2]=0;
Table_sommet[nbre_sommet-1][0]=0;Table_sommet[nbre_sommet-1][1]=penalite_max*(n+1);Table_sommet[nbre_sommet-1][2]=0;
# init donnees sur les autres sommets
for lin in range(0,m):
    for col in range(0,n):
        num_sommet=lin + m*col + 1
        fenetre_deb=0 #col; # numero de colonne juste pour ne pas mettre zero
        fenetre_fin=(col+1)*penalite_max;
        #if lin==0 : type_sommet=1 # on met station recharge sur la ligne 0
        #else: type_sommet=0 # sommet client
        # on met station recharge sur les sommets pairs en comptant Ã  partir de zero
        if num_sommet%ratio_station==0 : type_sommet=1 # on met station recharge
        else: type_sommet=0 # sommet client
        Table_sommet[num_sommet][0]=fenetre_deb
        Table_sommet[num_sommet][1]=fenetre_fin
        Table_sommet[num_sommet][2]=type_sommet

# init matrice adjacence
# sommets partant de 1 (0 en python)
for i in range(0,m):
  A[0][1+i]=1
# colonne 1 (0 en python)
col=0
for lin in range(0,m):
  num_sommet=lin + col*m + 1 # numero du sommet
  num_sommet1=lin + (col+1)*m + 1 # numero du sommet dans la colonne suivante (sur la mÃªme ligne)
  A[num_sommet][num_sommet1]=1
  if lin>0 : A[num_sommet][num_sommet1-1]=1
  if lin<m-1 : A[num_sommet][num_sommet1+1]=1
# colonnes suivantes sauf la derniere n (n-1 en python)
for col in range(1,n-1):
  for lin in range(0,m):
    num_sommet=lin + col*m + 1 # numero du sommet
    num_sommet1=lin + (col+1)*m + 1 # numero du sommet dans la colonne suivante
    A[num_sommet][num_sommet1]=1
    if lin>0 : A[num_sommet][num_sommet1-1]=1
    if lin<m-1 : A[num_sommet][num_sommet1+1]=1
    if lin>0 : A[num_sommet][num_sommet-1]=1 # sommet ligne dessous
    if lin<m-1 : A[num_sommet][num_sommet+1]=1 # sommet ligne dessus
# derniere colonne n (n-1 en python)
col=n-1
for lin in range(0,m):
  num_sommet=lin + col*m + 1 # numero du sommet
  num_sommet1=nbre_sommet-1 # numero du sommet final t
  A[num_sommet][num_sommet1]=1


for i in range(0,nbre_sommet):
  print (A[i])


print("\n")
germe=1 # germe du random
random.seed(germe) # initialise random avec germe
# matrice cout
for i in range(0,nbre_sommet):
    for j in range(0,nbre_sommet):
        if  A[i][j]==1 : #i!=0 and j!=nbre_sommet-1 and
            cout = random.randint(1,cost_max)
            C[i][j]=cout

for i in range(0,nbre_sommet):
  print (C[i])

print("\n")
# matrice penalite
for i in range(0,nbre_sommet):
    for j in range(0,nbre_sommet):
        if  A[i][j]==1 :  # i!=0 and j!=nbre_sommet-1 and
            cout = random.randint(1,penalite_max)
            D[i][j]=cout
            if test_penalite==1 : # on est en mode penalite
              D[i][j]=val_penalite

for i in range(0,nbre_sommet):
  print (D[i])


# impression dans fichier
mon_fichier = open(r"C:\Users\math-\OneDrive\Documents\Courses\MPRO\RORT\RORT\data\grille_test.txt",'w')
chaine="%i "%nbre_sommet # convertir integer en chaine
mon_fichier.write(chaine)

chaine="%i "%nbre_arc # convertir integer en chaine
mon_fichier.write(chaine)

chaine="%f "%taux_consom # convertir integer en chaine
mon_fichier.write(chaine)

chaine="%f "%taux_recharge # convertir integer en chaine
mon_fichier.write(chaine)

chaine="%i \n"%capacite # convertir integer en chaine
mon_fichier.write(chaine)

# on imprime les sommets, leur fenetre de temps (debut fin)
# et leur type client=0 station recharge=1
for i in range(0,nbre_sommet):
    chaine="%i "%(i+1) # convertir integer en chaine
    mon_fichier.write(chaine)
    chaine="%i "%Table_sommet[i][0] # convertir integer en chaine
    mon_fichier.write(chaine)
    chaine="%i "%Table_sommet[i][1] # convertir integer en chaine
    mon_fichier.write(chaine)
    chaine="%i \n"%Table_sommet[i][2] # convertir integer en chaine
    mon_fichier.write(chaine)

#chaine="%i \n"%nbre_arc_penalise # convertir integer en chaine
#mon_fichier.write(chaine)
for i in range(0,nbre_sommet):
    for j in range(0,nbre_sommet):
        if A[i][j]==1 :
            chaine="%i "%(i+1) # on renumerote les sommets Ã  partir de 1
            mon_fichier.write(chaine)
            chaine="%i "%(j+1) # on renumerote les sommets Ã  partir de 1
            mon_fichier.write(chaine)
            chaine="%i "%C[i][j] # convertir integer en chaine
            mon_fichier.write(chaine)
            chaine="%i \n"%D[i][j] # convertir integer en chaine
            mon_fichier.write(chaine)


mon_fichier.close()
