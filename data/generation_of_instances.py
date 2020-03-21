from random import randrange
import math

#Nombre de clients
n = 10
#Nombre d'arcs
#nb_arcs = 1
tps_max = 30

r = 0.5
q = 1
Q = 10

taille_grille = 10 # grille carrée taille_grille x taille_grille
nb_divisions = 2 # diviseur de taille_grille : nb_divisions**2 zones au total
ratio = int(taille_grille/nb_divisions)
n_tot = n + nb_divisions*nb_divisions

clients = []
for i in range (n) :
    clients.append([randrange(taille_grille), randrange(taille_grille)])

print(clients)

stations = []
for i in range (nb_divisions):
    for j in range (nb_divisions):
        stations.append([ratio*i + randrange(ratio), ratio*j + randrange(ratio)])

print(stations)

sommets = clients + stations

temps_stations = []
temps_clients = []
for i in range (nb_divisions**2):
    temps_stations.append([0, tps_max])
for i in range (n):
    temps_clients.append([0, tps_max])

arcs = [[0]*n_tot]*n_tot
for i in range(n_tot):
    for j in range(n_tot):
        if i != j :
            arcs[i][j] = 1

poids = [[0]*n_tot]*n_tot
for i in range(n_tot):
    for j in range(n_tot):
        poids[i][j] = 1

distances = [[0]*n_tot]*n_tot
for i in range(n_tot):
    for j in range(n_tot):
        if arcs[i][j] == 1 :
            distances[i][j] = math.sqrt((sommets[i][0] - sommets[j][0])**2 + (sommets[i][1] - sommets[j][1])**2)




mon_fichier = open(r"C:\Users\math-\OneDrive\Documents\Courses\MPRO\RORT\RORT\data\instance1.txt",'w')
chaine="%i "%n_tot
mon_fichier.write(chaine)
chaine="%i "%(n_tot**2 - n_tot)
mon_fichier.write(chaine)
chaine="%f "%r
mon_fichier.write(chaine)
chaine="%f "%q
mon_fichier.write(chaine)
chaine="%i \n"%Q
mon_fichier.write(chaine)

for i in range(n):
    chaine="%i "%(i+1)
    mon_fichier.write(chaine)
    chaine="%i "%temps_clients[i][0]
    mon_fichier.write(chaine)
    chaine="%i "%temps_clients[i][1]
    mon_fichier.write(chaine)
    chaine="%i \n"%0
    mon_fichier.write(chaine)

for i in range(nb_divisions*nb_divisions):
    chaine="%i "%(n+i+1)
    mon_fichier.write(chaine)
    chaine="%i "%temps_stations[i][0]
    mon_fichier.write(chaine)
    chaine="%i "%temps_stations[i][1]
    mon_fichier.write(chaine)
    chaine="%i \n"%1
    mon_fichier.write(chaine)

for i in range(n_tot):
    for j in range(n_tot):
        if arcs[i][j]==1 :
            chaine="%i "%(i+1)
            mon_fichier.write(chaine)
            chaine="%i "%(j+1)
            mon_fichier.write(chaine)
            chaine="%i "%poids[i][j]
            mon_fichier.write(chaine)
            chaine="%i \n"%distances[i][j]
            mon_fichier.write(chaine)


mon_fichier.close()
