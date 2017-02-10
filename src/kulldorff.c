#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>


double rapportVraisemblanceKulldorf(double nz, double muz, double N)
{
    double tmp1, tmp2, res;
    if (nz <= muz) {
	res=1.0;
	return(res);
    }
    tmp1 = R_pow(nz/muz, nz);
    tmp2 = R_pow((N-nz)/(N-muz), N-nz);
    res = tmp1*tmp2;
    return(res);
}

SEXP trouveCluster(SEXP voisinage, SEXP donnees, SEXP theorique, SEXP NparT)
{
    /* Déclaration des variables utilisées */
    int Nt, Npts, Ndiam, i, t, d, j1t, j2t, v, j,k, tma, dma;
    SEXP Dc, t1c, lambda, vec, j1, j2, nz, muz, lam, resu;
    double N, vma;

    /* Quelques variables utiles */
    Nt = length(donnees); /* Nombre de semaines (nb colonnes de donnees) */
    Npts = length(VECTOR_ELT(donnees,0)); /* Nombre de pixels dans la carte (nb 
					     lignes dans donnees) */
    Ndiam = length(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(voisinage,0),1),0)); /* Nombre de diamètres 
									    testes récupéré à
									    partir de voisinage */

    /* Initialisation des indices */
    j1t=0;
    j2t=0;

    /* Allocation de mémoire */
    PROTECT(Dc = allocVector(INTSXP, Npts)); /* Va contenir, pour chaque pixel, l'indice 
						du diamètre correspondant au cluster optimum */
    PROTECT(t1c = allocVector(INTSXP, Npts)); /* Va contenir, pour chaque pixel, l'indice
						 de la semaine correspondant au cluster optimum */
    PROTECT(lambda = allocVector(REALSXP, Npts)); /* Va contenir, pour chaque pixel, la valeur
						     du rapport de vraisemblance correspondant au
						     cluster optimum */
    PROTECT(nz = allocVector(REALSXP, Nt*Ndiam)); /* Va être renouvelé à chaque pixel: contiendra
						     le nombre observé de points pour un 
						     diamètre (lignes) et un nombre de 
						     semaines (colonnes) donnés
						   */
    PROTECT(muz = allocVector(REALSXP, Nt*Ndiam)); /* Va être renouvelé à chaque pixel: contiendra
						     le nombre théorique de points pour un 
						     diamètre (lignes) et un nombre de 
						     semaines (colonnes) donnés
						   */
    PROTECT(lam = allocVector(REALSXP, Nt*Ndiam)); /* Va être renouvelé à chaque pixel: contiendra
						      les rapports de vraisemblance pour un 
						     diamètre (lignes) et un nombre de 
						     semaines (colonnes) donnés
						   */
    PROTECT(resu = allocVector(VECSXP, 3)); /* Liste utilisée en sortie: première colonne:
					       le rapport de vraisemblance correspondant
					       à chaque pixel (lambda), deuxième colonne:
					       l'indice (R, donc commençant à 1) du diamètre
					       correspondant au cluster optimum pour ce pixel,
					       et troisième colonne l'indice R du nombre de 
					       semaines correspondant au cluster optimum 
					       pour ce pixel */
    
    
    /* Pour chaque pixel, calcul des clusters optimums */
    for (i = 0; i < Npts; i++) {

	/* On récupère le voisinage du point: vec contient l'indice des voisins iv
	   j1 est un vecteur de longueur Ndiam contenant l'indice, dans vec, du début
	   du subset des iv correspondant à ce diamètre, mais pas aux diamètres plus petits. 
	   j2 est un vecteur de longueur Ndiam contenant l'indice dans vec, de la fin de ce
	   subset.
	 */
	vec=VECTOR_ELT(VECTOR_ELT(voisinage, i), 0);
	j1=VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(voisinage, i), 1), 0);
	j2=VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(voisinage, i), 1), 1);
	
	/* Initialisation de la matrice des effectifs nz et muz, ainsi que
	   des rapports de vraisemblance à 0. En lignes: les diamètres, et en colonnes:
	   les semaines. La matrice fait donc Ndiam lignes et Nt colonnes
	*/
	for (j=0; j<Ndiam; j++) {
	    for (k=0; k<Nt; k++) {
		REAL(nz)[j+(k*Ndiam)] = 0.0;
		REAL(muz)[j+(k*Ndiam)] = 0.0;
/*		REAL(lam)[j+(k*Ndiam)] = 0.0;   */
	    }
	}
	
	
	/* Pour chaque semaine, */
	for (t = 0; t<Nt; t++) {
	    
	    /* Nombre total de cas observés sur la période */
	    N = REAL(NparT)[t];
	    
	    for (d = 0; d<Ndiam; d++) {
		
		/* Extraction des limites du subset des iv pour le diamètre d */
		j1t=INTEGER(j1)[d];
		j2t=INTEGER(j2)[d];
		
		
		/* Plus petit diamètre: pour les semaines suivant la première,
		   initialisation des effectifs du cluster correspondant à ce diamètre
		   sur le plus petit diamètre de la semaine d'avant. Pour 
		   la première semaine, les effectifs des semaines d'avant sont nuls
		   (initialisation hors boucle, donc c'est bon).
		*/
		if ((t>0)&&(d==0)) {
		    REAL(nz)[d+t*Ndiam] = REAL(nz)[d+(t-1)*Ndiam];
		    REAL(muz)[d+t*Ndiam] = REAL(muz)[d+(t-1)*Ndiam];
		}
		
		/* Pour toutes les semaines, pour les diamètres autres que le plus petit,
		   on initialise l'effectif sur le diamètre précédent de la même semaine
		*/
		if (d>0) {
		    REAL(nz)[d+t*Ndiam] = REAL(nz)[(d-1)+t*Ndiam];
		    REAL(muz)[d+t*Ndiam] = REAL(muz)[(d-1)+t*Ndiam];
		} else {
		    /* Et si l'on commence par le plus petit diamètre, on commence par 
		       ajouter les valeurs du pixel courant aux effectifs déjà présents */
		    REAL(nz)[d+t*Ndiam] = REAL(nz)[d+t*Ndiam] + REAL(VECTOR_ELT(donnees, t))[i];
		    REAL(muz)[d+t*Ndiam] = REAL(muz)[d+t*Ndiam]+ REAL(VECTOR_ELT(theorique, t))[i];
		}
			    
		/* Et celles des voisins pour le diamètre courant */
		for (v = j1t; v<=j2t; v++) {
		    REAL(nz)[d+t*Ndiam] = REAL(nz)[d+t*Ndiam] + REAL(VECTOR_ELT(donnees, t))[INTEGER(vec)[v]];
		    REAL(muz)[d+t*Ndiam] = REAL(muz)[d+t*Ndiam] + REAL(VECTOR_ELT(theorique, t))[INTEGER(vec)[v]];
		}
		/* Calcul du rapport de vraisemblance */
		REAL(lam)[d+t*Ndiam] = rapportVraisemblanceKulldorf(REAL(nz)[d+t*Ndiam],
								    REAL(muz)[d+t*Ndiam], N);
		
	    }
	}

	/* Période et diamètres les plus probables pour cette valeur de i */
	vma = REAL(lam)[0]; /* Plus forte valeur pour le rapport de vraisemblance */
	dma = 0; /* diamètre correspondant */
	tma = 0; /* semaine correspondante */
	for (d=0; d<Ndiam; d++) {
	    for (t=0; t<Nt; t++) {
		if (REAL(lam)[d+t*Ndiam] > vma) {
		    vma = REAL(lam)[d+t*Ndiam];
		    tma = t;
		    dma = d;
		}
	    }
	}

	/* stockage */
	REAL(lambda)[i] = vma;
	INTEGER(Dc)[i] = dma+1; /* indice dans R */
	INTEGER(t1c)[i] = tma+1; /* indice dans R */
    }

    /* Stockage des résultats */
    SET_VECTOR_ELT(resu, 0, lambda);
    SET_VECTOR_ELT(resu, 1, Dc);
    SET_VECTOR_ELT(resu, 2, t1c);
    
    /* Sorties */
    UNPROTECT(7);
    return(resu);
}

