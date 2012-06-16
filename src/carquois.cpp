/*
 * Copyright (c) 2007-2012, Joris Calabrese,
 *                          Grégoire Dupont, 
 *                          Matthieu Pérotin
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

 
#include "carquois.hpp"
#include "Exception.h"
/*
But: Constructeur
Entrée: un entier n, nombre de sommets du graphe associé
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données de l'objet sont allouées, la matrice d'incidence est initialisée à 0, pour tout i, pour tout j
*/
Carquois::Carquois(int n)
{
    int i,j;
    this->M=(int **)malloc(n*sizeof(int *));
    for(i=0;i<n;i++)
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            M[i][j]=0;
        }        
    this->n=n;
    this->graphAJour=0;
    semifree=0;
    nbVoisinsMax = -1;
    connexe = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Constructeur par recopie
Entrée: une référence sur un objet carquois
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont sont recopiées
*/

Carquois::Carquois(const Carquois &ca)
{
    int i,j;
    n=ca.n;
    if(ca.semifree == 0)
    {
        this->M=(int **)malloc(n*sizeof(int *));
        for(i=0;i<n;i++)
            (this->M)[i]=(int *)malloc(n*sizeof(int));
        for(i=0;i<n;i++)
            for(j=0;j<n;j++)
                M[i][j]=ca.M[i][j];
        this->semifree = 0;
    }
    else
    {
        this->semifree = 1;
    }
    this->graphAJour=ca.graphAJour;
    this->mutations = ca.mutations;
    this->score = ca.score;
    if(ca.graphAJour)
    {
        for(i=0;i<2*n;i++)
            nautyGC[i]=ca.nautyGC[i];
    }
    
    this->nbVoisinsMax = ca.nbVoisinsMax;
    this->connexe = ca.connexe;
    this->nextI = ca.nextI;
    this->nextJ = ca.nextJ;
    
}

/*
But: Constructeur par recopie et augmentation
Entrée: une référence sur un objet carquois, un entier k
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont recopiées dans un carquois de taille n+k
    /!\     Le Carquois construit n'est pas connexe !
        Il est nécessaire d'appeler connecter(i) pour obtenir un carquois connexe obtenu en reliant le sommet i à
            TOUS les autres.
*/

Carquois::Carquois(const Carquois &ca, int k)
{
    int i,j;
    n=ca.n+k;
    if(ca.semifree == 0)
    {
        this->M=(int **)malloc(n*sizeof(int *));
        for(i=0;i<n;i++)
            (this->M)[i]=(int *)malloc(n*sizeof(int));
        for(i=0;i<ca.n;i++)
            for(j=0;j<ca.n;j++)
                M[i][j]=ca.M[i][j];
        for(i=0;i<n;i++)
        {
            M[i][n-1]=0;
            M[n-1][i]=0;
        }
        this->semifree = 0;
    }
    else
    {
        this->semifree = 1;
    }

    this->graphAJour=0;
    semifree=0;
    nbVoisinsMax = -1;
    connexe = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Constructeur à partir d'une matrice
Entrée: une matrice
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données et les données de l'objet sont initialisées
*/

Carquois::Carquois(int ** mat_carquois, int taille, int indice)
{
    int i,j;
    n=taille;
    this->M=(int **)malloc(n*sizeof(int *));
    
    for(i=0;i<n;i++)
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            M[i][j]=mat_carquois[i][j];
    this->n=n;
    this->graphAJour=0;
    this->genScore();
    semifree=0;
    nbVoisinsMax = -1;
    connexe = -1;
    nextI=1;
    nextJ=0;
}

/* But construire des carquois types
 *     Entrée: type
 *     A 0
    D 1
    E 2
    ATILDE 3
    DTILDE 4
    ETILDE 5
    SPORADIQUE 6
    UNAMED 7
    nbSommets: le nombre de sommets à construire
    Sortie: Néant
    Précondition: type d'un type défini cf ci-dessus et .h
    PostCondition: les structures de données et les données de l'objet sont désallouées
*/
Carquois::Carquois(int type, int nbSommets, int orientation)
{
    int i;
    
    // Dans tous les cas le graphe ne sera pas à jour
    this->graphAJour=0;
    this->semifree=0;
    switch(type)
    {
        case A:
            n=nbSommets;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
            
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-1;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            
        break;
        case D:
            n=nbSommets;
            if(n<4)
                throw Exception("ERROR: not enough vertices");
            n=nbSommets;
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-2;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            M[n-3][n-1]=1;
            M[n-1][n-3]=-1;
            
        break;
        case E:
            n=nbSommets;
            switch(n)
            {
                case 6:
                case 7:
                case 8:
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    for(i=0;i<n-2;i++)
                    {
                        M[i][i+1]=1;
                        M[i+1][i]=-1;
                    }
                    M[n-4][n-1]=1;
                    M[n-1][n-4]=-1;
                break;
                default:
                    throw Exception("ERROR: bad vertex number asked");
            }    
        break;
        case ATILDE:
            n=nbSommets+orientation;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
                
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<nbSommets;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            for(i=1;i<orientation;i++)
            {
                M[nbSommets+orientation-i][nbSommets+orientation-(i+1)]=1;
                M[nbSommets+orientation-(i+1)][nbSommets+orientation-i]=-1;
            }
            M[0][n-1] = 1;
            M[n-1][0]=-1;
            #ifdef DEBUG
            this->affiche();
            #endif
        break;
        case ATILDEALT:
            n=nbSommets+orientation;
            if(n<2)
                throw Exception("ERROR: not enough vertices");
                
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=0;i<n-1;i++)
            {
                M[i][i+1]=1-(2*(i%2));
                M[i+1][i]=-1+(2*(i%2));
            }
            M[n-1][0]=-1;
            M[0][n-1]=1;
            #ifdef DEBUG
            this->affiche();
            #endif
        break;
        case DTILDE:
            if(nbSommets<3)
                throw Exception("ERROR: not enough vertices");
            n=nbSommets+1;
            this->M=(int **)calloc(n,sizeof(int *));
            for(i=0;i<n;i++)
                (this->M)[i]=(int *)calloc(n,sizeof(int));
            for(i=1;i<n-2;i++)
            {
                M[i][i+1]=1;
                M[i+1][i]=-1;
            }
            M[n-3][n-1]=1;
            M[n-1][n-3]=-1;
            M[0][2]=1;
            M[2][0]=-1;
        break;
        case ETILDE:
            switch(nbSommets)
            {
                case 6:
                    n=nbSommets+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));

                    M[0][2] = -1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][5] = -1;
                    M[3][1] = 1;
                    M[3][5] = -1;
                    M[4][5] = 1;
                    M[4][6] = -1;
                    M[5][2] = 1;
                    M[5][3] = 1;
                    M[5][4] = -1;
                    M[6][4] = 1;                
                    break;
                case 7:
                    n=nbSommets+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    
                    M[0][2] = -1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][4] = -1;
                    M[3][1] = 1;
                    M[3][5] = -1;
                    M[4][2] = 1;
                    M[4][7] = -1;
                    M[5][3] = 1;
                    M[5][7] = -1;
                    M[6][7] = 1;
                    M[7][4] = 1;
                    M[7][5] = 1;
                    M[7][6] = -1;
                    

                    break;
                case 8:
                    n=nbSommets+1;
                    this->M=(int **)calloc(n,sizeof(int *));
                    for(i=0;i<n;i++)
                        (this->M)[i]=(int *)calloc(n,sizeof(int));
                    M[0][2] = -1;
                    M[0][7] = 1;
                    M[1][3] = -1;
                    M[2][0] = 1;
                    M[2][4] = -1;
                    M[3][1] = 1;
                    M[3][6] = -1;
                    M[4][2] = 1;
                    M[4][6] = -1;
                    M[5][6] = 1;
                    M[6][3] = 1;
                    M[6][4] = 1;
                    M[6][5] = -1;
                    M[7][0] = -1;
                    M[7][8] = 1;
                    M[8][7] = -1;
                    
                break;
                default:
                    throw Exception("ERROR: bad vertex number asked");
            }
        break;
        case SPORADIQUE:
            n=nbSommets;
            if (n==3)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                for(i=0;i<n-1;i++)
                {
                    M[i][i+1] = 2;
                    M[i+1][i] = -2;
                }
                M[n-1][0] = 2;
                M[0][n-1] = -2;
            }
            else if (n==4)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                M[0][1] = -1;
                M[1][0] = 1;
                
                M[0][3] = -1;
                M[3][0] = 1;
                
                M[3][1] = 1;
                M[1][3]= -1;
                
                M[2][1] = 1;
                M[1][2] = -1;
                
                M[2][3] = 1;
                M[3][2] = -1;
                
                M[0][2] = 2;
                M[2][0] = -2;
                
            }
            else
            {
                throw Exception("ERROR: SPORADIQUE asked but n != 3 or 4");
            }
            
        break;
        case UNAMED:
            n=nbSommets;
            if (n>=5)
            {
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
                for(i=0;i<n-2;i++)
                {
                    M[i][i+1]=1;
                    M[i+1][i]=-1;
                }
                M[n-3][n-1]=1;
                M[n-1][n-3]=-1;
                M[0][1]=2;
                M[1][0]=-2;
            }
            else
                throw Exception("number of vertices too small");
        break;
        case  E_ELIPTIQUE:
            if(nbSommets>=6 && nbSommets<9)
            {
                n=nbSommets+2;
                this->M=(int **)calloc(n,sizeof(int *));
                for(i=0;i<n;i++)
                    (this->M)[i]=(int *)calloc(n,sizeof(int));
            }
            
            if(nbSommets==6)
            {
                M[0][1]=1;
                M[1][0]=-1;
                
                M[1][2]=1;
                M[2][1]=-1;
                            
                M[2][3]=2;
                M[3][2]=-2;
                
                M[1][3]=-1;
                M[3][1]=1;
                
                M[3][4]=1;
                M[4][3]=-1;
                
                M[5][3]=-1;
                M[3][5]=1;
                
                M[4][2]=1;
                M[2][4]=-1;
                
                M[5][6]=1;
                M[6][5]=-1;
                
                M[4][7]=1;
                M[7][4]=-1;
                
                M[5][2]=1;
                M[2][5]=-1;
    
                
            }
            else if(nbSommets==7)
            {
                
                M[0][1]=1;
                M[1][0]=-1;
                
                M[1][2]=1;
                M[2][1]=-1;
                            
                M[2][3]=2;
                M[3][2]=-2;
                
                M[1][3]=-1;
                M[3][1]=1;
                
                M[3][4]=1;
                M[4][3]=-1;
                
                M[5][3]=-1;
                M[3][5]=1;
                
                M[4][2]=1;
                M[2][4]=-1;
                
                M[5][6]=1;
                M[6][5]=-1;
                
                M[6][7]=1;
                M[7][6]=-1;
                
                M[5][2]=1;
                M[2][5]=-1;
                
                M[0][8]=1;
                M[8][0]=-1;
            }
            else if(nbSommets==8)
            {
                M[0][1]=1;
                M[1][0]=-1;
                M[1][2]=1;
                M[1][3]=-1;
                M[2][1]=-1;
                M[2][3]=2;
                M[2][4]=-1;
                M[2][5]=-1;
                M[3][1]=1;
                M[3][2]=-2;
                M[3][4]=1;
                M[3][5]=1;
                M[4][2]=1;
                M[4][3]=-1;
                M[5][2]=1;
                M[5][3]=-1;
                M[5][6]=1;
                M[6][5]=-1;
                M[6][7]=1;
                M[7][6]=-1;
                M[7][8]=1;
                M[8][7]=-1;
                M[8][9]=1;
                M[9][8]=-1;
            }
            else
            {
                throw Exception("Eliptique asked but wrong vertex number (must be 6,7 or 8)");
            }
            break;
        default:
            throw Exception("Bad type");
    } // end switch
    this->genScore();
    nbVoisinsMax = -1;
    // Ces carquois sont connexes
    connexe = 1;
    nextI=1;
    nextJ=0;

}

Carquois::Carquois(const char *file)
{
    std::string contenu,ligne;
    std::ifstream f(file);
    boost::char_separator<char> sep(",[] \t;");
    std::vector<int> val;
    std::istringstream *iss;
    int i,j;
    unsigned int n;
    bool qmu = false;
    semifree=0;
    if(!f)
        throw Exception("ERROR: cannot open file !");
    std::getline (f, ligne);
    if(ligne == "//Number of points")
    {
        std::cout << "qmu file format detected !" << std::endl;
        qmu = true;
        while(ligne != "//Matrix")
            std::getline(f,ligne);
        std::getline(f,ligne);
        std::getline(f,ligne);
    }
    
    while(!f.eof())
    {
        
        contenu+=ligne;
        std::getline(f,ligne);
        // Keller's new file format ends matrix definition by "Traffic lights"
        // the old file format ends it with "Points"
        // We keep both tests for retro compatibility
        if(qmu && (ligne == "//Traffic lights" || ligne == "//Points"))
        {
            break;
        }
        
    }
    tokenizer tokens(contenu, sep);
    for (tokenizer::iterator tok_iter = tokens.begin();
         tok_iter != tokens.end(); ++tok_iter)
    {
        iss = new std::istringstream(*tok_iter);
        *iss >> i;
        delete iss;
        val.push_back(i);
    }
    f.close ();
    n=(unsigned int)sqrt(val.size());
    if(n*n != val.size())
    {
        throw Exception("Bad file format !");
    }
    
    this->M=(int **)malloc(n*sizeof(int *));
    for(i=0;i<n;i++)
    {
        (this->M)[i]=(int *)malloc(n*sizeof(int));
    }
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            M[i][j]=val[i*n + j];
        }        
    }
    this->n=n;
    this->graphAJour=0;
    nbVoisinsMax = -1;
    connexe = -1;
    nextI=1;
    nextJ=0;
}

/*
But: Destructeur par défaut
Entrée: Néant
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données de l'objet sont désallouées
*/
Carquois::~Carquois()
{
    int i;
    if(semifree==0)
    {
        for(i=0;i<this->n;i++)
            free((this->M)[i]);
        free(this->M);
    }
}

void Carquois::semiDestroy()
{
    int i;
    if(semifree==0)
    {
        getNbVoisinsMax(); // Penser à faire ça... sinon...
        genScore();
        estConnexe();
        semifree=1;
        for(i=0;i<this->n;i++)
            free((this->M)[i]);
        free(this->M);
    }
}
/*
But: Constructeur par défaut
Entrée: Néant
Sortie: Néant
Précondition: Néant
PostCondition: les structures de données de l'objet sont allouées, la matrice d'incidence est initialisée avec une instance par défaut à 9 sommets

*/
Carquois::Carquois()
{

}

/*
But: Afficher la matrice d'incidence sur la sortie standard
Entrée: Néant
Sortie: Néant
Précondition: Néant
PostCondition: La matrice d'incidence est affichée sur la sortie standard
*/
void Carquois::affiche()
{
    int i,j;
    if(semifree==0)
    {
        for(i=0;i<this->n;i++)
        {
            for(j=0;j<this->n;j++)
            {
                std::cout << M[i][j] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    else
        std::cout << "The quiver has been semiFreed" << std::endl;
}

/*
But: Appliquer la fonction de mutation sur un sommet du graphe
Entrée: en entier k correspondant à un des sommets du graphe
Sortie: Néant
PréCondition: k est un sommet du graphe (=> k>=0 et k<n)
PostCondition: La fonction \mu_k est appliquée au carquois
*/

void Carquois::mutate(int k)
{
    int i,j;
    int dernierElement;
    // On ne fait rien si k ne correspond pas à un sommet du graphe
    if(k<0 || k>= this->n)
        return;
        
    //Création d'une matrice temporaire, copie de la matrice d'incidence
    /*int Mp[this->n][this->n];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            Mp[i][j]=M[i][j];
    */        
    // Application de la fonction mu        
    
    for(i=0;i<n;i++)
    {
        if (i==k) continue;
        for(j=0;j<n;j++)
        {
            if (j==k) continue;
            M[i][j] = M[i][j] + (valeurAbs(M[i][k])*M[k][j] + M[i][k]*valeurAbs(M[k][j]))/2;
        }
    }
    for(i=0;i<n;i++)
    {
        M[i][k]=-M[i][k];
        M[k][i]=-M[k][i];
    }

    if(this->infinite())
        throw Exception("Mutation class is infinite ! " + getMutations());
    if(this->graphAJour)
    {
        this->graphAJour=0;
    }

    /*
        On met à jour les mutations qui ont été appliquées sur le carquois
        Si la mutation appliquée est la même que la dernière qui avait été appliquée:
            alors on a appliqué deux fois la même, ce qui revient à ne pas l'appliquer, on l'efface de la liste
        Sinon
            on ajoute la mutation à la liste des mutations déjà appliquées
        
    */
    
    if(!mutations.empty())
    {
        dernierElement = mutations.back();
        if(dernierElement == k)
        {
            mutations.pop_back();
        }
        else
        {
            mutations.push_back(k);
        }
    }
    else
    {
        mutations.push_back(k);
    }
    this->genScore();
}

/*
But: Calculer la valeur absolue d'un entier
Entrée: un entier k
Sortie: un entier s
Précondition: Aucune
PostCondition: s <- |k|
*/

int Carquois::valeurAbs(int k)
{
    if(k < 0)
        return -k;
    else
        return k;
}

/*
But: Setter pour la matrice d'incidence
Entrée: 3 entiers i, j et val
Sortie: Aucune
Précondition: i et j compris entre 0 et n-1
PostCondition: M[i][j] <- val
*/
void Carquois::setM(int i, int j, int val)
{
    if(i<n && j < n && i>=0 && j>=0)
    {
        M[i][j]=val;
        if(this->graphAJour)
        {
            this->graphAJour = 0;
        }
        if(this->connexe)
        {
            this->connexe = -1;
        }
    }
    else
        throw Exception("DOMAIN_ERROR: setM");
    
    
}
/* Cette fonction génère les structures qui vont bien pour les appels à Nauty */
void Carquois::genGraph()
{
    int i,j,m,nbSommetsNauty;
    int lab1[MAXN],ptn[MAXN],orbits[MAXN];
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[5*MAXM];
    
    if(!this->graphAJour)
    {
        
        nbSommetsNauty = 2 * this->getN();
        m=(nbSommetsNauty + WORDSIZE - 1)/WORDSIZE;

        /* Si on trouve une valeur strictement positive dans la matrice d'incidence, alors on ajoute une arrête dans notre graphe */
        for(i=0;i<this->getN();i++)
        {
            gv=GRAPHROW(nautyG,i+this->getN(),m);
            EMPTYSET(gv,m);
            
            gv=GRAPHROW(nautyG,i,m);
            EMPTYSET(gv,m);
            /* On ajoute les fausses arrêtes entre le layer 0 et le layer 1 */
            ADDELEMENT(gv,i+this->getN());
            for(j=0;j<this->getN();j++)
            {
                /* multiplicité de 1 */
                if(this->getM(i,j)==1)
                {
                    gv=GRAPHROW(nautyG,i,m);
                    ADDELEMENT(gv,j);
                }
                else
                {
                    if(this->getM(i,j)==2)
                    {
                        gv=GRAPHROW(nautyG,i+this->getN(),m);
                        ADDELEMENT(gv,j+this->getN());
                    }
                }
            }
        }
        options.getcanon = TRUE;
        options.digraph = TRUE;
        options.defaultptn = FALSE;
        nauty_check(WORDSIZE,m,nbSommetsNauty,NAUTYVERSIONID);
        
        for(i=0;i<2*n;i++)
        {
            lab1[i]=i;
            ptn[i]=1;
        }
        ptn[n-1]=0;
        ptn[2*n-1]=0;
        
        
        nauty(nautyG,lab1,ptn,NULL,orbits,&options,&stats,
                                  workspace,5*MAXM,m,nbSommetsNauty,nautyGC);
        this->graphAJour=1;    
    }
}

/* Cette fonction teste deux propriétés du Carquois:
 *
 * - Si le carquois a une arrête au moins triple alors il n'est pas de mutation finie (Cas 0)
 * - Si un sommet x est connecté à un sommet v par une arrête double, et si v est
 *    aussi connecté par une arrête double à un autre sommet, alors le carquois n'est pas de
 *    mutation finie (Cas 1)
 *
 * Cette fonction renvoie vrai si l'un des deux cas présents ci-dessus est trouvé, et faux sinon.
 * Attention: ce n'est pas parceque cette fonction renvoie faux que le carquois est nécessairement de mutation finie !
 */
bool Carquois::infinite()
{
    int i,j,compteur;
    if(this->getN()>2)
    {
        for(i=0;i<this->getN();i++)
        {
            compteur = 0;
            for(j=0;j<this->getN();j++)
            {
                if(this->getM(i,j) >=2 || this->getM(i,j) <=-2)
                {
                    compteur++;
                }
                if(compteur == 2) /* Cas 1 le sommet i est connecté à 2 autres sommets par une arrête double */
                {
                    return true;
                }
                if(this->getM(i,j) >=3) /* Cas 0 */
                {
                    return true;
                }
            }
        }
    } 
    return false;
}



graph *Carquois::getNautyGraph()
{
    if(!this->graphAJour)
    {
        this->genGraph();
    }
    return  (graph *)&nautyGC;
}


void Carquois::printMutations()
{
    std::vector<int>::iterator i;
    if(mutations.empty()) std::cout << "-\n";
    else
    {
        for(i=mutations.begin();i!=mutations.end();i++)
        {
            std::cout << *i;
            std::cout <<".";
        }
    std::cout << "\n";
    }
}

/* Calcule le "score" d'un carquois, utile pour sortir un "bon" représentant
 * du groupe de mutation
 */
void Carquois::genScore()
{
    int i,j;
    score=0;
    
    #ifdef SCORE1
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (M[i][j]>0)
            {
                score-=M[i][j];
            }
    #else
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if (M[i][j]==2)
            {
                score-=1;
            }
    #endif
}

/* Effectue mutations mutations et regarde si le degré des arcs explose
 * 
 */
bool Carquois::testInfiniEmpirique(int mutations)
{
    int i;
    srand(time(NULL));
    // on travaille sur une copie
    Carquois t = *this;
    if(t.infinite())
        return true;
    for(i=0;i<mutations;i++)
    {
        t.mutate((int) (n * (rand() / (RAND_MAX + 1.0))));
        if(t.infinite())
            return true;
    }
    return false;
}

void Carquois::toFile(const char* filename)
{
    int i,j;
    std::ofstream fichierSortie(filename);
    if(!fichierSortie)
        throw Exception("ERROR: cannot open output file !");
    fichierSortie << "[";
    for(i=0;i<this->getN();i++)
    {
        fichierSortie << "[";
        for(j=0;j<this->getN();j++)
        {
            fichierSortie << this->M[i][j];
            if(j!=this->getN() - 1)
                fichierSortie << ",";
        }
        fichierSortie << "]"  ;
        if(i!=this->getN()-1)
            fichierSortie << "," << std::endl;
    }
    fichierSortie << "]" << std::endl;
    fichierSortie.close();
}

std::string Carquois::getMutations()
{
    std::vector<int>::iterator i;
    std::stringstream out;

    if(mutations.empty()) out << "-\n";
    else
    {
        for(i=mutations.begin();i!=mutations.end();i++)
        {
            out << *i;
            out << ".";
        }
    out << "\n";
    }
    return out.str();
}

/* Renvoie la valence maximale du carquois */

int Carquois::getNbVoisinsMax()
{
    int i,j;
    int nbVoisinsTemp;
    if(nbVoisinsMax == -1)
    {
        for(i=0;i<n;i++)
        {
            nbVoisinsTemp=0;
            for(j=0;j<n;j++)
            {
                if(this->M[i][j] != 0)
                {
                    nbVoisinsTemp += 1;
                }    
            }
            if(nbVoisinsTemp > nbVoisinsMax)
            {
                nbVoisinsMax = nbVoisinsTemp;
            }
        }    

    }
    return nbVoisinsMax;
}
/* Algorithme de Warshall, calcule la fermeture transitive du graphe*/
int Carquois::estConnexe()
{
    bool mat[n][n];
    int i,j,k;
    if(connexe != -1)
        return connexe;
    // On travaille sur une copie du graphe */
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if(M[i][j] != 0)
                mat[i][j] = true;
            else
                mat[i][j] = false;
#ifdef DEBUG
    for(i=0;i<this->n;i++)
    {
        for(j=0;j<this->n;j++)
        {
            std::cout << mat[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
#endif

    for(i = 0; i < n; i++)
       for(j = 0; j < n; j++)
          for(k = 0; k < n; k++)
             mat[j][k] = mat[j][k] || (mat[j][i] && mat[i][k]);
#ifdef DEBUG
    for(i=0;i<this->n;i++)
    {
        for(j=0;j<this->n;j++)
        {
            std::cout << mat[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";
#endif
/* Si la fermeture transitive est une clique, le graphe est connexe */
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            if(mat[i][j]!=true)
            {
                connexe = 0;
                return 0;
            }
    connexe = 1;
    return 1;
}

/* But: Dire si un sommet i a une arrête double ou non
 * Entrée: i un sommet du carquois
 * Sortie: vraie si i est connecté par une arrête double, faux sinon
 * Précondition: i est compris entre 0 et n
 * PostCondition: néant
 */
bool Carquois::aUneDouble(int sommet)
{
    int i;
    if(sommet < n && sommet >= 0)
    {
        for(i=0;i<n;i++)
        {
            if(M[i][sommet] == 2 || M[i][sommet] == -2)
                return true;
        }
        return false;
    }
    else
    {
        throw new Exception("aUneDouble: the argument is not a vertex");
    }
}

/*
But: Dire si trois sommets du graphe forment un 3-cycle orienté ou non
Entrée: 3 entiers correspondant à 3 sommets du carquois
Sortie: Vrai si les trois sommets forment un 3-cycle orienté, Faux sinon.
Précondition: i, j et k sont compris entre 0 et n et sont tous les trois distincts
PostCondition: néant
*/
bool Carquois::troisCycleOriente(int i, int j, int k)
{
    /* On vérifie que les entrées sont bien des sommets du graphe */
    if(i<n && j < n && k < n && i>=0 && j>=0 && k>=0)
    {
        if(i != j && i != k && j != k)
        {
            if( M[i][j] != 0 && M[j][k] != 0 && M[k][i] != 0 )
            {
                if( (M[i][j] > 0 && M[j][k] > 0 && M[k][i] > 0) || (M[i][j] < 0 && M[j][k] < 0 && M[k][i] < 0) )
                    return true; // Toutes les arrêtes sont dans le même sens
                else
                    return false;
            }
            else
            {
                // Les trois sommets ne forment pas un cycle !
                return false;
            }
        }
        else
        {
            // Les trois sommets ne sont pas différents !
            return false;
        }
    }
    else
    {
        throw new Exception("ERROR, Carquois::troisCycleOriente: One of the three arguments is not a vertex !");
    }
}

bool Carquois::cyclique()
{
    int visite[n];
    int i;
    for(i=0;i<n;i++)
    {
        visite[i] = 0;
    }
    for(i=0;i<n;i++)
    {
        if(visite[i] == 0) {
            if (exploreCycle(visite,i)) {
                return true;
            }
        }
    }
    return false;
}

bool Carquois::exploreCycle(int *visite,int i)
{
    int j;
    visite[i] = 1;
    for(j=0;j<n;j++)
    {
        if(M[i][j] > 0)
        {
            if(visite[j] == 1) {return true;}
            else 
            { 
                if (visite[j] == 0) 
                {
                    if(exploreCycle(visite,j))
                    {
                        return true;
                    }
                }
            }
        }
    }
    visite[i] = 2;
    return false;
}

std::vector<int> Carquois::getVoisins(int sommet)
{
    std::vector<int> voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] != 0)
            voisins.push_back(i);
    }
    return voisins;
}

std::vector<int> Carquois::getVoisinsDoubles(int sommet)
{
    std::vector<int> voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == 2 || M[i][sommet] == -2)
            voisins.push_back(i);
    }
    return voisins;
}

std::vector<int> Carquois::getVoisinsSimples(int sommet)
{
    std::vector<int> voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == 1 || M[i][sommet] == -1)
            voisins.push_back(i);
    }
    return voisins;
}

std::vector<int> Carquois::getVoisinsSimplesPredecesseurs(int sommet)
{
    std::vector<int> voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == -1)
            voisins.push_back(i);
    }
    return voisins;
}

std::vector<int> Carquois::getVoisinsSimplesSuccesseurs(int sommet)
{
    std::vector<int> voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == 1)
            voisins.push_back(i);
    }
    return voisins;
}

int Carquois::getNbVoisinsSimplesPredecesseurs(int sommet)
{
    int voisins;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == -1)
            voisins++;
    }
    return voisins;
}

int Carquois::getNbVoisinsSimplesSuccesseurs(int sommet)
{
    int voisins=0;
    int i;
    for(i=0;i<n;i++)
    {
        if(M[i][sommet] == 1)
            voisins++;
    }
    return voisins;
}

std::vector<int> Carquois::getSommetsArreteDoubleEntrante()
{
    std::vector<int> res;
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(M[i][j]==2)
            {
                res.push_back(j);
                break;
            }    
        }
    }
    return res;
}
std::vector<int> Carquois::getSommetsPasDArreteDouble()
{
    std::vector<int> res;
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(M[i][j]==2 || M[i][j] == -2)
            {
                break;
            }    
        }
        res.push_back(i);
    }
    return res;
}

/* Précondition: i est un sommet pourvu d'une arrête double entrante
 * Renvoie l'extrémité sortante de l'arrête
 * Renvoie -1 en cas d'erreur
 */
int Carquois::getSommetOrigineArreteDouble(int i)
{
    int j;
    for(j=0;j<n;j++)
    {
        if(M[j][i] == 2)
            return j;
    }
    return -1;
}
