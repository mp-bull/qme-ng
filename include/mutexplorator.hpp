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
#ifndef MUTEXPLORATOR
#define MUTEXPLORATOR

#include <iostream>
#include <utility>
#include <vector>
#include <limits.h>
#include "Exception.h"
#include "carquois.hpp"

typedef int ** matrice;

typedef struct infoMat derMat;
struct infoMat{
	matrice *matriceC;
	int taille;
	int sommet;
	int carquois;
};

class MutExplorator
{
	public:
		MutExplorator();
		~MutExplorator();
		virtual int isomorphismExplorator(Carquois depart,int granularite) = 0;
		void printArbre();
		void insertInList(Carquois *carquois);
		bool estDansC(Carquois *carquois);
		Carquois getRepresentant();
		void clearC();
		int getModeCmp();
		inline int getCSize()
		{
			return c.size();
		}
		
	protected:
		Carquois *carquois;
		std::vector<Carquois> c;
		void generateMutations(Carquois carquois);
		/* Fonction utilisant l'algorithme de Brendam McKay */
		bool myIsomorphismNauty(Carquois *a, Carquois *b);
		int modeComparaison;
		
};

#endif
