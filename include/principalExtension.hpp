/*
 * Copyright (c) 2011-2012, Grégoire Dupont, Matthieu Pérotin
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
#ifndef PRINCIPALEXT_H
#define PRINCIPALEXT_H

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>
#include <boost/tokenizer.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdint.h>
#include <gmpxx.h>
#include "Exception.h"
#include "carquois.hpp"
#include "nautinv.h"


typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

class PrincipalExtension
{
	public:
		PrincipalExtension(Carquois c);
		PrincipalExtension(const char* file);
		~PrincipalExtension();
		PrincipalExtension(const PrincipalExtension &ca);
		int mutate(int k, mpz_class p);
		void affiche();
		bool infinite(mpz_class p);
		void toFile(const char* filename);
		void printMutations(int s);
		void printMutationsE(int s);
		void genGraph();

		/* Getters et Setters */
		void setM(int i, int j, mpz_class val);
		/*
		But: Getter pour la matrice d'incidence
		Entrée: 2 entiers i et j
		Sortie: un entier
		Précondition: i et j compris entre 0 et n-1
		PostCondition: return M[i][j] si i et j conformes, sinon jette une exception

		*/
		inline mpz_class getM(int i, int j)
		{
			if(i<n && j < n && i>=0 && j>=0)
				return M[i*n+j];
			else
			{
				throw Exception("ERREUR_DOMAINE: getM");
			}
		}
		inline int getN()
		{
			return n;
		}
		inline int getNbSommetsNauty()
        {
            return nbSommetsNauty;
        }
		inline int lastMutation()
		{
			return (mutations.size()>0)?mutations[mutations.size()-1]:-1;
		}
		void generateSommetsVerts();
		int getNextSommetVert();
		int getRandomGreenVertex();
		void forceSommetVert(int s);
		inline bool getGraphAJour()
		{
			return graphAJour;
		}
		inline void setMultiplicity(uint64_t size, mpz_class mul)
		{
			this->multiplicity[size]=mul;
		}

		inline mpz_class getMultiplicity(uint64_t size)
		{
			if(this->multiplicity.find(size) == this->multiplicity.end())
            {
                return 0;
            }
            else
            {
			    return this->multiplicity[size];
            }
		}

		inline mpz_class incMultiplicity(uint64_t size)
		{
            this->multiplicity[size] += 1;
			return this->multiplicity[size];
		}

		inline mpz_class addMultiplicity(uint64_t size,mpz_class value)
		{
            this->multiplicity[size] += value;
			return this->multiplicity[size];
		}
		void addMultiplicity(PrincipalExtension &p);
        inline std::map<uint64_t, mpz_class> *getMultiplicityMap()
        {
            return &multiplicity;
        }
        
        inline std::map<mpz_class, mpz_class> *getMultiplicitiesMap()
        {
            return &multiplicities;
        }

        inline std::vector<int> getMutations(void)
        {
            return mutations;
        }
        inline std::string getMutationsString(void)
        {
            return mutationString;
        }
        inline uint64_t getMutationsSize(void)
        {
            return mutationsSize;
        }
		Carquois *getCarquois(void);
		graph *getNautyGraph();
        std::string mutationsToString();
        void shiftMultiplicities();
        void unshiftMultiplicities();
        void semiDestroy();
	private:
		mpz_class *M;
		int n;
		int nbSommetsNauty;
		std::map<uint64_t,mpz_class> multiplicity;
        std::map<mpz_class,mpz_class> multiplicities;
		std::vector<int> mutations;
		std::vector<int> sommetsVerts;
        graph nautyG[MAXN*MAXM];
		graph nautyGC[MAXN*MAXM];
		set *gv;
		bool graphAJour;
        bool semiFreed;
        std::string mutationString;
        uint64_t mutationsSize;

};
#endif
