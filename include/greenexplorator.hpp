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
#ifndef GREENEXPLORATOR
#define GREENEXPLORATOR

#include <boost/unordered_map.hpp>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <limits.h>
#include "Exception.h"
#include "enum.hpp"
#include "principalExtension.hpp"
#include <stdint.h>
#include <map>
#include <sstream>
#include "greenSizeHash.hpp"


class GreenExplorator
{
	public:
		GreenExplorator();
		~GreenExplorator();
		void printArbre();
		int insertInList(PrincipalExtension &pe);
		int insertInList(PrincipalExtension &pe, std::list<PrincipalExtension> &c);
		void clearC();
		void greenExploration(PrincipalExtension);
        inline void setIsomorphTest(bool value) { isomorphTest = value;};
        inline void setDumpTruncated(bool value) { dumpTruncated = value;};
        inline void setP(mpz_class value) { p = value;};
        inline void setMaxDepth(int value) { max_depth = value;};
		
	protected:
		Carquois *carquois;
		std::list<PrincipalExtension> c;
		std::list<PrincipalExtension> cemetary;
		int generateMutations(PrincipalExtension &pe);
	private:
		uint64_t numGreen;
		uint64_t minLength;
		uint64_t maxLength;
		uint64_t infCut;
		uint64_t depthCut;
		mpz_class p;
		int max_depth;
        int truncated;
        bool isomorphTest;
        bool dumpTruncated;
		bool myIsomorphismNauty(PrincipalExtension &a, PrincipalExtension &b);
        std::map<uint64_t,mpz_class> sizes;	
        GreenSizeHash gsh;
       
		
};

#endif
