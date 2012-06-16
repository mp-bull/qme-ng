/*
 * Copyright (c) 2012, Grégoire Dupont, Matthieu Pérotin
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
#ifndef GREENSIZEHASH
#define GREENSIZEHASH

#include <boost/unordered_map.hpp>
#include <map>
#include <stdint.h>
#include <vector>
#include <gmpxx.h>
#include <iostream>
#include <string>
#include <sstream>
#include "mapHasher.hpp"


// Two typedefs to begin, the types manipulated
// by this class are kind of complicated...
typedef boost::unordered_map<std::string, const std::map<uint64_t,mpz_class> *> strhash;
typedef boost::unordered_map<std::map<uint64_t,mpz_class> , int, MapHasher<std::map<uint64_t,mpz_class> > > mul_hash;

class GreenSizeHash
{

    public:
        GreenSizeHash();
        ~GreenSizeHash();
        void increment(std::vector<int> & m, uint64_t s);
        void addSizes(std::vector<int> & m, std::map<uint64_t, mpz_class> &s);
        void dumpSizeToScreen();
        void dumpMulToScreen();

        uint64_t GreenSizesSetSize(std::string);
        mpz_class GreenSize(std::string,uint64_t);
        
        mpz_class MultiplicitiesRefs(std::map<uint64_t,mpz_class> &);

        inline strhash getGreen_size() {return green_size;}
        inline std::map<uint64_t,mpz_class> *getGreenSizes(std::string s) {return const_cast<std::map<uint64_t,mpz_class> *>(green_size[s]);}
        inline mul_hash getMultiplicities() {return multiplicities;}




    private:
        strhash green_size;
        mul_hash multiplicities;
};
#endif
