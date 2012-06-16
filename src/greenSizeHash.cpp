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

#include "greenSizeHash.hpp"
GreenSizeHash::GreenSizeHash() {}
GreenSizeHash:: ~GreenSizeHash() {}

void GreenSizeHash::increment(std::vector<int> &mutations, uint64_t size)
{
    std::map<uint64_t,mpz_class> temp;
    temp[size] = 1;
    addSizes(mutations,temp);
}

void GreenSizeHash::addSizes(std::vector<int> &mutations, std::map<uint64_t,mpz_class>& sizes)
{
    uint64_t i;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    std::map<uint64_t,mpz_class> temp;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::string mutations_str;
    strhash::iterator strhash_it;
    mul_hash::iterator mul_hash_it;
    std::pair<mul_hash::iterator, bool> ins_it;
    // 1. For each mutation subchain
    for(i=0;i<mutations.size();i++) {
        // Empty the ss object;
        ss.clear();ss.str("");
        ss << mutations_str << mutations[i] << " ";
        mutations_str = ss.str();
        // 2. Lookup the mutations_str 
        strhash_it = green_size.find(mutations_str);
        if(strhash_it != green_size.end())
        {
            // The subchain already has a multiplicities list
            temp = *(strhash_it->second);
            mul_hash_it = multiplicities.find(temp);
            mul_hash_it->second -= 1;
            if(mul_hash_it->second == 0)
            {
                multiplicities.erase(temp);
            }
        }
        // If the subchain does not already have a multiplicities subchain
        // We leave temp in an unitialized state
        // Do the math
        for(map_it=sizes.begin();map_it!=sizes.end();map_it++)
        {
            temp[map_it->first]+=map_it->second;
        }
        //3. Lookup the newly created map
        mul_hash_it = multiplicities.find(temp);
        if(mul_hash_it != multiplicities.end())
        {
            // The created map already exists !
            // a. increment refs
            mul_hash_it->second += 1;
            // b. update
            if(strhash_it != green_size.end())
            {
                strhash_it->second = &(mul_hash_it->first);
            }
            else
            {
                green_size[mutations_str] = &(mul_hash_it->first);
            }
        }
        else
        {
            // The created map does not yes exist !
            // a. insert
            ins_it = multiplicities.insert(std::make_pair(temp,1));
            // b. update
            green_size[mutations_str] = &(ins_it.first->first);
        } 
        temp.clear();
    }
}

uint64_t GreenSizeHash::GreenSizesSetSize(std::string s)
{
    strhash::iterator strhash_it;
    strhash_it = green_size.find(s);
    if(strhash_it == green_size.end())
    {
        return 0;
    }
    else
    {
        return strhash_it->second->size();
    }

}

mpz_class GreenSizeHash::GreenSize(std::string s,uint64_t v)
{
    strhash::iterator strhash_it;
    strhash_it = green_size.find(s);
    if(strhash_it == green_size.end())
    {
        return 0;
    }
    else
    {
        return strhash_it->second->at(v);
    }

}
        
mpz_class GreenSizeHash::MultiplicitiesRefs(std::map<uint64_t,mpz_class> &rmap)
{
    mul_hash::iterator it;
    it = multiplicities.find(rmap);
    if(it == multiplicities.end())
    {
        return -1;
    }
    else
    {
        return it->second;
    }
}

void GreenSizeHash::dumpSizeToScreen()
{
    strhash::iterator strhash_it;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::map<uint64_t,mpz_class> tmp;
    
    for(strhash_it = green_size.begin(); strhash_it != green_size.end() ; strhash_it++)
    {
        std::cout << strhash_it->first << std::endl;
        std::cout << "\t";
        tmp = *(strhash_it->second);
        for(map_it=tmp.begin() ; map_it != tmp.end() ; map_it++)
        {
            std::cout << "(" << map_it->first << ";" << (map_it->second).get_str() << ")" ;
        }
        std::cout << std::endl;
    }

}


void GreenSizeHash::dumpMulToScreen()
{
    mul_hash::iterator mul_hash_it;
    std::map<uint64_t,mpz_class>::iterator map_it;
    std::map<uint64_t,mpz_class> tmp;
    
    for(mul_hash_it = multiplicities.begin(); mul_hash_it != multiplicities.end() ; mul_hash_it++)
    {
        tmp = mul_hash_it->first;
        for(map_it=tmp.begin() ; map_it != tmp.end() ; map_it++)
        {
            std::cout << "(" << map_it->first << ";" << (map_it->second).get_str() << ")" ;
        }
        std::cout <<std::endl;
        std::cout << "\t" << mul_hash_it->second << std::endl;
    }
}
