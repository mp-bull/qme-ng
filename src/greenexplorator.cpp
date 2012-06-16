/*
 * Copyright (c) 2007-2012, Grégoire Dupont, Matthieu Pérotin
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
#include "greenexplorator.hpp"

GreenExplorator::GreenExplorator()
{
    numGreen=0;
    minLength=std::numeric_limits<uint64_t>::max();
    maxLength=0;
    truncated = 0;
    isomorphTest = true;
    dumpTruncated = false;
    infCut = 0;
    depthCut = 0;
}
GreenExplorator::~GreenExplorator()
{
}
void GreenExplorator::printArbre()
{
    std::list<PrincipalExtension>::iterator i;
    int j=0;
    for(i=c.begin();i!=c.end();i++)
    {
        std::cout << j++ << ":";
        i->printMutations(0);
    }
}

void GreenExplorator::clearC()
{
    c.clear();
}
        
int GreenExplorator::generateMutations(PrincipalExtension &pe)
{
    int i;
    int ret;
    int sommet;
    int create=0;
    unsigned int size;
    std::string mutations_str="";
    std::string filename = "";
    std::vector<int> mutations_v;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    PrincipalExtension p = pe;
    strhash::iterator iit;
    std::map<uint64_t,mpz_class>::iterator mult_it;
    std::map<uint64_t,mpz_class> *mult;
    #ifdef DEBUG
    std::cout << "Travail avec "; pe.printMutationsE(0);
    #endif
    sommet = pe.getNextSommetVert();
    if(sommet == -1) {
        return 1;
    }
    if(sommet==p.lastMutation())
    {
        // No need to mutate twice on the same vertex
        return 1;
    }
    if(p.getMutationsSize() >= this->max_depth) {
        return 4;
    }
    ret = p.mutate(sommet, this->p);
    if (ret == 0) {
        // if mutate returned 0, then infinity was detected
        return 3;
    }
    if (ret == 1) {
        // a Green quiver was detected !
        size = p.getMutationsSize();
        mult = p.getMultiplicityMap();
        for(mult_it = mult->begin(); mult_it!=mult->end(); mult_it++) { 
            sizes[mult_it->first] += mult_it->second;
        }
        if (isomorphTest) {
            mutations_v = p.getMutations();
            gsh.increment(mutations_v,size);
        }
        if(size > maxLength) { 
            maxLength = size; 
            std::cerr << "M:";
            p.printMutationsE(1);
            std::cerr << "Q:";
            p.printMutationsE(0);
        }
        if(size < minLength) { 
            minLength = size; 
            std::cerr << "m:";
            p.printMutationsE(1);
            std::cerr << "Q:";
            p.printMutationsE(0);
        }

        numGreen+=1;
        if((numGreen % 1000000) == 0) {
            std::cerr << "S (" << numGreen << "):";
            p.printMutationsE(0);
        }
        if(isomorphTest) {insertInList(p,cemetary);}     
        return 0;
    }
    if (isomorphTest) {
        return insertInList(p);
    }
    else {
        c.push_back(p);
        return 2;
    }
    #ifdef DEBUG
    std::cout << "Fin du travail avec "; p.printMutations(0);    
    #endif
}

int GreenExplorator::insertInList(PrincipalExtension &pe)
{
    std::list<PrincipalExtension>::iterator ri; 
    std::list<PrincipalExtension>::reverse_iterator rxi;
    std::map<uint64_t,mpz_class>::iterator it_mul;
    std::map<uint64_t,mpz_class> *mul_map;
    std::map<uint64_t,mpz_class> tmp;
    std::map<uint64_t,mpz_class> *green_sizes;
    std::map<uint64_t,mpz_class>::iterator it;
    strhash::iterator iit;
    std::string mutations_str="";
    std::vector<int> mutations_v;
    int i;
    uint64_t pe_size,n_size;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    // ri is an  iterator, it browses the list from the beginning
    for(ri=c.begin();ri!=c.end();ri++)
    {
        if(this->myIsomorphismNauty(pe,*ri))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "est isomorphe à  "; 
                (*ri).printMutations(0); 
                std::cout << "\n";
            #endif
            // This principal extension is still to be considered
            // increase its multiplicity
            (*ri).addMultiplicity(pe);
            break;
        }
    }
    if( ri != c.end()) { return 0;}
    for(rxi=cemetary.rbegin();rxi!=cemetary.rend();rxi++)
    {
        if(this->myIsomorphismNauty(pe,*rxi))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "est isomorphe à (C) "; 
                (*rxi).printMutations(0); 
                std::cout << "\n";
            #endif
            // This quiver is isomorph to an already considered
            // quiver

            mutations_str = (*rxi).getMutationsString();
            if(gsh.GreenSizesSetSize(mutations_str) != 0) {
                // The quiver is isomorph to a quiver leading to a green max suite
                // Update de length list !

                // 2. For all the quiver sizes attainable with the quiver
                green_sizes=gsh.getGreenSizes(mutations_str);
                for(it=green_sizes->begin();it!=green_sizes->end();it++) {
                    // For all the sizes multiplicities
                    mul_map = pe.getMultiplicityMap();
                    for(it_mul=mul_map->begin();it_mul != mul_map->end();it_mul++) {
                        pe_size=it_mul->first;
                        n_size=it->first-(*rxi).getMutationsSize()+pe_size;
                        sizes[n_size]+=pe.getMultiplicity(it_mul->first)* it->second;
                        #ifdef DEBUG
                        std::cout << "Add size: " <<  n_size << 
                                     " Mul " << "("<<pe.getMultiplicity(it_mul->first) << "*" <<it->second<<") ="
                                             << pe.getMultiplicity(it_mul->first) *it->second 
                                             << std::endl;
                        #endif
                    }
                    n_size=it->first-(*rxi).getMutationsSize()+pe.getMutationsSize();
                    tmp[n_size] = it->second;
                }
                
                // 3. Update the quiver list
                mutations_v = pe.getMutations();
                gsh.addSizes(mutations_v,tmp);
            }
            break;
            
        }
    }
    // if ri went all the way through the end, then the quiver is not
    // ismomorph to any quiver in already in the list, so we add it
    if(rxi==cemetary.rend())
    {
        c.push_back(pe);
        #ifdef DEBUG
            std::cout << "Ajout de";
            pe.printMutations(0);
        #endif
        // Insertion done, return 2;
        return 2;
        
    }
    // No insertion done, return 0
    return 0;
}

int GreenExplorator::insertInList(PrincipalExtension &pe, std::list<PrincipalExtension> &c)
{
    std::list<PrincipalExtension>::iterator ri; 
    std::list<PrincipalExtension>::reverse_iterator rxi;

    // ri is an  iterator, it browses the list from the beginning
    for(ri=c.begin();ri!=c.end();ri++)
    {
        if(this->myIsomorphismNauty(pe,*ri))
        {
            #ifdef DEBUG
                pe.printMutations(0); 
                std::cout << "est isomorphe cim à  "; 
                (*ri).printMutations(0); 
                std::cout << "\n";
            #endif
            break;
            
        }
    }
    // if ri went all the way through the end, then the quiver is not
    // ismomorph to any quiver in already in the list, so we add it
    if(ri==c.end())
    {
        pe.semiDestroy();
        c.push_back(pe);
        #ifdef DEBUG
            std::cout << "Ajout de (C)";
            pe.printMutations(0);
        #endif
        // Insertion done, return 2;
        return 2;
        
    }
    // No insertion done, return 0
    return 0;
}

void GreenExplorator::greenExploration(PrincipalExtension pe)
{
    int index = 0;
    int ret;
    std::map<uint64_t,mpz_class>::iterator it;
    mpz_class total=0;
    std::list<PrincipalExtension>::iterator pei;
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    std::string filename;
    uint64_t cutPending = 0;
    // Initial population of the list
    insertInList(pe);
    pei = c.begin();
    while(generateMutations(*pei) != 1) {pei=c.begin();};
    if (isomorphTest) {
        insertInList(*pei,cemetary);
    }
    c.erase(pei);
    pei = c.begin();
    for(index=c.size();index>=1;index--) {
        while(generateMutations(*pei) != 1) {pei=c.begin();};
        if (isomorphTest) {
            insertInList(*pei,cemetary);
        }
        c.erase(pei);
        pei = c.begin();
    }
    pei=c.end();pei--;
    // Main loop
    while(c.size()!=0) {
    #ifdef DEBUG
    std::cout << "C.size: " << c.size() << " Cem.size: " << cemetary.size() << "\t\t";
    std::cout << "Travail avec "; (*pei).printMutations(0);
    (*pei).affiche();
    #endif
        ret = generateMutations(*pei);
        switch(ret) {
            case 4:
                // Branch cut
            case 3:
                // Infinity detected on the branch
                // Cut the branch !
                if(dumpTruncated) {
                    ss.clear();
                    if(ret == 3) {
                        ss << "DInf_" << infCut << ".quiv";
                    }
                    if(ret == 4) {
                        ss << "DTrunc_" << depthCut << ".quiv";
                    }
                    ss >> filename;
                    (*pei).toFile(filename.c_str());
                }
                if(ret == 3) {
                    infCut++;
                }
                if(ret == 4) {
                    this->truncated = 1;
                    depthCut ++;
                }
                if(isomorphTest) {
                    if(((*pei).getMultiplicityMap())->size() > 0) {
                        cutPending++; 
                    }
                }
                
            case 1:
                if (isomorphTest) {
                    insertInList(*pei,cemetary);
                }
                c.erase(pei);
            case 2:
                pei=c.end();pei--;
                break;
                // Nothing to be done for case 0
                // a green quiver was detected
                // and the list is not empty... some
                //mutations remains to be explored !
        }
    }
    // Print results
    for ( it=sizes.begin() ; it != sizes.end(); it++ ) {
        std::cout << (*it).first << "\t=>\t" << (*it).second << std::endl;
        total +=(*it).second;
    }
    std::cout << "Total: " << total << std::endl;
    if(this->truncated == 1) {
        std::cout << "Exploration Truncated at depth: " << max_depth << std::endl;
    }
    if(this->infCut > 0) {
        std::cout << "Num branches cut because of excessive Mul: " << infCut << std::endl;
    }
    if(this->depthCut > 0) {
        std::cout << "Num branches cut because of excessive depth: " << depthCut << std::endl;
    }
    if(cutPending > 0) {
        std::cout << "WARNING: "<< cutPending << " branches were cut with pending isomorphs." << std::endl;
        std::cout << "All max green suite < "<< max_depth << " may not have been found." << std::endl;
    }
}

bool GreenExplorator::myIsomorphismNauty(PrincipalExtension &a, PrincipalExtension &b)
{
    int i;
    int n = a.getN();
    int nbNautyVert;
    std::map<mpz_class, mpz_class> *mul_a;
    std::map<mpz_class, mpz_class> *mul_b;
    std::map<mpz_class, mpz_class>::iterator ita,itb,mulend;

    graph *c1;
    graph *c2;


    // These two calls must be placed before hand
    // They are responsible for initializing all the other variables of
    // objects (nbSommetsNauty, multiplicities...)
    c2 = (graph *)a.getNautyGraph();
    c1 = (graph *)b.getNautyGraph();

    // This is very different from the number of vertices
    // This is the number of vertices in the nauty graph
    // which depends on the edge multiplicities
    nbNautyVert = a.getNbSommetsNauty();
    if(nbNautyVert!=b.getNbSommetsNauty()) { 
        return false;
    }

    // If here, the number of verticies is the same
    // We must ensure that the actual multiplicity map is
    // the same
  
    mul_a = a.getMultiplicitiesMap();
    mul_b = b.getMultiplicitiesMap();
    
    if(mul_a->size() != mul_b->size())
    { 
        return false;
    }
    for(ita=mul_a->begin(),itb=mul_b->begin(),mulend=mul_a->end();ita!=mulend;ita++,itb++)
    {
        if(ita->first  != itb->first)  {
            return false;
        }
        if(ita->second != itb->second) {
            return false;
        }
    }
    
    // Now that we are certain the multiplicities are the same,
    // compare nauty map 
    for(i=0;i<nbNautyVert;i++)
    {
        if(c1[i] != c2[i])
        {
            return false;
        }
    }

    // All possibilities are exhausted
    // The two graphs are Isomorphs !
    return true;
}
