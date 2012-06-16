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
#include "greenfinder.hpp"

GreenFinder::GreenFinder(PrincipalExtension pe, mpz_class p=0, int min_depth=0, int max_depth=-1)
{
    depart = new PrincipalExtension(pe);
    this->p = p;
    this->min_depth = min_depth;
    this->max_depth = max_depth;
    srand(time(NULL));
}

void GreenFinder::find(uint64_t tries)
{
    int ret;
    int i;
    int sommet;
    int totalTries=tries;
    uint64_t cutMax=0, cutMin=0, cutNoMore=0, cutInf=0;
    PrincipalExtension *pe;
    // Main loop
    while((tries != 0) && (ret !=1)) {
        pe = new PrincipalExtension(*depart);
        tries--;
        std::cout <<  (char)27 << "[2K" << (char)27 << "E" 
                  << "Try " << totalTries-tries << "/" << totalTries 
                  << " (" << round((totalTries-tries)*10000.0/totalTries)/100
                  << "%)" << std::flush;
        while(true) {
            sommet = pe->getRandomGreenVertex();
            if(sommet == -1) {
                cutNoMore++;
                break;
            }
            if(sommet==pe->lastMutation())
            {
                // No need to mutate twice on the same vertex
                continue;
            }
            if(pe->getMutationsSize() >= this->max_depth) {
                cutMax++;
                break;
            }
            ret = pe->mutate(sommet, this->p);
            if (ret == 0) {
                // if mutate returned 0, then infinity was detected
                cutInf++;
                break;
            }
            if (ret == 1) {
                if(pe->getMutationsSize() < this->min_depth) {
                    cutMin++;
                    ret=-1;
                    break;
                }
                // a Green quiver was detected !
                    std::cout << std::endl << "Found !" << std::endl;
                    std::cout << "Size: "; pe->printMutationsE(1);
                    std::cout << "Suite: ";pe->printMutationsE(0);
                break;
            }
        }
        delete pe;
    }
    if(tries == 0)
    {
        std::cout << std::endl << "Not Found !" << std::endl;
    }
    std::cout << "Branches Cut:" << std::endl;
    std::cout << "\t>= max_depth: " << cutMax << std::endl;
    std::cout << "\t<  min_depth: " << cutMin << std::endl;
    std::cout << "\tno more Greens: " << cutNoMore << std::endl;
    std::cout << "\tInfinity: " << cutInf << std::endl;
    std::cout << "Total Cut: " << cutInf + cutMin + cutMax + cutNoMore << std::endl;
}
