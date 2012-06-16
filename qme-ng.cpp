/*
 * Copyright (c) 2005-2012, Grégoire Dupont, Matthieu Pérotin
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

#include <boost/program_options.hpp>
#include "greenexplorator.hpp"
#include "greenfinder.hpp"
#include "carquois.hpp"
#include "mutexploratorSeq.hpp"
#include "mutexplorator.hpp"
#include "time.h"
#include <stdint.h>

int main(int argc, char **argv)
{
    int **test;
    mpz_class p;
    int max_depth;
    int min_depth;
    int mutationClassSize;
    bool iso=true;
    int random_tries;
    std::vector<int> size;
    std::string type;
    Carquois *carquois;
    GreenExplorator ge;
    GreenFinder *gf;
    MutExploratorSeq *explorator;
    PrincipalExtension *pt;
    std::map<std::string,int> type_carquois;
    std::map<std::string,int>::const_iterator it_map;
    /* Initializing type_carquois map */
    type_carquois["A"] = A;
    type_carquois["D"] = D;
    type_carquois["E"] = E;
    type_carquois["ATILDE"] = ATILDE;
    type_carquois["ATILDEALT"] = ATILDEALT;
    type_carquois["DTILDE"] = DTILDE;
    type_carquois["ETILDE"] = ETILDE;
    type_carquois["SPORADIQUE"] = SPORADIQUE;
    type_carquois["UNAMED"] = UNAMED;
    type_carquois["E_ELIPTIQUE"] = E_ELIPTIQUE;

    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Prints this help")
        ("file", boost::program_options::value<std::string>(), "Quiver file")
        ("pefile", boost::program_options::value<std::string>(), "Principal Extension file")
        ("type", boost::program_options::value<std::string>(), "Quiver type (A, D, E, ATILDE, DTILDE, ETILDE, SPORADIQUE, UNAMED, E_ELIPTIQUE)")
        ("size", boost::program_options::value< std::vector<int> >(),"Quiver size (must be used with type)")
        ("green", "Green exploration")
        ("one", boost::program_options::value<int>(&random_tries)->default_value(0), "Find one green suite, give number of tries")
        ("p", boost::program_options::value<mpz_class>(&p)->default_value(0),"P param")
        ("max_depth", boost::program_options::value<int>(&max_depth)->default_value(INT_MAX),"Max exploration depth")
        ("min_depth", boost::program_options::value<int>(&min_depth)->default_value(0),"Min exploration depth")
        ("no-iso", "Isomorph discrimination")
        ("dump-class", boost::program_options::value<std::string>(), "Dump Mutation Class")
        ("dump-trunk", "Dump truncated quivers")
    ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    ge.setP(p);
    ge.setMaxDepth(max_depth);
    if(vm.count("no-iso"))
    {
        ge.setIsomorphTest(false);
    }
    if(vm.count("dump-trunk"))
    {
        ge.setDumpTruncated(true);
    }
    try
    {
        if(vm.count("file") && vm.count("pefile")) {throw Exception("Only one type of file can be given !");}
        if(vm.count("file"))
        {
                    carquois = new Carquois(vm["file"].as<std::string>().c_str());
        }
        else if(vm.count("type"))
        {
            if(vm.count("size"))
            {
                size = vm["size"].as< std::vector<int> >();
                type = vm["type"].as<std::string>();
                it_map = type_carquois.find(type);
                if(it_map == type_carquois.end())
                {
                    // Unknown type
                    std::cerr << desc << "\n";
                    return 1;
                }
                else
                {
                    if(size.size() == 1) {
                        carquois = new Carquois(it_map->second,size[0]);
                    }
                    if(size.size() == 2) {
                        carquois = new Carquois(it_map->second,size[0],size[1]);
                    }
                }    
            }
            else
            {
                std::cerr << desc << "\n";
                return 1;
            }
        }
        else
        {
            if(!vm.count("pefile")) {
                std::cerr << desc << "\n";
                return 1;
            }
        }
    } catch (Exception e)
    {
        std::cout << e.m_Msg << "\n";
        return 1;
    }

    if(vm.count("green"))
    {
        if(vm.count("pefile")) { pt = new PrincipalExtension(vm["pefile"].as<std::string>().c_str());}
        else  {
            pt = new PrincipalExtension(*carquois);
        }
        pt->affiche();
        pt->generateSommetsVerts();
        try {
            if(random_tries == 0) {
                ge.greenExploration(*pt);
            }
            else
            {
                gf = new GreenFinder(*pt, p, min_depth, max_depth);
                gf->find(random_tries);
                delete gf;
            }
        } catch (Exception e)
        {
            std::cout << e.m_Msg << "\n";
            delete carquois;
            return 1;
        }
    }
    else
    {
        explorator = new MutExploratorSeq();
        try {
        mutationClassSize = explorator->isomorphismExplorator(*carquois,5000);
        std::cout << "Set Size:" << mutationClassSize << "\n";
        explorator->getNbVoisinsMax();
        if(vm.count("dump-class")) {
            explorator->dumpFiles((vm["dump-class"].as<std::string>()).c_str());
        }
        if(explorator->acyclique()) {
            std::cout << "The mutation class is acyclic !\n";
        }
        else
        {
            std::cout << "The mutation class is not acyclic !\n";
        }
        } catch (Exception e)
        {

            std::cout << e.m_Msg << "\n";
            delete explorator;
            return 1;
        }
        delete explorator;
    }
    if(!vm.count("pefile")) {
        delete carquois;
    }
    delete pt;
    return 0;
}
