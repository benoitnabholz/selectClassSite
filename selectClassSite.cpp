/**
g++ -std=c++14 --static -g ~/Dropbox/Progr/C++/seq/selectClassSite.cpp -o ~/bin/selectClassSite \
 -DVIRTUAL_COV=yes -Wall \
 -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib \
 -lbpp-popgen -lbpp-seq -lbpp-core 

 
g++ --static -g ~/pCloudDrive/Progr/C++/seq/selectClassSite.cpp -o ~/bin/selectClassSite \
 -DVIRTUAL_COV=yes -Wall \
 -I$HOME/local/bpp/dev/include  -L$HOME/local/bpp/dev/lib \
 -lbpp-popgen -lbpp-seq -lbpp-core 

strip ~/bin/selectClassSite
  
**/

//
// Created by: Benoit Nabholz
//

/*

 
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for sequences analysis.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently readAlignmenting this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/



#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/Seq/GeneticCode/StandardGeneticCode.h>
#include <Bpp/Seq/GeneticCode/VertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/InvertebrateMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/EchinodermMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/AscidianMitochondrialGeneticCode.h>
#include <Bpp/Seq/GeneticCode/MoldMitochondrialGeneticCode.h>

#include <Bpp/PopGen/SequenceStatistics.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>

#include <Bpp/Numeric/VectorTools.h>

#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <iostream>



using namespace std;
using namespace bpp;


/************************************************/

bool isTransversion(int i, int j){

  if(!(((i==0 && j==2) || (i==2 && j==0)) ||
           ((i==1 && j==3) || (i==3 && j==1)))) {
    return true;
  }else{
    return false;
  }

}
/************************************************/
Site sansGap(Site site){

  Site sansGapSite = Site(site.getAlphabet());

  for (unsigned int j = 0; j < site.size(); j++){
    // if not gap count
    if (!(site.getAlphabet()->isGap(site.getValue(j))) && !(site.getAlphabet()->isUnresolved(site.getValue(j)))){
      sansGapSite.addElement(site.getValue(j));
    }
  }
  return sansGapSite;
}



/******************************************************************************/
/******************************************************/
unsigned int L4(int i, const GeneticCode & gc){

  const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
  if(gc.isStop(i)) return 0;

  unsigned int nbs4FoldSite = 0;
  int nbsynmut = 0.0;

  vector<int> codon = ca->getPositions(i);
  int acid = gc.translate(i);

  
  // test all the substitution on third codon position
  for (int an=0; an < 4; an++) {
    if (an == codon[2]) continue;
    vector<int> mutcodon = codon;
    mutcodon[2] = an;
    int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
    if(gc.isStop(intcodon)) continue;
    int altacid = gc.translate(intcodon);
    if (altacid == acid) { //if synonymous
      nbsynmut ++;
    }
  }

  if(nbsynmut == 3){
     nbs4FoldSite = 1;
  }
  return nbs4FoldSite;
}
/******************************************************************************/
bool isFourFoldDegenerated(const Site& site, const GeneticCode& gc)
{
  if (!SiteTools::isConstant(site, true))
  {
    /** If non-synonymous mutation **/
    if (!(CodonSiteTools::isSynonymousPolymorphic(site, gc)))
      return false;

    for (unsigned int i = 0; i < site.size(); i++)
    {
      if (!(gc.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  else
  {
    for (unsigned int i = 0; i < site.size(); i++)
    {
      if (!(gc.isFourFoldDegenerated(site.getValue(i))))
      {
        return false;
      }
    }
  }
  return true;
}

/******************************************************************************/
double L2v (int i, const GeneticCode & gc) throw(Exception)
{
  try {
    const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
    if(ca->getName(ca->intToChar(i))=="Stop") return 0;

    vector<int> codon = ca->getPositions(i);
    int acid = gc.translate(i);

    double nbsyn2V = 0.0;

    for (int pos=0; pos < 3; pos++) {

      double nbsynS = 0;
      double nbsynV = 0;
      double nbnonsynS = 0;
      double nbnonsynV = 0;

      for (int an=0; an < 4; an++) {

        if (an == codon[pos]) continue;

        vector<int> mutcodon = codon;
        mutcodon[pos] = an;
        int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);

        if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;

        int altacid = gc.translate(intcodon);
        if (altacid == acid) { //if synonymous
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbsynV ++;
          } else { //if transition
            nbsynS ++;

          }
        }else{
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbnonsynV ++;
          } else { //if transition
            nbnonsynS ++;

          }
        }
      }


      // if nbsynS == 2 the site is 4 fold degenerated
      if(nbsynV > 0 && (nbsynV+nbsynS < 3)){
        // special case of arginine and isoleucine fold degenerated site
        if(nbsynV == 1 && (nbsynV+nbsynS+nbnonsynV+nbnonsynS) == 3){
          nbsyn2V += nbsynV/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
        }else{
          nbsyn2V += 1;
        }
      }

    }
//     cout << nbsyn2V<<endl;
    return nbsyn2V;

  } catch (...) {} // !!!!! en cas d'exception, plante! il faudrait forwarder l'exception
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/
double posL2v (int i, const GeneticCode & gc, int position) throw(Exception)
{
  try {
    const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
    if(ca->getName(ca->intToChar(i))=="Stop") return 0;

    vector<int> codon = ca->getPositions(i);
    int acid = gc.translate(i);

    double nbsyn2V = 0.0;

    int pos = position;

    double nbsynS = 0;
    double nbsynV = 0;
    double nbnonsynS = 0;
    double nbnonsynV = 0;

    for (int an=0; an < 4; an++) {

      if (an == codon[pos]) continue;

      vector<int> mutcodon = codon;
      mutcodon[pos] = an;
      int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);

      if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;

      int altacid = gc.translate(intcodon);
      if (altacid == acid) { //if synonymous
        if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
          nbsynV ++;
        } else { //if transition
          nbsynS ++;

        }
      }else{
        if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
          nbnonsynV ++;
        } else { //if transition
          nbnonsynS ++;
        }
      }
    }


    // if nbsynS == 2 the site is 4 fold degenerated
    if(nbsynV > 0 && (nbsynV+nbsynS < 3)){
      // special case of arginine and isoleucine fold degenerated site
      if(nbsynV == 1 && (nbsynV+nbsynS+nbnonsynV+nbnonsynS) == 3){
        nbsyn2V += nbsynV/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
      }else{
        nbsyn2V += 1;
      }
    }
//     cout << nbsyn2V<<endl;
    return nbsyn2V;

  } catch (...) {} // !!!!! en cas d'exception, plante! il faudrait forwarder l'exception
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/

double L2s(int i, const GeneticCode & gc) throw(Exception)
{
  try {
    const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
    if(ca->getName(ca->intToChar(i))=="Stop") return 0;

    vector<int> codon = ca->getPositions(i);
    int acid = gc.translate(i);

    double nbsyn2S = 0;


    for (int pos=0; pos < 3; pos++) {

      unsigned int nbsynS = 0;
      unsigned int nbsynV = 0;
      double nbnonsynS = 0;
      double nbnonsynV = 0;

      for (int an=0; an < 4; an++) {

        if (an == codon[pos]) continue;

        vector<int> mutcodon = codon;
        mutcodon[pos] = an;
        int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);

        if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;

        int altacid = gc.translate(intcodon);
        if (altacid == acid) { //if synonymous
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbsynV ++;
          } else { //if transition
            nbsynS ++;

          }
        }else{
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbnonsynV ++;
          } else { //if transition
            nbnonsynS ++;

          }
        }
      }

      // if nbsynV == 2 the site is 4 fold degenerated
      if(nbsynS > 0 && (nbsynV+nbsynS < 3) ){
        // special case of arginine and isoleucine fold degenerated site
        if(nbsynV+nbsynS == 2 && nbsynV == 1){
          nbsyn2S += nbsynS/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
        }else{
          nbsyn2S += 1;
        }
      }

    }

    return nbsyn2S;

  } catch (...) {} // !!!!! en cas d'exception, plante! il faudrait forwarder l'exception
  // This line is never reached but sends a warning if not there:
  return 0.;
}

/******************************************************************************/
double posL2s(int i, const GeneticCode & gc, int position) throw(Exception)
{
  try {
    const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
    if(ca->getName(ca->intToChar(i))=="Stop") return 0;

    vector<int> codon = ca->getPositions(i);
    int acid = gc.translate(i);

    double nbsyn2S = 0;

    int pos = position;

    unsigned int nbsynS = 0;
    unsigned int nbsynV = 0;
    double nbnonsynS = 0;
    double nbnonsynV = 0;

    for (int an=0; an < 4; an++) {

      if (an == codon[pos]) continue;

      vector<int> mutcodon = codon;
      mutcodon[pos] = an;
      int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);

      if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;
      int altacid = gc.translate(intcodon);
      if (altacid == acid) { //if synonymous
        if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
         nbsynV ++;
        } else { //if transition
          nbsynS ++;
        }
      }else{
        if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
          nbnonsynV ++;
        } else { //if transition
          nbnonsynS ++;

        }
      }
     }

    // if nbsynV == 2 the site is 4 fold degenerated
    if(nbsynS > 0 && (nbsynV+nbsynS < 3) ){
      // special case of arginine and isoleucine fold degenerated site
      if(nbsynV+nbsynS == 2 && nbsynV == 1){
        nbsyn2S += nbsynS/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
      }else{
        nbsyn2S += 1;
      }
    }

    return nbsyn2S;

  } catch (...) {} // !!!!! en cas d'exception, plante! il faudrait forwarder l'exception
  // This line is never reached but sends a warning if not there:
  return 0.;
}
/******************************************************************************/
double L0(int i, const GeneticCode & gc) throw(Exception)
{
  const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
  if(ca->getName(ca->intToChar(i))=="Stop") return 0;

  double nbs0FoldSite = 0;

  vector<int> codon = ca->getPositions(i);
  int acid = gc.translate(i);

    for (int pos=0; pos < 3; pos++) {

      unsigned int nbsynS = 0;
      unsigned int nbsynV = 0;
      double nbnonsynS = 0;
      double nbnonsynV = 0;

      for (int an=0; an < 4; an++) {
        if (an == codon[pos]) continue;
        vector<int> mutcodon = codon;
        mutcodon[pos] = an;
        int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
        if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;
        int altacid = gc.translate(intcodon);
        if (altacid == acid) { //if synonymous
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbsynV ++;
          } else { //if transition
            nbsynS ++;

          }
        }else{
          if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
            nbnonsynV ++;
          } else { //if transition
            nbnonsynS ++;

          }
        }
      }

      if(nbnonsynS > 0 || nbnonsynV > 0 ){
        // special case of arginine and isoleucine fold degenerated site
        if(nbsynV == 1 && (nbsynV+nbsynS+nbnonsynV+nbnonsynS) == 3){
          nbs0FoldSite += (nbnonsynS+nbnonsynV)/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
        }
        if(nbsynV+nbsynS == 0){
          nbs0FoldSite += 1;
        }
      }
    }

  return nbs0FoldSite;

}
/******************************************************************************/
double posL0(int i, const GeneticCode & gc, int position) throw(Exception)
{
  const CodonAlphabet * ca = dynamic_cast<const CodonAlphabet *>(gc.getSourceAlphabet());
  if(ca->getName(ca->intToChar(i))=="Stop") return 0;

  double nbs0FoldSite = 0;

  vector<int> codon = ca->getPositions(i);
  int acid = gc.translate(i);

  int pos = position;

  unsigned int nbsynS = 0;
  unsigned int nbsynV = 0;
  double nbnonsynS = 0;
  double nbnonsynV = 0;

  for (int an=0; an < 4; an++) {
    if (an == codon[pos]) continue;
    vector<int> mutcodon = codon;
    mutcodon[pos] = an;
    int intcodon = ca->getCodon(mutcodon[0], mutcodon[1], mutcodon[2]);
    if(ca->getName(ca->intToChar(intcodon))=="Stop") continue;
    int altacid = gc.translate(intcodon);
    if (altacid == acid) { //if synonymous
      if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
        nbsynV ++;
      } else { //if transition
        nbsynS ++;
      }
    }else{
      if(isTransversion(codon[pos],mutcodon[pos])) { // if it is a transversion
        nbnonsynV ++;
      } else { //if transition
        nbnonsynS ++;

      }
    }
  }

  if(nbnonsynS > 0 || nbnonsynV > 0 ){
    // special case of arginine and isoleucine fold degenerated site
    if(nbsynV == 1 && (nbsynV+nbsynS+nbnonsynV+nbnonsynS) == 3){
      nbs0FoldSite += (nbnonsynS+nbnonsynV)/(nbsynV+nbsynS+nbnonsynV+nbnonsynS);
    }
    if(nbsynV+nbsynS == 0){
      nbs0FoldSite += 1;
    }
  }

  return nbs0FoldSite;
}

  /***************         MAIN             *****************************/

int main (int argc[], char *argv[]){
try{


if(argv[1] == NULL){
  cout << "\n##################\nselectClassSite file_name fomat 'Standard or VertebrateMitochondrial or InvertebrateMitochondrial' 'L4 or third or independent_codon' " << endl;
  cout << "\n\tformat : fasta or phylip" << endl;
  cout << "\tthird: Extract third codon position" << endl;
  cout << "\tL4: Extract fourfold degenrate site" << endl;
  cout << "\tindependent_codon: print two files with independent codon (selected randomly)" << endl;
  cout << endl;
  cout << "####################\nPlease cite Bio++ ( http://biopp.univ-montp2.fr/ ;  Guéguen et al. 2013 MBE) if you use this program )\n" << endl;
  return 0;
}

string nameSeq = argv[1];
string format = argv[2];
string type = argv[3];
string typeMut = argv[4];

// StringTokenizer::StringTokenizer st =  StringTokenizer(nameSeq, "./");

cout << nameSeq << endl;

string nOut = "";

if(typeMut == "L4")
  nOut = nameSeq+".L4.fst";

if(typeMut == "third")
  nOut = nameSeq+".3.fst";

Fasta Fst(10000);
Phylip Ph(true,true);


const NucleicAlphabet *alpha = new DNA();
const GeneticCode *geneticcode = NULL;
const CodonAlphabet *CodonAlpha = new CodonAlphabet(alpha);
const ProteicAlphabet Prot = ProteicAlphabet();



if(type == "Standard"){
  geneticcode = new StandardGeneticCode(alpha);
}
if(type == "VertebrateMitochondrial"){
  geneticcode = new VertebrateMitochondrialGeneticCode(alpha);
}

if(type == "InvertebrateMitochondrial"){
  geneticcode = new InvertebrateMitochondrialGeneticCode(alpha);
}

if(type != "Standard" && type != "VertebrateMitochondrial" && type != "InvertebrateMitochondrial"){
  cout << "type should be Standard or VertebrateMitochondrial or InvertebrateMitochondrial!!" << endl;
  return 0;
}

SequenceContainer * tmp = NULL;
if(format == "fasta"){
 tmp = Fst.readAlignment(nameSeq , alpha);
}
if(format == "phylip"){
 tmp = Ph.readAlignment(nameSeq , alpha);
}
OrderedSequenceContainer * tmp1 = new VectorSequenceContainer(*tmp);
VectorSiteContainer * sequences = new VectorSiteContainer(*tmp1);
delete tmp;

// convert alphabet
SequenceContainerTools st;
VectorSiteContainer * sc1 = new VectorSiteContainer(CodonAlpha);
st.convertAlphabet(*sequences, *sc1);

VectorSiteContainer * scOut = new VectorSiteContainer(sequences->getSequencesNames(),CodonAlpha);

/** start at 1 because remove ATG **/
unsigned int start, end;
start = end = 0;

const Site sitef = sc1->getSite(0);
const Site sitel = sc1->getSite(sc1->getNumberOfSites()-1);
Site sgf = sansGap(sitef);
Site sgl = sansGap(sitel);

if(sgf.size() > 0){
  if(Prot.intToChar(geneticcode->translate(sgf.getValue(0))) == "M"){
    start = 1;
  }
}

if(sgl.size() > 0){
  for(unsigned int i = 0; i < sgl.size(); i++){
    if(geneticcode->isStop(sgl.getValue(i))){
      end = sc1->getNumberOfSites()-1;
      break;
    }else{
      end = sc1->getNumberOfSites();
    }
  }
}else{
  end = sc1->getNumberOfSites();
}

ofstream actualPosition;
actualPosition.open ("Actual_Position_L4.txt");
for(unsigned int i = start; i < end; i++){

  //ApplicationTools::displayGauge(i, end, '=', "select classe sites");

  const Site site = sc1->getSite(i);
  Site SiteSg = sansGap(site);


  if(SiteSg.size() == 0){
// 		cout << (start*3)+(i*3) << endl;
    continue;
  }

  /** test if site has non-synonymous change **/
  if(!(SiteTools::isConstant(SiteSg,true))){
    if(!(CodonSiteTools::isSynonymousPolymorphic(SiteSg,*geneticcode)) && typeMut == "L4"){
      continue;
    }
  }
  
  if(typeMut == "L4"){
    if(L4(SiteSg.getValue(0), *geneticcode) > 0){
	  actualPosition << i << endl;
      scOut->addSite(site);
    }
  }
  
  if(typeMut == "third"){
      scOut->addSite(site);
  }

}

VectorSiteContainer * sc2 = new VectorSiteContainer(alpha);
st.convertAlphabet(*scOut, *sc2);
VectorSiteContainer * scOut3pos = new VectorSiteContainer(sequences->getSequencesNames(),alpha);


/** take third codon position **/
if(typeMut == "third" || typeMut == "L4"){
	
  int compt = 0;
  for(unsigned int i = 0; i < sc2->getNumberOfSites(); i++){
    ApplicationTools::displayGauge(i, sc2->getNumberOfSites(), '=', "creat output");
    compt ++;
    if(compt == 3){
      compt = 0;
      const Site site = sc2->getSite(i);
      scOut3pos->addSite(site);
    }
  }
  cout << endl << "Number of sites: " << scOut3pos->getNumberOfSites() << endl;
  Fst.writeAlignment(nOut,*scOut3pos);
}


/** print two file with independent codon **/
if(typeMut == "independent_codon"){
  VectorSiteContainer * scOut1 = new VectorSiteContainer(sequences->getSequencesNames(),CodonAlpha);
  VectorSiteContainer * scOut2 = new VectorSiteContainer(sequences->getSequencesNames(),CodonAlpha);

  for(unsigned int i = 0; i < sc1->getNumberOfSites(); i++){
    const Site site = sc1->getSite(i);
    if(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0) > 0.5){
      scOut1->addSite(site);
    }else{
      scOut2->addSite(site);
    }
  }
  nOut =  nameSeq+".f1.fst";
  cout << scOut1->getNumberOfSites() << endl;
  Fst.writeAlignment(nOut,*scOut1);
  nOut =  nameSeq+".f2.fst";
  cout << scOut2->getNumberOfSites() << endl;
  Fst.writeAlignment(nOut,*scOut2);

  delete scOut1;
  delete scOut2;
}


delete alpha;
delete CodonAlpha;
delete geneticcode;
delete scOut;

}
catch(exception & e){
  cout << e.what() << endl;
}


return 0;
}


