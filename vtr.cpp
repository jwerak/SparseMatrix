#include "vtr.h"
#include <fstream>
#include <cstdio>

// ==============================================================================
                        //KONSTRUKTORY A DESTRUKTORY
// ==============================================================================
/// PRVNI KONSTRUKTOR, ALOKUJE PAMET PRO POCET PRVKU A NAPLNI JE NULAMA
VTR::VTR(long n)
{
    length = n;
    val = new double [length];
    for (long i=0;i<length;i++)
        val[i]=0;
    allocated = true;
}


// ==============================================================================
/**
  \brief Upravi posloupnost sloupcu ve vektoru podle permutace P
  \param[in] P - permutace zmeny
*/
void VTR::apply_permutation(const long *P)
{
    long i;
    double* new_val = new double[length];
    for(i=0;i<length;i++)
        new_val[P[i]]=val[i];
    delete[] val;
    val = new_val;
}


// ==============================================================================
                                //OSTATNI METODY
// ==============================================================================


// ==============================================================================
/**
  \brief Naplni vektor jednickami dle vzoru scilabu
*/
void VTR::ones()
{
    for(long i=0;i<length;i++)
        val[i]=1;
}


// ==============================================================================
/**
  \brief Vytiskne vektor
*/
void VTR::print()const
{
    for (long i=0; i<length; i++)
        printf("%lf\t",val[i]);
    printf("\n");
}


// ==============================================================================
/**
  \brief Nacte vektor ze souboru. Ocekava, ze prvni udaj je pocet prvku
*/
void VTR::load_from_file(string cesta)
{
    fstream vstup;
    vstup.open(cesta.c_str());//otevri soubor
    vstup>>length;//prvni hodnota je pocet hodnot a zaroven delka vektoru
    initialize(length);
    for(long i=0;i<length;i++)
        vstup>>val[i];
}


// ==============================================================================
/**
  \brief Zapise vektor do souboru s nazvem cesta. Na prvni misto zapise pocet prvku
*/
void VTR::print_to_file(string cesta)
{
     FILE* out;
     out = fopen(cesta.c_str(),"w");
     fprintf (out, "%li\n",length);
     for(long i=0;i<length;i++)
         fprintf (out, "%2.20lf\n",val[i]);
     fclose(out);
}


// ==============================================================================
/**
  \brief
*/
void VTR::operator=(const VTR& b)
{
    if (length < b.length){
        dealloc();
        alloc(b.length);
    }
    length = b.length;

    for(long i=0;i<b.length;i++){
        val[i]=b.val[i];
    }
}


// ==============================================================================
//   ===============================================================================
//                                  PRIVATE
//   ===============================================================================


// ==============================================================================
/**
  \brief Alokuje pole.
*/
void VTR::alloc(const long& n)
{    
    if(!allocated){
        val = new double[n];
    }else{
        dealloc();
        alloc(n);
        return;
    }
    length = n;
}


// ==============================================================================
/**
  \brief Naplni vektor nulama
*/
void VTR::fill_zeros()
{
    for(long i=0;i<length;i++)
        val[i]=0;
}


// ==============================================================================
/**
  \brief Nakopiruje do sebe hodnoty o delce n z pole v.
  Nezajima ho jestli je pole do ktereho kopiruje dost dlouhe
*/
void VTR::copy_in(double* v, long n)
 {
     for (long i=0;i<n;i++)
         val[i]=v[i];
 }
