#include <fstream>
#include "vtr_sparse.h"

VTR_SPARSE::VTR_SPARSE()
{
    allocated = false;
    value_number = 0;
}


// ==============================================================================
/**
  \brief Konstruktor
  \param n - delka vytvareneho vektoru
  */
VTR_SPARSE::VTR_SPARSE(const long& n){
    allocated = false;
    initialize(n);
}


// ==============================================================================
/**
  \brief Ziska dalsi hodnotu v poradi
  \param r - Vystupni hodnota
  */
void VTR_SPARSE::get_next_ac(long& r)
{
    ac = index[ac];
    r = ac;
}


// ==============================================================================
/**
  \brief Zajisti vektor o delce l - melo by se pouzit jen jednou
  \param l - pozadovana delka vektoru
  */
void VTR_SPARSE::initialize(const long &l)
{
    long i;
    if(!allocated){
        length = l;
        alloc();
    }else if(length<l){
        dealloc();
        length = l;
        alloc();
    }
    for(i=0; i<length; i++){//Nastav vychozi hodnoty
        values[i] = 0;
        index[i] = -1;
    }
    start = -2;
}


// ==============================================================================
/**
  \brief Pretvori rikdou strukturu z v a ind do plne struktury pouzite ve VTR_SPARSE
  \param v - Ridke pole hodnot
  \param ind - Ridky zapis sloupcu kde jsou zapsany hodnoty
  \param vn - Value number - kolik nenulovych hodnot je treba zapsat
  \param l - length - plna delka radku
*/
void VTR_SPARSE::set_val_from_sparse(const double* v, const long* ind, const long& vn)
{
    long i, n;
    //PRIPRAV OBJEKT PRO ZAPIS
    anul_val();

    //ZAPIS HODNOTY DO POLI
    start = ind[0];

    value_number = vn;

    for(i=0;i<vn;i++){
        n = ind[i];
        values[n] = v[i];
        index[n] = ind[i+1];
    }
    index[n] = -2;// Posledni index oznamuje, ze je posledni
}


// ==============================================================================
/**
  \brief Pricte hodnoty z ridkeho zapisu do VTR_SPARSE.
  \param v - Ridke pole hodnot
  \param ind - Ridky zapis sloupcu kde jsou zapsany hodnoty
  \param vn - Value number - kolik nenulovych hodnot je treba zapsat
  \param m - multiply - Volitelny parametr (vych. hod. = 1), roznasobi se jim hodnoty, ktere se pricitaji
*/
void VTR_SPARSE::add_val_from_sparse(const double* v, const long* ind, const long& vn, const double& m)
{
    long i,l,n;
    for (i=0; i<vn; i++){
        n = ind[i];
        if(index[n] != -1)
            values[n] += v[i]*m;
        else{//Hodnota jeste nebyla zadana
            l = find_val_left(n);
            index[n] = index[l];
            index[l] = n;
            values[n] = v[i]*m;
            value_number++;
        }
    }
}


// ==============================================================================
/**
  \brief Vytiskni pomoci printf.
*/
void VTR_SPARSE::print_full()
{
    for(long i=0; i<length; i++){

        if(index[i] == -1)
        {
            printf("0 \t");
        }else{
            printf("%1.6f \t",values[i]);
        }
    }
    printf("\n");
}


// ==============================================================================
/**
  \brief Nastavi vektor na zacatek - vsechny hodnoty val = 0 a indexy = -1
  \warning Zde pada pri pruchodu z prikazove radky
*/
void VTR_SPARSE::anul_val()
{
    long n = start;
    long j;
    while(n != -2){
        j = index[n];
        values[n] = 0;
        index[n] = -1;
        n=j;
    }
    value_number = 0;
}


// ==============================================================================
/**
  \brief Alokuje potrebna pole
*/
void VTR_SPARSE::alloc()
{
    values      = new double[length];
    index       = new long [length];
    allocated   = true;
}


// ==============================================================================
/**
  \brief Dealokuje potrebna pole
*/
void VTR_SPARSE::dealloc()
{
    delete values;
    delete index;
    allocated = false;
}

