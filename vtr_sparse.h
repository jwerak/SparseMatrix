#ifndef VTR_SPARSE_H
#define VTR_SPARSE_H

#include "row.h"

/**
  \brief Kompromis mezi vtr a row. Je to plny radek ve kterem je zapsano -1 dokud dokud na danem miste
    neni zapsana hodnota.
  */
class VTR_SPARSE
{
    friend class MTX;
    friend class SPARSE;
public:
    VTR_SPARSE();
    VTR_SPARSE(const long& n);

    //============================================================
    //====================== PUBLIC METHODS ======================
    //============================================================
    void            set_row(const ROW& r);
    void            get_next_ac(long& r);

    long            get_first_ac(){ ac = start; return ac; }
    void            initialize(const long& l);
    void            set_val_from_sparse(const double* v, const long* ind, const long& vn);
    void            add_val_from_sparse(const double* v, const long* ind, const long& vn, const double& m=1);
    void            print_full();

private:
    //============================================================
    //===================== PRIVATE VARIABLES ====================
    //============================================================
    long            length;//delka radku
    long            value_number; //POCET NENULOVYCH HODNOT - Bude pouzivat i jako castecne ridky radek
    double*         values;        //JEDNOTLIVE HODNOTY VEKTORU
    long*           index;    //Indexy hodnot, ktere jsou vpravo od dane hodnoty
    bool            allocated;  //Informuje o tom, jestli jsou pole alokovana
    long            ac;//Jaky je aktualni sloupec pri rozkladani
    long            start;//Ktery sloupec je prvni v radku//posledni je index[posledni]=-2;

    //============================================================
    //===================== PRIVATE METHODS ======================
    //============================================================.
    /// \brief Najde index pozice, kde je nejblizsi hodnota vlevo od i
    long            find_val_left(long i){
        i--;
        while(index[i] == -1)
            i--;
        return i;
    }

    /// \brief Najde index pozice, kde je nejblizsi hodnota vpravo od i
    long            find_val_right(long i){
        i++;
        while(index[i] == -1 && i < (length-1))
            i++;
        return i;
    }

    void            anul_val();
    void            alloc();
    void            dealloc();
};

#endif // VTR_SPARSE_H
