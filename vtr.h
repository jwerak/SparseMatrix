#ifndef VTR_H
#define VTR_H
#include <iostream>
#include "math.h"

using namespace std;

//   ===============================================================================
//   ===============================================================================
//                                  CLASS VTR
//   ===============================================================================
//   ===============================================================================
/**
  Slouzi pro uchovavani plneho radku. Pouziva se kdyz je trida ROW zbytecne komplikovana
  */
class VTR
{
    friend class    SPARSE;
    friend class    MTX;
    friend class    ROW;
    //   ===============================================================================
    //                                  PUBLIC
    //   ===============================================================================
public:
    //    ===================== PUBLIC METHODS ================================
    VTR(double* v, long n):length(n),val(v),allocated(true){}  // VYTVORI VEKTOR Z POLE
    VTR():length(0),allocated(false){}
    VTR(long);
    ~VTR(){if(allocated)
        delete[] val;}

    void    apply_permutation(const long *P);
    long    get_length()const{return length;}// POCET HODNOT VEKTORU

    //    ===================== OSTATNI METODY ================================
    void    initialize(long n){alloc(n); fill_zeros();}
    double  get_element(long misto) const {return val[misto];}
    void    add_element(long misto, double hodnota) {val[misto]+=hodnota;}
    void    set_element(long misto, double hodnota) {val[misto]=hodnota;}
    void    ones();
    void    print()const;
    void    load_from_file(string cesta);
    void    print_to_file(string cesta);

    //    ===================== PUBLIC OPERATORY ================================
    void    operator=(const VTR& b);
//   ===============================================================================
//                                  PRIVATE
//   ===============================================================================
private:
    //    ===================== PRIVATE VARIABLES ================================
    long            length;     //  POCET CLENU VEKTORU
    long            val_number; //POCET NENULOVYCH HODNOT - Bude pouzivat i jako castecne ridky radek
    double*         val;        //  JEDNOTLIVE HODNOTY VEKTORU
    bool            allocated;

    //    ===================== PRIVATE METHODS ================================
    void            alloc(const long& n);
    void            dealloc(){delete[] val; allocated=false;}
    void            fill_zeros();
    void            copy_in(double*, long n);
    long            get_index_of_next(long i){
        /**
          \brief Najde index dalsi hodnoty ve vektoru.
          \warning Nelze pouzit pro posledni prvek,pretekl by. Nebudu osetrovat, nechci pouzivat
          */
        while(val[++i] == 0){}
        return i;
    }
};

#endif // VTR_H
