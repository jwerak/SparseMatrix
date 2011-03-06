#ifndef MTX_H
#define MTX_H

#include <stdio.h>
#include <iostream>
#include <cstring>

#include "row.h"
#include "vtr_sparse.h"

#define RECURSIVE_PERMUTATION
class SPARSE;

class MTX
{
    //   ===============================================================================
    //                                  PUBLIC
    //   ===============================================================================
public:
    //    ===================== PUBLIC CONSTRUCTORS ================================

    MTX(): M_aloc( false), addit_values(false), addit_double_number(0), addit_long_number(0), addit_long_bool(false),  addit_double_bool(false), pa(false), sa(false){}
    ~MTX(){if(pa){delete P; delete S;}
        dealloc_addit();
        dealloc_addit_d();
        dealloc_addit_l();
    }
    //    ===================== PUBLIC METHODS ================================
    // ==============================================================================
    void                    solve_Axb(MTX&, VTR& x, VTR& b);
    void                    make_permutation_ND();
    virtual void            make_symbolic_factorization()=0;

    void                    make_ldu(MTX& ldu);// vytvori ldu rozklad z aktualni matice
    virtual void            loadFromFile(string cesta) = 0;
    virtual void            print_full_mtx() = 0;
    virtual void            write_log(double t)=0;


    virtual void            get_mtx_L(MTX& A);
    virtual void            get_mtx_D(MTX& A);
    virtual void            get_mtx_U(MTX& A);

    const long&             get_column_number()const; ///< Vrati pocet sloupcu
    const long&             get_row_number()const; ///< Vrati pocet radku
    virtual const long&     get_values_number()const = 0; ///< Vrati pocet zapsanych prvku (nenulovych)
    virtual ROW             get_row(const long r, const long c=0)const = 0;
    virtual void            set_vtr_sparse(const long r, VTR_SPARSE& vtr)const = 0;

    virtual void            set_row(const ROW& a, const long& r) = 0;
    virtual void            set_row(VTR_SPARSE& a, const long& r) = 0;
    virtual void            set_row(const double* v,const long* i,const long& l, const long& r) = 0;
    virtual double          get_diagonal(const long& r)const = 0;///< Vrati diagonalni prvek na r-tem radku


    //    ===================== PUBLIC OPERATORS ================================
    void                    multiply_matrices(MTX& B, MTX& R);
    virtual void            operator=(const MTX& B) = 0;
    const VTR               multiply(const VTR& b);

    void            reverse_permutation_vector();
    //   ===============================================================================
    //                                  PRIVATE
    //   ===============================================================================
private:
    //    ===================== PRIVATE METHODS ================================

#ifdef RECURSIVE_PERMUTATION
    virtual void    find_sep(long p, long& ap, long* S, long** fl, long** fr, long* TMP)=0;
#else
    virtual void    find_sep(long p, long* n, long& in, long& ap, long* S, long** fl, long** fr)=0;
#endif

    virtual void    apply_permutation()=0;

    virtual void    initialize(const long& nor, const long& noc, const long& nov, const int& rm=0) = 0;
    virtual void    update_ldu_row(long ac, long j, ROW& ldu) = 0;
    virtual void    update_ldu_row(long ac, VTR_SPARSE& ldu) = 0;
    void            solve_Lzb(const VTR& b, VTR& z);
    void            solve_Dyz(const VTR& z, VTR& y);
    void            solve_Uxy(const VTR& y, VTR& x);
    virtual ROW     get_rowL(const long& r)const=0;
    virtual ROW     get_rowU(const long& r)const=0;
    //ZBYTECNE JEN TESTOVACI:
    virtual void    print_sparse_row(const long& N)=0;

    //    ===================== PRIVATE VARIABLES ================================
    //   ===============================================================================
    //                                  PROTECTED
    //   ===============================================================================
protected:
    //    ===================== PROTECTED VARIABLES ================================
    const static int    addit_fields = 2; ///< pocet doplnkovych poli o delce jednoho radku
    long        row_number, column_number, value_number; ///< pocet danych hodnot(row, column se nesmi menit)
    bool        M_aloc; ///< Boolean, ktery sdeluje jestli jsou maticovy pole alokovany, nebo ne.
    bool        addit_values;///< Informuje o tom jestli jsou pole addit_val a addit_col alokovana
    double*     addit_val[addit_fields];///< pole hodnot o delce radku 1, pomocne pole pro ruzne vypocty
    long*       addit_col[addit_fields];///< pole hodnot o delce radku 1, pomocne pole pro ruzne vypocty
    long        addit_double_number;///< Kolik poli je alokovano
    long        addit_long_number;///< Kolik poli je alokovano
    long**      addit_long;///< Pole longu, pocet poli bude promenny
    double**    addit_double;///< Pole longu, pocet poli bude promenny
    bool        addit_long_bool;///< je alokovan?
    bool        addit_double_bool;///< je alokovan?
    long*       P;///< permutation vector
    long*       S;///< Helping vector for permutation
    bool        pa;///< permutation allocated
    bool        sa;///< Field S allocated?

    //    ===================== PROMENNE PRO SPOCTENI OPERACI ================================
    long        value_number_orig;
    long        na,ns,nm,nd;//Number of - addition, subtraction, multipl., division.


    //    ===================== PROTECTED METHODS ================================
    void    BB_sortl(const long n, long* a, long* k); //setridi dane pole longu
    void    BB_sortd(const long n, double* a); //setridi dane pole doublu
    void    alloc_addit();
    void    alloc_addit_l(long n);
    void    alloc_addit_d(long n);
    void    dealloc_addit();
    void    dealloc_addit_l();
    void    dealloc_addit_d();

    void    update_field_L_by_key(const long n, long* a, long* k);
    void    update_field_D_by_key(const long n, double* a, long* k);
};

#endif // MTX_H

