#ifndef SPARSE_H
#define SPARSE_H

#include "mtx.h"
#include <iostream>

#define RECURSIVE_PERMUTATION
class SPARSE : public MTX
{
    //   ===============================================================================
    //                                  PUBLIC
    //   ===============================================================================
public:
    SPARSE():value_length_max(0){M_aloc = false; addit_values = false;}
    SPARSE(const SPARSE& A){copy_values(A.get_row_number(), A.get_column_number(), A.get_values_number(), A.value_length, A.value_start, A.column_start, A.diagonal);}
    ~SPARSE (){dealloc();}

    //============================================================
    //===================== Virtualni metody =====================
    //============================================================
    void            operator=(const MTX& B);
    ROW             get_row(const long r, const long c=0)const;
    void            set_vtr_sparse(const long r, VTR_SPARSE& vtr)const;
    void            set_row(const ROW& a, const long& r);
    void            set_row(VTR_SPARSE& a, const long& r);
    void            set_row(const double* v,const long* i,const long& l, const long& r);
    double          get_diagonal(const long& r)const{return *(diagonal[r]);}
    void            loadFromFile(string cesta);
    void            print_full_mtx();
    void            write_log(double);
    const long&     get_values_number()const{return value_number;} ///< Vrati pocet zapsanych prvku (nenulovych)

    //    ===================== PUBLIC OPERATORS ================================
    //**************************************************************
    //********************** Virtualni metody **********************
    //**************************************************************

    //   ===========================================================
    //                                  PRIVATE
    //   ===========================================================
private:
    //============================================================
    //================== PRIVATE METHODS =======================
    //============================================================
    void    alloc(const long& r, const long& nov, const int& rm=0); //alokuj pouzitou matici
    void    dealloc();

    void    copy_values(const long& nor, const long& noc, const long& nov, const long* v_length, const double* const* v_start, const long* const* c_start, const double* const* d);
    void    copy_row_in(const double* v, const long* i,  const long& r, const long& l);
    long    enlarge_field(const long& ph, const long& rm){return (ph + ph*(rm/100.));} //ph-puvodni pocet hodnot rm - rozsireni[%]
    void    order_values_by_row(long* c, long* r, double* v, long** rs, const long rn);
    long    find_edges(long** f, long* ff, long iv, long* S, long* rv, long& nw);


    //============================================================
    //================= PERMUTACE A SPOL. ========================
    //============================================================
    long    make_wave(long* S, long& ap, const long* IW, long* OW, const long& NOI, const long& ST, const long& SF, const long& CT);
    long    make_symbolic_wave(const long& LV, const long& AR, const long& IC, const long& SNV, const long& OL, long* S);


    //============================================================
    //===================== SPRAVA PAMETI ========================
    //============================================================
    //Metody pro spravu zapsane matice
    double* M_get_free_place(const long& l);// najde nove misto v poli value, jestli ne, tak ho doalokuje...
    double* M_find_free_place(const long& l, bool& found); //najde nove misto v poli value pro dalsi radek napriklad.
    void    M_update_empty_pointers(const long& l){M_empty_start[allocated_field]  = (M_empty_start[allocated_field]+l);}//prehazi pointery starajici se o volne misto
    bool    M_alloc_new_field();//zalokuj nove pole v poradi
    void    M_shake_down(const int& rm);//vsechny radky nasype poporade do jednoho pole
    void    initialize_value_length(const long* length);//nastavi pointery na zacatky radku
    
    //============================================================
    //================== PRIVATE VARIABLES =======================
    //============================================================
    // PRACE S MATICI
    const static int    extra_values=10;///< pocet poli, ktera muzou byt doalokovana pri preteceni
    const static int    extra_field_enlargement=50;///< Delka rozsireni poli value[>0], column[>0] v %
    long                extra_field_length;///< Delka poli value[>0], column[>0]
    long*               column[extra_values];///< pole poli hodnot
    double*             value[extra_values];///< pole hodnot, muze byt doalokovana nova pole s pameti
    int                 allocated_field; ///< kolikate pole (*column[], *value[]) se pouziva

    // SPRAVA PAMETI
    long            M_allocated_field_length[extra_values]; ///< rozsiruje pole sloupcu a hodnot=>informuje o delce poli
    double*         M_empty_start[extra_values];///< Ukazatel na **value. Ukazuje na prvni volne misto v kazdem poli
    double*         M_empty_end[extra_values];///< Ukazatel na **value. Ukazuje na posledni volne misto v kazdem poli

    double**        value_start;///< pole pointeru ukazujicich na value, znaci zac radku
    long**          column_start; ///< pole pointeru na cislo sloupce na zacatku radku, obdoba value_start, delka se precte z value_length
    double**        diagonal; ///< pole pointeru na diagonalni hodnoty value. Kdyz se objevi 0 na diagonale, bude obsazena take v poli value \see{value}
    long*           value_length; ///< pocet hodnot v jednotlivych radcich. 1 prvek = 1
    long            value_length_max; ///< maximalni pocet hodnot v jednotlivych radcich. 1 prvek = 1

    //============================================================
    //===================== Virtualni metody =====================
    //============================================================
#ifdef RECURSIVE_PERMUTATION

    void            find_sep(long p, long& ap, long* S, long** fl, long** fr, long* TMP);
#else
    void            find_sep(long p, long* n, long& in, long& ap, long* S, long** fl, long** fr);
#endif
    void            apply_permutation();
    void            make_symbolic_factorization();
    void            initialize(const long& nor, const long& noc, const long& nov, const int& rm=0);

    void            update_ldu_row(long ac, long j, ROW& ldu);
    void            update_ldu_row(long ac, VTR_SPARSE& ldu);
    ROW             get_rowL(const long& r)const;
    ROW             get_rowU(const long& r)const;
    //ZBYTECNE JEN TESTOVACI:
    void            print_sparse_row(const long& N);
    //**************************************************************
    //********************** Virtualni metody **********************
    //**************************************************************
};

#endif // SPARSE_H
