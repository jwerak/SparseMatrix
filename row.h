#ifndef ROW_H
#define ROW_H

#include "vtr.h"

/**
  \brief Trida urcena pro radkove operace mezi maticemi.

  Pri praci s private hodnotamy se pri kazdem dalsim zapisu
  pise do hodnot orig. proto se pouziva metoda set_field(const long& n) .
  */
class ROW
{
    friend class MTX;
    friend class SPARSE;

    //   ===============================================================================
    //                                  PUBLIC
    //   ===============================================================================
public:
     //    ===================== PUBLIC CONST/DEST ================================
    ROW(const long& b){initialize(b);}
    ROW(const ROW& b){copy_row(b);}
    ROW():a(0),na(1),value_number(0){allocated[a]=false;allocated[na]=false;length[a]=0;length[na]=0;}
    ~ROW(){dealoc();}
    //    ===================== PUBLIC METHODS ================================
    void print_full();//vytiskne vsechny prvky tak jak jdou za sebou vcetne nul.
    void initialize(const long& l);//pripravi radek tak aby se do nej dalo znovu zapisovat

    //   ===============================================================================
    //                                  PRIVATE
    //   ===============================================================================
private:
    //    ===================== PRIVATE CONSTRUCTORS ================================

    ROW(double* v, long* l, long v_num, long r_len);

     //    ===================== PRIVATE OPERATORS ================================
    void        operator+= (const ROW& b);
    void        operator=  (const ROW& b){copy_row(b);}
    const ROW   operator*  (const double& b);
    double      operator*  (const ROW& b);
    double      operator*  (const VTR& b);
    ROW         operator/  (const double& b);
    void        operator/= (const double& b);

    //    ===================== PRIVATE METHODS ================================
    void alloc(long n);///< Alokuje jedno dostatecne dlouhe pole orig
    void dealoc();
    void prepare_field(const long& n);///< Nastavi vse tak aby se dalo zapsat n prvku do pole orig
    void switch_fields();///< Stara se o prohozeni orig a addit poli a pointeru okolo
    void copy_row(const ROW& b);///< Zkopiruje hodnoty objektu ROW do sebe
    void merge_fields(const double* v1, const double* v2, const long* i1, const long* i2, const long& i1_num, const long& i2_num, double*& v_out, long*& i_out, long& i_sum);

    //    ===================== PRIVATE VARIABLES ================================
    static const double enlarge = 1; ///< O kolik bude zvetseno pole orig kdy nebude dostatecne velke v ROW::set_field(const long& n)
    static const int nfields = 2;///< Pocet poli do kterych se budou zapisovat hodnoty

    //uloziste
    double* val[nfields];///< Udrzuje hodnoty radku
    double* values;///< Pointer na aktualni pole val[a]
    long*   in[nfields];///< Udrzuje indexy sloupcu
    long*   index;///< Pointer na aktualni pole index[a]
//    double* diag;///< pointer na diagonalni prvek

    int a;///< Actual field
    int na;///< Not actual field

    //bool
    bool allocated[nfields]; ///< Jsou alokovana pole orig?

    //informacni
    long value_number;///< Aktualni pocet uchovavanych hodnot
    long full_length;///< Teoreticka maximalni delka radku
    long length[nfields];///< Potencionalni plna delka alokovaneho poli orig
};
#endif // ROW_H
