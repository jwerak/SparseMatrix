#include "row.h"
#include <stdio.h>


//   ===============================================================================
//                                  PUBLIC
//   ===============================================================================

//    ===================== PUBLIC OPERATORS ================================
// ==============================================================================
/**
  Vytiskne vsechny prvky tak jak jdou za sebou vcetne nul
  */
void ROW::print_full()
{
    long ac=0;
    for (long i=0;i<full_length;i++){
        if (index[ac] == i && ac<value_number){
            printf("%f\t",values[ac]);
            ac++;
        }else{
            printf("0\t");
        }
    }
    printf("\n");
}


// ==============================================================================
/**
  \brief Ziska alokovane pole, bud najde dost. dlouhe, nebo zalokuje nove a nastavi ho jako vychozi.
  Pouzivat jen kdyz nechci pouzivat prvky, ktere byly v radku.
  */
void ROW::initialize(const long& l)
{
    if(allocated[a]){//je alokovane prvni pole?
        if(length[a] > l){//je pole dostatecne dlouhe?
            value_number=0;
            return;
        }else{//neni dost dlouhe... Co to druhe?
            if(allocated[na]){//je vubec druhe alokovane?
                if (length[na] > l){
                    switch_fields();
                    value_number=0;
                    return;
                }
            }else{//neni alokovane, zalokuj
                switch_fields();
                alloc(l);
                return;
            }
        }
    }else{//neni jeste alokovane prvni pole
        alloc(l);
        return;
    }
}


// ==============================================================================
/**
  */
void ROW::operator+=(const ROW& b)
{
    prepare_field((b.value_number + value_number));//orig pole bude volne k zapisu
    merge_fields(b.values, val[na], b.index, in[na], b.value_number, value_number, values, index, value_number);
}


// ==============================================================================
/**
  \param b - hodnota, kterou se roznasobi dany radek
  */
const ROW ROW::operator*(const double& b)
{
    double* vv = new double[value_number];
    long*   ii = new long[value_number];

    for (long i=0; i<value_number; i++)
    {
        vv[i] = values[i] * b;
        ii[i] = index[i];
    }
    return ROW(vv, ii, value_number, full_length);
}


// ==============================================================================
/**
*/
double ROW::operator*(const ROW& a)
{
    double result = 0;
    long ia=0;
    long i=0;

    while(ia<a.value_number && i<value_number){//udelej pro kazdou hodnotu nebo dokud nedojdes na konec "a"
        while(a.index[ia] < index[i])//popojed v acku na sloupec, ktery neni mensi
            ia++;
        if(ia >= a.value_number)
            break;

        if(a.index[ia]==index[i]){//neni nahodou i stejne velky?
            result += a.values[ia]*values[i];
            ia++;
        }
        i++;
    }
    return result;
}


// ==============================================================================
/**
  \brief Roznasobi radek s vektorem, pouzivam napriklad u podfunkci solve_Axb
  */
double ROW::operator* (const VTR& b)
{
    double result = 0;
    long i;

    for(i=0;i<value_number;i++)
        result += values[i] * (b.val[index[i]]);

    return result;
}


// ==============================================================================
/**
  */
void ROW::operator/=(const double& b)
{
    for (long i=0; i<value_number; i++)
        values[i] /= b;
}


//   ===============================================================================
//                                  PRIVATE
//   ===============================================================================
//    ===================== PRIVATE CONSTRUCTORS ================================
//   ===============================================================================
/**
  \brief prevezme pole a vytvori z nich novy radek
  \param v - pole hodnot
  \param l - pole indexu sloupce
  \param v_num - pocet hodnot v radku
  \param r_len - celkovy pocet sloupcu v radku
  \param d - pozice diagonalniho prvku v ridkem radku
  */
ROW::ROW(double* v, long* l, long v_num, long r_len):
a(0),na(1),value_number(v_num)
{
    allocated[a]    = true;
    allocated[na]   = false;
    full_length     = r_len;
    length[a]       = v_num;
    length[na]      = 0;
    values = val[a] = v;
    index  = in[a]  = l;
}


// ==============================================================================
//    ===================== PRIVATE METHODS ================================
// ==============================================================================

// ==============================================================================
/**
  \brief Pokud neni k dispozici dostatecne dlouhe pole, alokuje nove.
  */
void ROW::alloc(long n)
{
    if(allocated[a]){//je prvni pole alokovane?
        if(length[a]<n){//neni moc kratke?
            length[a] = n*enlarge;
            delete[] values;
            delete[] index;
            values = val[a] = new double[length[a]];
            index  = in[a]  = new long[length[a]];
        }
    }else{//poprve zalokuj pole
        length[a]=n*enlarge;
        values = val[a] = new double[length[a]];
        index  = in[a]  = new long[length[a]];
        allocated[a]  = true;
    }
}


// ==============================================================================
/**
  */
void ROW::dealoc()
{
    if(allocated[a]){
        delete[] val[a];
        delete[] in[a];
        allocated[a]=false;
    }

    if(allocated[na]){
        delete[] val[na];
        delete[] in[na];
        allocated[na]=false;
    }
}


// ==============================================================================
/**
\brief zajisti aby nove bylo k dispozici pole orig o minimalni velikosti n, zachova puvodni hodnoty ve druhem poli val[na]
\param n - vyzadovana delka pole(max. pocet ktery budu potrebovat)
  */
void ROW::prepare_field(const long& n)
{
    switch_fields();
    alloc(n);
}


// ==============================================================================
/**
  \brief Prohodi pole orig a addit
  */
void ROW::switch_fields()
{
    switch (a){
    case 0:
        a=1;
        na=0;
        break;
    case 1:
        a=0;
        na=1;
        break;
    default:
        throw "ROW.switch_fields - Je nastaveno nespravne pole";
    }
    values = val[a];
    index  = in[a];
}


// ==============================================================================
/**
  \brief Vytvori presnou kopii dat ze vstupniho radku.
  */
void ROW::copy_row(const ROW& b)
{
    alloc(b.value_number);

    full_length = b.full_length;
    value_number = b.value_number;

    for (long i=0; i<value_number; i++){
        values[i]=b.values[i];
        index[i]=b.index[i];
    }
}


// ==============================================================================
/**
  Spoj pole v1,v2 do v_out a i1, i2 do i_out a pocet setridenych hodnot i_num.
  Spojuje podle i1 a i2, pocita, ze pole jsou setridena
  \return v_out-spojene double hodn
  \return i_out indexy sloupcu v_out
  \return i_sum pocet spojenych hodnot (delka pole)
  \return k-celkovy pocet hodnot ve spojenem poli
*/
void ROW::merge_fields(const double* v1, const double* v2, const long* i1, const long* i2, const long& i1_num, const long& i2_num, double*& v_out, long*& i_out, long& i_sum)
{
    long i = 0; /// \param i - akt umisteni indexu [0]
    long j = 0; /// \param j - akt umisteni indexu [1]
    long k = 0; /// \param k - akt pocet vystupu

    // Dokud nedojdes na konec obou poli indexu...
    while ((i<i1_num && j<i2_num)){
        //======================================================
        //???Ktera hodnota z indexu na ktere se divam je vetsi???
        if(i1[i] < i2[j]){// Prvni index bude zapsan
            i_out[k] = i1[i];
            v_out[k] = v1[i];
            i++;
        }else{
            if(i2[j] < i1[i]){// Druhy index bude zapsan
                i_out[k] = i2[j];
                v_out[k] = v2[j];
                j++;
            }else{// Oba indexy jsou stejne -> secist values
                i_out[k] = i1[i];
                v_out[k] = v1[i]+v2[j];
                i++; j++;
            }
        }
        k++;
    }

    while (i<i1_num){
        i_out[k] = i1[i];
        v_out[k] = v1[i];
        i++; k++;
    }
    while (j<i2_num){
        i_out[k] = i2[j];
        v_out[k] = v2[j];
        j++; k++;
    }
    i_sum = k;
}


