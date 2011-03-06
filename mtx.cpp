#include "mtx.h"
#include <stdio.h>

#define RECURSIVE_PERMUTATION
#define LDU_VEKTOR

// ==============================================================================
// ==============================================================================
                                //PUBLIC
// ==============================================================================
// ==============================================================================


// ==============================================================================
/**
  \brief Resi rovnici A*x=b
  \return x - vektor neznamych
  \param[in,out] LDU - matice do ktere bude zapsan ldu soucin
  \param[in] b - prava strana
*/
void MTX::solve_Axb(MTX& LDU, VTR& x, VTR& b){
    printf("\nMatice pred permutaci\n");
//    print_full_mtx();

    printf("\nPermutuji matici\n");
//    make_symbolic_factorization();
    //===========================
    make_permutation_ND();

    printf("\nAplikuji permutaci\n");

    apply_permutation();
    printf("\nPermutace aplikovana na matici\n");
    b.apply_permutation(MTX::P);
    printf("\nPermutace aplikovana na vektor b\n");
    //===========================
    LDU = *this;
    printf("\nMatice po permutaci\n");
//    print_full_mtx();
    printf("\nRozkladani matice\n");

//    try{
        make_ldu(LDU);
//    }catch(...){
//        std::cerr<<"Chyba: LDU Rozklad \n";
//        throw;
//    }
    printf("\n\nMatice rozlozena na LDU\n");
//    LDU.print_full_mtx();
    printf("\nVstupni vektor b:\n");
//    b.print();
    VTR z,y;
    LDU.solve_Lzb(b,z);
    printf("\nVysledny vektor z:\n");
//    z.print();
    LDU.solve_Dyz(z,y);
    printf("\nVysledny vektor y:\n");
//    y.print();
    LDU.solve_Uxy(y,x);
    printf("\nVysledny vektor x:\n");
    reverse_permutation_vector();
    apply_permutation();
    b.apply_permutation(MTX::P);
    x.apply_permutation(MTX::P);
//    x.print();
}


// ==============================================================================
/**
 \brief Najde permutaci podle Nested Dissection Method a prepermutuje matici.

 \brief Po kazdem rozpuleni si zapise dalsi pocatecni body pro pulky beg[a], ktere vytvoril a neskonci dokud
 neprojde vsechny pulky z predchozicho kroku beg[b].

\brief Toto se vykonava dokud neni zaplnen vektor permutaci.
 \in <b>Podfunkce pouzivaji</b> MTX::addit_col

Promenna <b>S</b>(oznamuje status jednotlivych vrcholu)
    - 0 - Pocatecni hodnota, jeste s ni nebylo manipulovano
    - Ve funkci findsep:
        - 1 - Je prohledavana zleva
        - -1 - Je prohledavana zprava
    - Ve funkci findedge:
        - 2 - hleda krajni bod poprve*k
        - -2 - hleda krajni bod podruhe*k
    - 3 - Vrchol jiz byl zadan k permutaci
*/
void MTX::make_permutation_ND()
{
    //Promenne pro funkce find_sep a find_edges - alokovany zde aby byly alokovany jen jednou
    long* fl[2];//pro find_sep
    long* fr[2];//pro find_sep
    long* temp;//Docasne alokovane pole, pouzivano uvnitr find_sep

    long ap = row_number-1;//aktualni permutace - pozice kam se ma zapsat nalezena permutace
    long i;//index

    if(!sa)
        S = new long[row_number];

#ifndef RECURSIVE_PERMUTATION
    long *beg[2];//beginning - zapisou se vsechny pocatecni hodnoty, ktere maji byt v pristim kroku zapsany
    long n[2];//number of beginnings - kolik vrcholu se naslo pri predchozim puleni a kolik ted
    long a  = 0;//index pole do ktereho se zapisuje
    long b  = 1;//index pole do ktereho se zapisovalo minule kolo a ted se jenom cte
    beg[0]  = new long[row_number];
    beg[1]  = new long[row_number];
#endif

    for (i=0;i<row_number;i++)
        S[i]=0;

    //Alokuj pole permutace jestli neni alokovano
    if(!pa){pa = true;  P = new long[row_number]; S = new long[row_number];}

    //Priprav pole pro MTX::find_sep(). jedno volani si vyzada 4 pole long o delce "row_number"
    alloc_addit_l(5);
    fl[0] = addit_long[0];
    fl[1] = addit_long[1];
    fr[0] = addit_long[2];
    fr[1] = addit_long[3];
    temp  = addit_long[4];

#ifdef RECURSIVE_PERMUTATION
    find_sep(0, ap, S, fl, fr,temp);
#else
    //Priprav hodnoty na prvni zapis:
    n[a]=0;
    n[b]=1;
    beg[b][0]=0;
    find_sep(0, beg[a], n[a], ap, S, fl, fr);

    while (ap >= 0){//Hlavni smycka pro nalezeni permutace:
        for(i=0;i<n[b];i++){//proved puleni na kazde nerozpulene (dostatecne velke) casti. Ne rekurzivne!!!
            find_sep(beg[b][i], beg[a], n[a], ap, S, fl, fr);
        }
        (a==0)?a=1:a=0;
        (b==0)?b=1:b=0;
        n[a] = 0;
    }

    delete beg[0];
    delete beg[1];
#endif
    printf("\nHODNOTY, KTERE NEBYLY OZNACENE:\n");
    for(i=0;i<row_number;i++){
        if(S[i]!=3){
            printf("NEOZNACENA HODNOTA: %li\n",S[i]);
//            printf("NAVAZUJICI VRCHOLY: \n");
//            print_sparse_row(i);
//            P[i]=ap--;
        }
    }
    delete[] S;
    printf("\n\n Zbyvajicich hodnot k prepermutovani: %li\n", ap+1);
//    printf("Tisk permutace: \n\n");
//    for (i=0;i<row_number;i++)
//        printf("%li ",P[i]);
}


// ==============================================================================
/**
  \brief Rozklada se na LDU soucin
  \return Rozlozena matice ldu
  \param[in,out] LDU - matice do ktere bude zapsan ldu soucin
  \todo Predelat aby se matice LDU vytvorila zde a nemusela vytvaret externe - vytvor MTX copy(MTX)
*/
void MTX::make_ldu(MTX& LDU)
{
    LDU.value_number_orig = LDU.value_number;
    LDU.na = LDU.ns = LDU.nm = LDU.nd = 0;
    long j=0;//chodi po ridkem, prave vytvarenem, radku
    long ac; //Actual column/index
    double d=0;  // Diagonalni prvek aktualniho radku

//    double l;  // Hodnota l pri nasobeni radku - jen pro prehlednost
//    ROW u; //nacte se z LDU

#ifdef LDU_VEKTOR
    VTR_SPARSE ldu(LDU.column_number);//radek, ktery se nacte z matice a na konci upravy zapise do LDU

    for(long i=0; i<LDU.row_number; i++){//pres radky
//        if(i%1000==0)
            printf("\nJe zpracovavan radek: %li",i);

        LDU.set_vtr_sparse(i,ldu);
//        printf("\n\nRADEK PO NACTENI:\n");
//        ldu.print_full();

        ac = ldu.get_first_ac();

//        try{
            while(ac<i && ac!=-2){//VYTVOR CAST L
                LDU.update_ldu_row(ac, ldu);
                ldu.get_next_ac(ac);
            }
//        }catch(...){
//            std::cerr<<"Chyba: cast L";
//            throw;
//        }

        //JSI NA DIAGONALNIM PRVKU
        d = ldu.values[i]; //PREDELAT

        //VYTVOR CAST U
        if(ac != -2)
            ldu.get_next_ac(ac);

        try{
            while(ac != -2 /*|| ac < (ldu.value_number)*/){
                ldu.values[ac] /= d;
                ldu.get_next_ac(ac);
            }
        }catch(...){
            std::cerr<<"Chyba: cast U";
            throw;
        }

        LDU.nd += ldu.value_number-j;//POCET OPERACI
        LDU.set_row(ldu,i);
//        printf("RADEK PO ROZLOZENI:\n");
//        ldu.print_full();

//        LDU.print_full_mtx();
        //======================
//        //SPOCTI POCTY OPERACI
//        LDU.nd++;
//        LDU.nm++;
////            LDU.nm += u.value_number;

#else
    ROW ldu;//radek, ktery se nacte z matice a na konci upravy zapise do LDU // PRO ROW::

    for(long i=0; i<LDU.row_number; i++){//pres radky
        std::cout<<"\nJe zpracovavan radek: " << i<<std::endl;

        ldu = get_row(i); // PRO ROW::
//        ldu.print_full();
        j = 0;

        ac = ldu.index[j];
        while(ac<i){//VYTVOR CAST L

            //PUVODNI METODA
//            l = ldu.values[j] /= LDU.get_diagonal(ac);
//
//            u = LDU.get_row(ac,ac+1);
//            d = LDU.get_diagonal(ac);
//            ldu += u*(-l*d);
//            j++;
            //===============

            //NOVY ZPUSOB - NEVYTVARI ROW "U"
            LDU.update_ldu_row(ac, j, ldu);
            j++;
            ac = ldu.index[j];
//            ldu.print_full();
            //===============

            //======================
            //SPOCTI POCTY OPERACI
            LDU.nd++;
            LDU.nm++;
//            LDU.nm += u.value_number;
        }

        //JSI NA DIAGONALNIM PRVKU
        d = ldu.values[j];
        j++;

        //VYTVOR CAST U

        LDU.nd += ldu.value_number-j;//POCET OPERACI

        while(j<ldu.value_number){
            ldu.values[j] /= d;
            j++;
            //======================
            //SPOCTI POCTY OPERACI
        }
//        ldu.print_full();
        LDU.set_row(ldu,i);
#endif
    }
}



// ==============================================================================
/**
  \brief Nacte cast L z matice A
 */
void MTX::get_mtx_L(MTX& A)
{
    long i,j;
    ROW r;
    initialize(A.get_row_number(),A.get_column_number(),A.get_values_number());

    for(i=0;i<A.get_row_number();i++){
        j=0;
        r = A.get_row(i,0);

        while(r.index[j] < i)
            j++;
        r.values[j]=1;//nastav diagonalu
        r.value_number= ++j;//zkrat radek at se zbytek neopisuje

        set_row(r,i);
    }
}


// ==============================================================================
/**
  \brief Nacte cast D z matice A
 */
void MTX::get_mtx_D(MTX& A)
{
    long i;
    double* d = new double[1];/// pro pouziti do funkce SPARSE::set_row(const double* v,const long* i,const long& l, const long& r)
    long*   ii= new long[1];/// pro pouziti do funkce SPARSE::set_row(const double* v,const long* i,const long& l, const long& r)
    initialize(A.get_row_number(),A.get_column_number(),A.get_values_number());

    for(i=0;i<A.get_row_number();i++){
        d[0] = A.get_diagonal(i);
        ii[0]=i;
        set_row(d,ii,1,i);
    }
}


// ==============================================================================
/**
  \brief Nacte cast U z matice A
 */
void MTX::get_mtx_U(MTX& A)
{
    long i;
    ROW r;
    initialize(A.get_row_number(),A.get_column_number(),A.get_values_number());

    for(i=0;i<A.get_row_number();i++){
        r = A.get_row(i,i);
        r.values[0]=1;
        set_row(r,i);
    }
}


// ==============================================================================
/**
*/
const long& MTX::get_column_number()const
{
    return column_number;
}


// ==============================================================================
/**
*/
const long& MTX::get_row_number()const
{
    return row_number;
}


// ==============================================================================
//                              OPERATORS
// ==============================================================================
// ==============================================================================
/**
  \brief Roznasobeni dvou matic, jen pro testovaci ucely. Velmi pomaly zpusob. Doma nezkouset.
  Projede radky a roznasobi hodnotu leve matice v danem sloupci c s c-tym radkem prave matice
*/
void MTX::multiply_matrices(MTX& B, MTX& R)
{
    long i,j,ac; // ac - actual column
    double v;//value for multiplication of row
    ROW l;//row of actual left matrix
    ROW r;//row of right matrix
    ROW b;//resultant row - will write into this
    R = B;

    for (i=0;i<row_number;i++){//projdi kazdy radek
//        printf("\nRadek: %li\.n",i);
        j=0;
        l = get_row(i,0);
        while(j < l.value_number){//projdi vsechny nenulove hodnoty v radku
            if(j==0){
                ac = l.index[j];
                r = B.get_row(ac,0);
                v = l.values[j];
                b = r*v;
            }else{
                ac = l.index[j];
                r = B.get_row(ac,0);
                v = l.values[j];
                b += r*v;
            }
            j++;
        }
        R.set_row(b,i);
    }
}


// ==============================================================================
/**
  \brief
*/
const VTR MTX::multiply(const VTR& b)
{
    ROW a;
    double* r = new double[row_number];//return vektor

    for(long i=0;i<row_number;i++){
        a = get_row(i,0);
        r[i] = a*b;
    }
    return VTR(r,row_number);
}

// ==============================================================================
// ==============================================================================
                                //PRIVATE
// ==============================================================================
// ==============================================================================


// ==============================================================================
/**
  \brief Prestavi MTX::P tak, aby se dalsi permutaci vse vratilo do puvodniho stavu
*/
void MTX::reverse_permutation_vector()
{
    long* new_P = new long[row_number];
    for(long i=0;i<row_number;i++){
        new_P[P[i]] = i;
    }
    delete[] P;
    P = new_P;
}


// ==============================================================================
/**
  \brief Resi soustavu L*z = b
  \param [in]  b - prava strana
  \param [out] z - Vysledny vektor
*/
void MTX::solve_Lzb(const VTR& b, VTR& z)
{
    double s;//semiresult
    ROW l;
    long i;
    z.initialize(column_number);
    for(i=0;i<row_number;i++){
        l=get_rowL(i);
        s=l*z;
        z.val[i]=(b.val[i]-s);
    }
}


// ==============================================================================
/**
  \brief Resi soustavu D*y = z
  \param z - ziskano z solve_Lzb(const VTR& b, VTR& z)
*/
void MTX::solve_Dyz(const VTR& z, VTR& y)
{
    y.initialize(row_number);
    for(long i=0;i<row_number;i++)
        y.val[i]=z.val[i]/get_diagonal(i);
}


// ==============================================================================
/**
  \brief Resi soustavu U*x = y
  \param y - ziskano z solve_Dyz(const VTR& z, VTR& y)
*/
void MTX::solve_Uxy(const VTR& y, VTR& x)
{
    double s;//semiresult
    ROW l;
    long i;
    x.initialize(column_number);
    for(i=(row_number-1); i>=0; i--){
        l=get_rowU(i);
        s=l*x;
        x.val[i]=(y.val[i]-s);
    }
}


// ==============================================================================
/**
    \brief Alokuje dodatecna pole addit_val a addit_col a vynuluje hodnoty
    \warning Je obezni, co nejdriv se ji zbavit a tam, kde byla pouzivana nahradit MTX::alloc_addit_d(), nebo MTX::alloc_addit_l()
*/
void MTX::alloc_addit()
{
    if(!addit_values){
        for (int i=0;i<addit_fields; i++){
            addit_val[i] = new double[column_number];
            addit_col[i] = new long[column_number];
        }
        addit_values=true;
    }
}


// ==============================================================================
/**
    \brief Alokuje n poli addit_long o velikosti MTX::row_number
*/
void MTX::alloc_addit_l(long n)
{
    if(n<=addit_long_number && addit_long_bool)//jestli jsou pozadovana pole dostupna
        return;
    if(addit_long_bool)//nejsou dostupna, ale jsou alokovana
        dealloc_addit_l();

    addit_long = new long* [n];
    for(long i=0;i<n;i++)
        addit_long[i] = new long[row_number];
    addit_long_bool = true;
    addit_long_number = n;
}


// ==============================================================================
/**
    \brief Alokuje n poli addit_double o velikosti MTX::row_number
*/
void MTX::alloc_addit_d(long n)
{
    if(n<=addit_double_number && addit_double_bool)//jestli jsou pozadovana pole dostupna
        return;
    if(addit_double_bool)//nejsou dostupna, ale jsou alokovana
        dealloc_addit_d();

    addit_double = new double* [n];
    for(long i=0;i<n;i++)
        addit_double[i] = new double[row_number];
    addit_double_bool = true;
    addit_double_number = n;
}


// ==============================================================================
/**
    \brief Odalokuje dodatecna pole addit_val a addit_col.
*/
void MTX::dealloc_addit()
{
    if(addit_values){
        for (int i=0;i<addit_fields; i++){
            delete addit_val[i];
            delete addit_col[i];
        }
        M_aloc = false;
    }
}


// ==============================================================================
/**
    \brief Odalokuje dodatecna pole addit_long
*/
void MTX::dealloc_addit_l()
{
    if(addit_long_bool){
        for(int i=0;i<addit_long_number; i++)
            delete addit_long[i];
        delete[] addit_long;
        addit_long_bool=false;
    }
}


// ==============================================================================
/**
    \brief Odalokuje dodatecna pole addit_double
*/
void MTX::dealloc_addit_d()
{
    if(addit_double_bool){
        for(int i=0;i<addit_double_number; i++)
            delete addit_double[i];
        delete[] addit_double;
        addit_double_bool=false;
    }
}


// ==============================================================================
// ==============================================================================
                                //PROTECTED
// ==============================================================================
// ==============================================================================

// ==============================================================================
/** \brief Bubble sort of a ordinary long*
  Shakersort from book Algorithm and Data Structures by N. Wirth (1995) applied
  \param n - Delka pole ktere se ma setridit
  \param a - Pole, ktere se ma setridit
  \param key - pole do ktereho se zaznamenaji zmeny
  \warning Neni radne otestovan
    */
void MTX::BB_sortl(long n, long* a, long* key)
{
    long L = 0; /// \param L -
    long R = n-1; ///  \param R -
    long k = R; /// \param k - index of last changed value
    long i;

    for(i=0;i<n; i++)
        key[i]=i;

    while (L<R){
        for (i=R; i>L; i--){
            if (a[i-1]>a[i]){
                swap(a[i-1],a[i]);
                swap(key[i-1],key[i]);
                k = i;
            }
        }
        L=k;
        for(i = L; i<=R; i++){
            if (a[i-1]>a[i]){
                swap(a[i-1],a[i]);
                swap(key[i-1],key[i]);
                k = i;
            }
        }
        R=k-1;
    }
}


// ==============================================================================
/** \brief Bubble sort of a ordinary double*
  Shakersort from book Algorithm and Data Structures by N. Wirth (1995) applied
  \warning Neni radne otestovan
    */
void MTX::BB_sortd(long n, double* a)
{
    long L = 0; /// \param L -
    long R = n-1; ///  \param R -
    long k = R; /// \param k - index of last changed value
    long i;

    while (L<R){
        for (i=R; i>L; i--){
            if (a[i-1]>a[i]){
                swap(a[i-1],a[i]);
                k = i;
            }
        }
        L=k+1;
        for(i = L; i<=R; i++){
            if (a[i-1]>a[i]){
                swap(a[i-1],a[i]);
                k = i;
            }
        }
        R=k-1;
    }
}


// ==============================================================================
/** \brief Setridi pole *a o velikosti n podle klice *k. Funguje jen na pole mensi nez plna delka radku
*/
void MTX::update_field_L_by_key(const long n, long* a, long* k)
{
    int i;
    alloc_addit_l(1);

    for (i=0;i<n;i++)//srovnej hodnoty do jineho pole
        addit_long[0][k[i]]=a[i];

    for (i=0;i<n;i++)//zkopiruj hodnoty do puvodniho pole
        a[i] = addit_long[0][k[i]];
}


// ==============================================================================
/** \brief Setridi pole *a o velikosti n podle klice *k. Funguje jen na pole mensi nez plna delka radku
*/
void MTX::update_field_D_by_key(const long n, double* a, long* k)
{
    long i;

    alloc_addit_d(1);
    for (i=0;i<n;i++)//srovnej hodnoty do jineho pole
        addit_double[0][i] = a[k[i]];

    for (i=0;i<n;i++)//zkopiruj hodnoty do puvodniho pole
        a[i] = addit_double[0][i];
}

