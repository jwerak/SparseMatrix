#include <fstream>

#include "sparse.h"
#include "vtr.h"

#define LOAD_1 //ZNACI, ZE VSTUPNI SOUBOR ZACINA AZ 1. RADKEM
#define RECURSIVE_PERMUTATION
//#define WRITE_DETAILS
//#define MAKE_WAVE_THICK
#define NEW_FINDSEP


// ==============================================================================
// ==============================================================================
                                //PUBLIC
// ==============================================================================
// ==============================================================================

// ==============================================================================
// ==============================================================================
                            //VIRTUAL
// ==============================================================================
// ==============================================================================


// ==============================================================================
/**
 */
void SPARSE::operator=(const MTX& B)
{
    SPARSE* A = (SPARSE*)&B;
    copy_values((*A).row_number,(*A).column_number,(*A).value_number,(*A).value_length,(*A).value_start,(*A).column_start,(*A).diagonal);
}


// ==============================================================================
/**
  \return Objekt tridy ROW
  \param r-cislo radku, ktery ma byt nacten
  \param c-index prvniho sloupce, ktery bude nacten
  \param e-index posledniho sloupce, ktery bude nacten
  \return Objek typu ROW. radek pozadovanych vstupnich parametru
 */
ROW SPARSE::get_row(const long r, const long c)const
{
    //zjisti ktery sloupec se ma jako prvni predat
    long i=0;///< chodi po radku matice
    long j=0;///< chodi po novem radku, nemusi byt stejny jako i
    long e;///< index posledniho nacteneho clena
    //nastav zacatek
    while(column_start[r][i] < c)
        i++;
    //nastav konec
    e=value_length[r];

    double* tempv = new double[e];
    long*   tempc = new long[e];

    for (j=0; i<e; i++){
        tempv[j]=value_start[r][i];
        tempc[j]=column_start[r][i];
        j++;
    }
    return ROW(tempv, tempc, j, column_number);
}


// ==============================================================================
/**
  \brief Ziskej z matice radek, ktery zapises plne do radku VTR_SPARSE. ostatni hodnoty budou nula
  \param r - kolikaty radek se bude zapisovat
  \param [out] vtr - Do ktereho objektu se bude zapisovat
 */
void SPARSE::set_vtr_sparse(const long r, VTR_SPARSE& vtr)const
{
    vtr.set_val_from_sparse(value_start[r],column_start[r],value_length[r]);
}


// ==============================================================================
/**
  \brief Ziskej z matice radek ktery odpovida radku z casti L ldu rozkladu
  \param r-pozadovane cislo radku
 */
ROW SPARSE::get_rowL(const long& r)const
{
    long i=0;
    double* v = new double[value_length[r]];
    long* c   = new long[value_length[r]];
    while(column_start[r][i]<r){
        v[i]=value_start[r][i];
        c[i]=column_start[r][i];
        i++;
    }
    v[i]=1;
    c[i]=r;
    i++;
    return ROW(v,c,i,column_number);
}


// ==============================================================================
/**
  \brief Ziskej z matice radek ktery odpovida radku z casti U ldu rozkladu
  \param r-pozadovane cislo radku
 */
ROW SPARSE::get_rowU(const long& r)const
{
    long i=0;
    long j=0;
    //vystupni pole
    double* v = new double[value_length[r]];
    long* c   = new long[value_length[r]];
    //najdi zacatek radku U
    i = diagonal[r]-value_start[r];
    //zapis 1 na diagonalu
    v[j]=1;
    c[j]=r;
    i++;
    //zkopiruj zbyvajici prvky
    for(;i<value_length[r];i++){
        j++;
        v[j]=value_start[r][i];
        c[j]=column_start[r][i];
    }
    j++;
    return ROW(v,c,j,column_number);
}

// ==============================================================================
/**
  \brief Vytiskne N-ty radek
  \param N - kolikaty radek chci tisknout
 */
void SPARSE::print_sparse_row(const long& N)
{
    for(long i=0;i<value_length[N];i++)
        printf("%li ",column_start[N][i]);
    printf("\n");
}


// ==============================================================================
/**
  \brief Ulozi radek ROW do r-teho radku
 */
void SPARSE::set_row(const ROW& a, const long& r)
{
    if(a.value_number > value_length[r]){
        value_start[r] = M_get_free_place(a.value_number);
        long s = value_start[r]-value[allocated_field];
        column_start[r] = column[allocated_field]+s;
    }
    copy_row_in(a.values,a.index,r,a.value_number);
}


// ==============================================================================
/**
  \brief Ulozi radek VTR_SPARSE do r-teho radku
 */
void SPARSE::set_row(VTR_SPARSE& a, const long& r)
{
    long i, ac;
    long e = a.value_number - value_length[r];//prodlouzeni radku po pricteni hodnot
    double *v;//pointer na novy zacatek radku

    if(a.value_number > value_length[r]){
        v = M_get_free_place(a.value_number);
        value_start[r] =v;
        long s = value_start[r]-value[allocated_field];
        column_start[r] = column[allocated_field]+s;
    }
    if(a.value_number>value_length_max)
        value_length_max = a.value_number;
    value_length[r] = a.value_number;
    value_number += e;

    i = 0;
    ac = a.get_first_ac();
    while(ac != -2){
        value_start[r][i] = a.values[ac];
        column_start[r][i] = ac;
        if(ac == r){//Nasel jsi diagonalu?
            diagonal[r] = &value_start[r][i];
        }
        i++;
        a.get_next_ac(ac);
    }
}

// ==============================================================================
/**
  \param v - pole hodnot
  \param i - pole indexu
  \param l - delka radku
  \param r - cislo radku do ktereho budu zapisovat
 */
void SPARSE::set_row(const double* v,const long* i,const long& l, const long& r)
{
    if(l > value_length[r]){
        value_start[r] = M_get_free_place(l);
        long s = (value_start[r]-value[allocated_field]);
        column_start[r] = (column[allocated_field]+s);
    }
    copy_row_in(v,i,r,l);
}


// ==============================================================================
/**
Load matrix from file
 */
/** \todo       - Zjistovat jestli je matice dost velka???
*/
void SPARSE::loadFromFile(string cesta){
    long i,j,k;//indexy pro chozeni po radcich, polich, atd...

    fstream vstup;
    vstup.open(cesta.c_str());

    // ==========================================
    //zjisti rozmery a pocty
    vstup>>value_number;
    vstup>>row_number;
    vstup>>column_number;

    alloc(row_number, value_number);
    for(i=0;i<row_number;i++)
        value_length[i]=0;

    // ==========================================
    //zatim budou ukazovat na vytvorene pole indexu radku
    long** row_start = new long*[row_number+1];
    long* row_ind = new long[value_number+1];
    long* column_ind = new long[value_number+1];

    // ==========================================
    //nakopiruj hodnoty a zjisti teor. delky radku
    for (i=0; i<value_number; i++){
        vstup>>row_ind[i];
        vstup>>column_ind[i];
#ifdef LOAD_1
        row_ind[i]--;
        column_ind[i]--;
#endif
        vstup>>value[0][i];
        value_length[row_ind[i]]++;
    }

    row_start[0] = row_ind;
    value_start[0] = value[0];
    for (i=1; i<=row_number; i++){
        row_start[i]   = row_start[i-1] + value_length[i-1];
        value_start[i] = value_start[i-1] + value_length[i-1];
    }

    // ==========================================
    //setrid hodnoty ve stejnem radku za sebe
    order_values_by_row(column_ind, row_ind, value[0], row_start, row_number);
    
    // ==========================================
    //serad radky podle sloupcu v jednotlivych radcich
    alloc_addit();
    for (i=0; i<row_number; i++){
        BB_sortl(value_length[i], &column_ind[row_start[i]-row_ind],addit_col[0]);
        update_field_D_by_key(value_length[i], value_start[i],addit_col[0]);
    }

    // ==========================================
    //secti opakujici se hodnoty
    long cp,rp; //column, row - previous
    value_start[0]  = value[0];
    column[0][0] = column_ind[0];
    column_start[0] = column[0];
    cp = column[0][0]= column_ind[0];
    rp = 0;
    if(column[0][0]==0)
        diagonal[0] = &value_start[0][0];
    else
    {
        printf("loadFromFile: Neni diagonalni prvek na pozici 0,0 !!!");
        return;
    }
    if (row_ind[0] != 0)
    {
        printf("loadFromFile: Matice nema prvni radek, nebo neni spravne setridena");
        return;
    }

    for(i=1,j=0,k=0;;){//i-chodi po puvodnich vstupech, j-chodi po radcich, k-chodi po sloupcich v radku
        while(cp==column_ind[i]){
            value_start[j][k] += value[0][i];
            i++;
            if(i==value_number)
                break;
        }
        cp = column_ind[i];
        k++;//v danem miste uz nic zapisovat nebudu
        if(rp != row_ind[i]){//prehoupnul jsem se na dalsi radek
            if(i==value_number)
                break;
            k = 0;
            j = rp = row_ind[i];
            value_start[j] = value_start[j-1]+value_length[j-1];
            column_start[j] = column_start[j-1]+value_length[j-1];
        }
        value_start[j][k] = value[0][i];
        column_start[j][k] = column_ind[i];
        if(j==column_start[j][k])
            diagonal[j] = &value_start[j][k];
        i++;
        if(i==value_number)
            break;
    }
    for(i=0;i<row_number;i++)
        if(value_length_max<value_length[i])
            value_length_max = value_length[i];
    delete[] row_start;
    delete[] row_ind;
    delete[] column_ind;
}


// ==============================================================================
/**
\brief Tisk na standardni vystup, vytiskne plnou matici
*/
void SPARSE::print_full_mtx()
{
    long ar=0;// aktualni radek - ridky
    long ac = 0; //aktualni sloupec v radku - ridky
    printf("\n");
    for (long i=0; i<row_number; i++){
        for (long j=0; j<column_number;j++){
            if (column_start[ar][ac]== j && ac < value_length[ar]){
                printf("%1.3f \t",value_start[ar][ac]);
                ac++;
            }else{
                printf("0 \t");
            }
        }
        printf("\n");
        ar++;
        ac = 0;
    }
}


// ==============================================================================
/**
\brief
*/
void SPARSE::write_log(double t)
{
    FILE* out;
    string cesta = "log_Axb.out";
    out = fopen(cesta.c_str(),"w");
    fprintf (out,"Zapis z reseni soustavy A*x=b\n");
    fprintf (out,"Matice prerozdelena poomoci metody Nested dissection a rozlozena na ldu soucin\n");
    fprintf (out,"Pocet nenulovych prvku pred rozkladanim:%li\n",value_number_orig);
    fprintf (out,"Pocet nenulovych prvku po rozlozeni:%li\n",value_number);
    fprintf (out,"Pocet operaci:\n NASOBENI: %li\n DELENI: %li\n SCITANI: %li\n ODCITANI: %li\n",nm,nd,na,ns);
    fprintf (out,"Run time: %1.3lf min\n", t);
    fclose(out);
}

//********************************************************************
//********************************************************************
                        //VIRTUALNI
//********************************************************************
//********************************************************************


// ==============================================================================
// ==============================================================================
                                //PRIVATE
// ==============================================================================
// ==============================================================================

// ==============================================================================
//                              PRACE S MATICI
// ==============================================================================
/**
Procedura - alokuje pole pouzita v matici. Neni blbuvzdorna - jestli byla alokovana, smaze se a prealokuje
\param r/c    - row, column, two fields
\param ph     - pocet hodnot
\param rm     - rozsireni matice
*/
void SPARSE::alloc(const long& r, const long& nov, const int& rm)
{
    if (M_aloc == true) dealloc();
    M_aloc = true;
    M_allocated_field_length[0]     = enlarge_field(nov, rm);
    extra_field_length              = (long)M_allocated_field_length[0]*extra_field_enlargement/100;
    allocated_field                 = 0; // pri kopirovani se vytvori jenom jedno pole do ktereho se vse zkopiruje

    value[allocated_field]          = new double [M_allocated_field_length[0]];
    column[allocated_field]         = new long  [M_allocated_field_length[0]];
    value_start                     = new double* [r+1];
    column_start                    = new long* [r+1];
    value_length                    = new long  [r];
    diagonal                        = new double* [r];
    addit_values                    = false; //jsou alokovana doplnkova pole addit_val a addit_col?

    M_empty_start[allocated_field]  = &value[0][0];
    M_empty_end[allocated_field]    = &value[0][M_allocated_field_length[0]-1];
}


// ==============================================================================
/**
  \brief Najde pas, ktery puli graf matice na dve temer stejne casti.
Nalezenoy pas zapise do pole MTX::P.
  \brief Permutace zapisuje primo zde, protoze pouziva jenom jeden procesor, jinak by musel jako
  navratovou hodnotu pouzivat hodnoty a delku ziskaneho oddelovace a vse se muselo jeste jednou zkopirovat

  \param[in] p - Vychozi vrchol casti, ktera se ma pulit
  \param[in,out] n - next - Pole vychozich prvku, rozdelenych casti - zde se zadaji max 2, min 0
  \param[in,out] in - Index of n - kolik hodnot je jiz zadano/ na kterou se zapisuje dalsi nova
  \param[in,out] ap - Aktualni permutace - na ktere misto se zapisuje nalezeny oddelovac
  \param[in,out] S  - Informuje o stavu vrcholu viz. MTX::make_permutation_ND()
  \param[in] fl,fr - Alokovana pole aby nedochazelo ke zbytecnemu prealokovani
  \param[in] fr, fl - mistni uchovava informace o vrcholech ktere je treba projit- found left/right
 */
#ifdef RECURSIVE_PERMUTATION
void SPARSE::find_sep(long p, long& ap, long* S, long** fl, long** fr, long* TMP){
#else
void SPARSE::find_sep(long p, long* n, long& in, long& ap, long* S, long** fl, long** fr){
#endif
    const long min_num = 6;//Jestli je v poli mene prvku nez min_num, budou vrcholy zapsany do P
    const long min_nw  = 6;//Jestli pri prochazeni je mene nez min_nw, tak uz pulku dal neprochazej
    const long sf = 3;//Zadej kdyz najdes vlnu z druhe strany
    const long ct = 2;//Jeste jsem nenarazil na zadnou vlnu
    long str = 1;//Zadej pri vyhledavani stredu - zadava se pri pruchodu zprava
    long stl = -1;//Zadej pri vyhledavani stredu - zadava se pri pruchodu zleva
    long nvl[2];// pocet hodnot zleva
    long nvr[2];// pocet hodnot zprava
    long a = 0;// index fl a fr do kterych se zapisuje
    long b = 1;// index fl a fr ze kterych se cte(byli ziskany v minule vlne)
    long rv[2];// Do tohoto pole budou zapsany vystupni vrcholy z find_edges()
    long nv;//number of vertices - pocet hodnot ktere jsou v casti, ktera zacina v iv
    long nw;//number of waves - kolik vln se vytvorilo
    long i;//obecny index

#ifdef RECURSIVE_PERMUTATION
    long in=0;//Pocet rozpuleni
#endif

    //===================================
    //Najdi krajni "nejvzdalenejsi" vrcholy
    //===================================
    nv = find_edges(fl, TMP, p, S, rv, nw);
    if(nv<=min_num || nw<= min_nw){
        for(i=0; i<nv; i++){
            P[TMP[i]] = ap--;
            S[TMP[i]] = 3;
        }
        return;
    }
    //===================================
    //Rozdeleni casti
    //===================================
    //Inicializace:
    nvl[a]  = 0;
    nvl[b]  = 1;
    nvr[a]  = 0;
    nvr[b]  = 1;

    fl[b][0]= rv[0];//nalezene prvky
    fr[b][0]= rv[1];

    S[fl[b][0]] = stl;
    S[fr[b][0]] = str;

    //Hlavni smycka - pracuje dokud jsou ve vlne nalezeny nove vrcholy
    //POSLI NAPROTI SOBE DVE VLNY
    while(nvl[b]>0 || nvr[b]>0){
        nvl[a] = make_wave(S, ap, fl[b], fl[a], nvl[b], stl, sf, ct);
        nvr[a] = make_wave(S, ap, fr[b], fr[a], nvr[b], str, sf, ct);
        //Zamen pole
        swap(a,b);
        nvl[a]=nvr[a]=0;
    }
    //ZJISTI JESTLI NEKTERA Z HODNOT NEZUSTALA VE SLEPE ULICCE, JESTLI ANO, OZNAC JI JAKO STRED
    for(i=0;i<nv;i++){
        if(fabs(S[TMP[i]]) == ct){
            S[TMP[i]] = sf;
            P[TMP[i]] = ap--;
        }
    }
    //Po nalezeni dalsich polovin zavolej na kazdou polovinu pulici funkci
    in = 2;
#ifdef RECURSIVE_PERMUTATION
    for(i=in-1; i>-1; i--){
        find_sep(rv[i],ap,S,fl,fr,TMP);
    }
#endif
}


// ==============================================================================
/**
  \brief Procedura - Zprehazi sloupce a radky podle permutace MTX::P
*/
void SPARSE::apply_permutation()
{
    long i,j;
    size_t s;//velikost pole pro kopirovani pomoci memcpy
    long** new_col_start;
    long*  new_value_length;
    double** new_val_start;
    long* key = new long[row_number];//Klic podle ktereho se spravne seradi radek

    //===================================
    //Zamen sloupce
    //===================================
    for(i=0;i<row_number;i++){
        for(j=0;j<value_length[i];j++){
            key[j] = P[column_start[i][j]];
        }
        s = sizeof(long)*value_length[i];
        memcpy(column_start[i],key,s);
    }

    //===================================
    //Prohazej radky
    //===================================
    //Alokuj nove pole pointeru na zacatky radku
    new_col_start = new long* [row_number];
    new_val_start = new double* [row_number];
    new_value_length = new long [row_number];
    //zapis do nich prepermutovane zacatky radku
    for(i=0; i<row_number; i++){
        new_col_start[P[i]]    = column_start[i];
        new_val_start[P[i]]    = value_start[i];
        new_value_length[P[i]] = value_length[i];
    }

    //Smaz stare zacatky
    delete[] value_start;
    delete[] column_start;
    delete[] value_length;

//    //Prehazej nove vytvorene na spravna mista
    value_start = new_val_start;
    column_start = new_col_start;
    value_length = new_value_length;

    //===================================
    //Setrid hodnoty v radku podle sloupce
    //===================================
    for(i=0;i<row_number;i++){
        BB_sortl(value_length[i],column_start[i], key);
        update_field_D_by_key(value_length[i],value_start[i],key);
        for(j=0;j<value_length[i];j++)
            if(column_start[i][j]==i)
                diagonal[i]=&value_start[i][j];
    }
    delete[] key;
}



// ==============================================================================
/**
  \brief Najde dva body v poli ktere jsou od sebe "nejvzdalenejsi"
  \param[in]    f - Pole pro zapisovani hodnot z kazde vlny. Alokovano je jinde aby se nemuselo zbytecne alokovat
  \param[out]   ff - Pro zapisovani platnych vrcholu ktere pouziji - seznam vrcholu v dane casti
  \param[in]    iv - Initial value. Odkud zacne prvni vyhledavani
  \param[in,out] S - Informace o stavu vrcholu viz. MTX::make_permutation_ND()
  \param[out]   rv[2] - Resulting value - nalezene nejvzdalenejsi body
  \param[out]   nw - Number of waves - kolikrat jsem projizdel graf
  \return       nv - pocet hodnot ktere jsou v dane casti grafu
*/
long SPARSE::find_edges(long** f, long* ff, long iv, long* S, long* rv, long& nw)
{
    long nvv[2]= {0};//Number of vertices - kolik hodnot nacetl do f - v jedne vlne
    long nv;//Number of vertices - V celem jednom prohledani
    long ar,ac;//Actual row/column
    long a = 0;
    long b = 1;
    long ii[2]={-2,2};//Identifikator pro plneni pole S
    long iia = 0;//Index pro prave zapisovanou hodnotu do S
    long iib = 1;//Index pro minule zapisovanou hodnotu do S
    long i,j,k;//Obyc indexy
    const long nop = 2; //Number of passes

    if(S[iv] == 3)
        printf("\nVstup do find edges je 3-ka\n");

    //===================================
    //Main loop:
    //===================================
    for(i=0;i<nop;i++){//Projed nop krat
        //Initialize:
        nv      = 0;
        nw      = 0;
        ff[nv++]= iv;
        f[b][0] = iv;
        nvv[b]  = 1;
        S[iv]   = ii[iia];
        //Loop:
        while(nvv[b]>0){//Prohledej vlnu za vlnou dokud se v predchozi vlne nasli nove vrcholy
            nw++;
            nvv[a]=0;
            for(j=0;j<nvv[b];j++){//Pro kazdy vrchol, ktery byl nalezen v predchozi vlne
                ar = f[b][j];
                for(k=0;k<value_length[ar];k++){//Najdi vrcholy, ktere tvori nasledujici vlnu
                    ac = column_start[ar][k];
                    if(S[ac] != ii[iia] && S[ac]!=3){//Bude patrit k dalsi vlne
                        f[a][nvv[a]++] = ac;
                        ff[nv++] = ac;
                        S[ac]=ii[iia];
                    }
                }
            }
            //Indexy pro chozeni po f
            swap(a,b);
        }
        //Indexy pro zapis do S
        swap(iia,iib);
        //Dalsi vlna nebyla jiz vytvorena -> v aktualnim poli jsou nejvzdalenejsi vrcholy
        if(i < (nop-1))//Bude-li dochazet k dalsimu prohledani, zmen pocatecni hodnotu
            iv = ff[nv-1];//Zapis novou pocatecni hodnotu
    }
    rv[0]=iv;
    rv[1]=ff[nv-1];
    if(rv[0]==rv[1]){
        printf("\n !!!!!!!!!!!!!!\n rv[0]==rv[1]\n !!!!!!!!!!!!!!\n");
        printf("Pocet vln: %li\n",nw);
        printf("Pocet hodnot: %li\n",nv);
    }
    return nv;
}


// ==============================================================================
/**
 \brief Ze vstupnich vrcholu vytvori jednu vlnu. Kdyz narazi na vlnu z druhe strany,
tak zapise puleni.
 \param S - Pole do ktereho se zapisuji indexy o aktualnim stavu vrchlu viz. MTX::make_permutation_ND()
 \param ap - aktualni cislo permutace - pozice kam se ma zapsat nalezena permutace
 \param IW - Input wave - vstupni vrcholy predchozi vlny
 \param OW - Output wave - vystupni vrcholy zde vytvorene vlny
 \param NOI - Num of inputs - pocet vrcholu ze kterych se vytvari nove vlny
 \param ST - SET if True - cislo, ktere se zada na misto nove vytvarene vlny
 \param SF - SET if False - cislo, ktere se zada jestli nova vlna neni vytvorena (nejspis 3)
 \param CT - Check true - jestli je absolutni hodnota tohoto cisla, zadej ST
 \return Pocet nalezenych vrcholu ve vytvorene vlne
*/
long SPARSE::make_wave(long* S, long& ap, const long* IW, long* OW, const long& NOI, const long& ST, const long& SF, const long& CT)
{
    long i,j;
    long ar,ac; //aktualni radek/column
    long nov = 0;//num. of values
    //PRO VSECHNY VRCHOLY VE VSTUPNI VLNE
    for(i=0;i<NOI;i++){
        ar = IW[i];
        if(S[ar]!=SF){
            //PRO KAZDY SPOJENY VRCHOL
            for(j=0;j<value_length[ar];j++){
                ac = column_start[ar][j];
                //JESTLI JE VRCHOL VOLNY K ZAPISU NOVE VLNY
                if(fabs(S[ac]) == CT){
                    S[ac] = ST;
                    OW[nov++] = ac;
                }else{
                    //JINAK JESTLI JE TO VHODNE, ZAPIS TAM STREDOVY PROUZEK
#ifndef MAKE_WAVE_THICK
                    if(S[ac]==(-1*ST)){
                        S[ac] = SF;
                        P[ac] = ap--;
                    }
#else
                    //ZAPOVEZENO, ZPOMALUJE TO ROZKLAD VIC NEZ BEZ PERMUTACE
                    if((S[ac]==(-1*ST) || S[ac]==SF) && S[ar]!=SF){
                        S[ar] = SF;
                        P[ar] = ap--;
                    }
#endif
                }
            }
        }
    }
    return nov;
}


// ==============================================================================
/**
  \brief Rekurzivni funkce - prohledava graf matice a zjisti kolik radku a kterych bude v kterem radku
  \param LV - Least value - nejmensi hondota v radku
  \param AR - Actual row - cislo aktualniho radku
  \param IC - Input column - vrchol ze ktereho jsem byl sem poslan
  \param SNV - Set new value
  \param OL - Input length - jak byl puvodne dlouhy radek
  \param S - Rika jestli byl jiz vrchol prohledavan(S[IC]==SNV), nebo ne.
*/
long SPARSE::make_symbolic_wave(const long& LV, const long& AR, const long& IC, const long& SNV, const long& OL, long* S)
{
    long rl = 0;
//    if(IC == SNV)
//        return rl;
    long i,j,ac;
    //Pro kazdy vrchol
#ifdef WRITE_DETAILS
    printf("\nVSTUPUJI DO FUNKCE make_symbolic_wave:\nLV\tAR\tIC\tOL\n");
    printf("%li\t%li\t%li\t%li\n",LV,AR,IC,OL);
#endif
    for(i=0; i<OL; i++){
        //Jestli ten vrchol uz nebyl v teto funkci
        ac = column_start[IC][i];
#ifdef WRITE_DETAILS
        printf("\nac = %li\n",ac);
#endif
        if(S[ac]!=SNV){
            S[ac] = SNV;
            //Byl jiz vrchol symbolicky eliminovan?
            if(ac >= AR){//NEBYL
                rl++;
#ifdef WRITE_DETAILS
                printf("Hodnota NEBYLA eliminovana, jen se pricita\n");
#endif
            }else{//BYL ELIMINOVAN
                //Vyl eliminova, jestli je mensi, tak se pres nej jenom sklouznu na vrcholy s nim spojene
#ifdef WRITE_DETAILS
                printf("Hodnota BYLA drive eliminovana\n");
#endif
                if(ac >= LV){
                    rl++;
#ifdef WRITE_DETAILS
                    printf("Hodnota se pricita a spousti dalsi funkce\n");
#endif
                }
                for(j=0;j<value_length[ac];j++){
                    rl += make_symbolic_wave(column_start[AR][0],AR,column_start[ac][j],SNV,value_length[ac],S);
                }
            }
        }else{
#ifdef WRITE_DETAILS
            printf("Vrchol jiz byl zadan\n",ac);
#endif
        }
#ifdef WRITE_DETAILS
        printf("\n");
#endif
    }
    return rl;
}


// ==============================================================================
/**
  \brief Provede symbolicky rozklad matice
*/
void SPARSE::make_symbolic_factorization()
{
    long i,ol;
    long snv = 4;//set new value - vice o ni v make_symbolic_wave
    long rl = 0;//resulting length
    printf("\nSYMBOLICKA FAKTORIZACE\n");

    if(!pa){pa = true;  P = new long[row_number]; S = new long[row_number];}

    for(i=0;i<row_number;i++){
//        for(j=0;j<value_length[i];j++){
//            S[column_start[i][j]] = snv;
//        }
//        for(j=0;j<value_length[i];j++){
//            if(column_start[i][j]!=i)
#ifdef WRITE_DETAILS
        printf("\n==================================\nSpoustim radek: %li\n==================================\n",i);
#endif
        ol = value_length[i];
        rl = make_symbolic_wave(column_start[i][0],i,i,snv,ol,S);
        value_length[i] = rl;
//        }
        snv++;
//        make_symbolic_wave(value_start[i][0],i,value_start[i][0],sfv,S);
        printf("Value length %li = %li\n",i,value_length[i]);
    }
}


// ==============================================================================
/**
 \brief Allocate arrays and set matrix properties
 \param nor/c - num of rows/columns/values
 \param nov - kolik hodnot je ocekavano, ze bude matice obsahovat.
 \param rm - rozsireni matice v %, implicitne je nastaveno na 0% == 100% velikost nov. pouzije se v enlarge_field(const long& ph, const int& rm).
 */
void SPARSE::initialize(const long& nor, const long& noc, const long& nov, const int& rm)
{
    alloc(nor, nov, rm);
    row_number = nor;
    column_number = noc;
    value_number = 0;
    value_length_max = 0;
    for (long i=0;i<nor; i++){
        value_length[i]=0;
    }
}


// ==============================================================================
/**
    \brief Nahrazuje ldu += u*(-l*d); tak aby nedochazelo k nacitani u do noveho radku => Chaby pokus o zrychleni
    \param[in] ac - actual column
    \param[in] jj - Na kolikate pozici v radku jsem
    \param[in,out] ldu - Radek rozkladane matice, ktery se updatuje
    \warning - Zastarale - ROW nahrazeno jinou tridou
*/
void SPARSE::update_ldu_row(long ac, long jj, ROW& ldu)
{

    double d = *diagonal[ac];
    double l = ldu.values[jj] /= d;

    //PRIPRAV HONDOTY PRO PREDELANE MERGE FIELDS
    //IN_1 - PUVODNI HODNOTA RADKU ldu
    const double* v1  = ldu.values;
    const long*   i1  = ldu.index;
    const long i1_num = ldu.value_number;

    //IN_2 - HODNOTA U Z LDU MATICE
    const double*   v2      = (diagonal[ac]+1);
    long            s       = v2 - value_start[ac];//Vzdalenost zacatku v_2 od prvniho prvku v radku
    const long*     i2      = &column_start[ac][s];
    const long      i2_num  = (&value_start[ac][value_length[ac]]-v2);

    ldu.prepare_field(ldu.full_length);//orig pole bude volne k zapisu

    //OUT - VYSLEDNY RADEK ldu
    double* v_out = ldu.values;
    long*   i_out = ldu.index;

    //=========================
    //PREDELANE MERGE_FIELDS
    //=========================
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
                v_out[k] = v2[j]*(-l*d);
                j++;
            }else{// Oba indexy jsou stejne -> secist values
                i_out[k] = i1[i];
                v_out[k] = v1[i]+v2[j]*(-l*d);
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
        v_out[k] = v2[j]*(-l*d);
        j++; k++;
    }
    ldu.value_number = k;
    //*************************
    //PREDELANE MERGE_FIELDS
    //*************************
}


// ==============================================================================
/**
    \brief Nahrazuje ldu += u*(-l*d); tak aby nedochazelo k nacitani u do noveho radku => Chaby pokus o zrychleni
    \param[in] ac - actual column
    \param[in] jj - Na kolikate pozici v radku jsem
    \param[in,out] ldu - Radek rozkladane matice, ktery se updatuje
*/
void SPARSE::update_ldu_row(long ac, VTR_SPARSE& ldu)
{
    //DRIVE: lepe prehledne
//    double d = *diagonal[ac];
//    double l = ldu.values[ac] /= d;
//
//    const double*   v2      = (diagonal[ac]+1);
//    long            s       = v2 - value_start[ac];//Vzdalenost zacatku v2 od prvniho prvku v radku
//    const long*     i2      = &column_start[ac][s];
//    const long      i2_num  = (&value_start[ac][value_length[ac]]-v2);
//    ldu.add_val_from_sparse(v2, i2, i2_num, (-l*d));

    //NYNI:

    ldu.values[ac] /= *diagonal[ac];

    ldu.add_val_from_sparse((diagonal[ac]+1), (&column_start[ac][(diagonal[ac]+1) - value_start[ac]]), (&value_start[ac][value_length[ac]]-(diagonal[ac]+1)), (-(ldu.values[ac])*(*diagonal[ac])));
}


// ==============================================================================
/**
    \brief Dealokuje pole pouzita v matici.
*/
void SPARSE::dealloc()
{
    if (M_aloc){
        for (int i=0; i<=allocated_field; i++){
            delete[] value[i];
            delete[] column[i];
        }
        delete[] value_start;
        delete[] column_start;
        delete[] value_length;
        delete[] diagonal;
        M_aloc = false;
    }

    if(addit_values){
        for (int i=0;i<addit_fields; i++){
            delete addit_val[i];
            delete addit_col[i];
        }
        addit_values=false;
    }
}


// ==============================================================================
/**
\brief Zkopiruje vstupni hodnoty do matice
\param nor/c/v- num of rows/columns/values
\param v/c_pointer    - value, column pointers
\todo Predelat - bude dostavat jednotlive pole pro kazdy radek
 */
void SPARSE::copy_values(const long& nor, const long& noc, const long& nov, const long* v_length, const double* const* v_start, const long* const* c_start, const double* const* d)
{
    row_number      = nor;
    column_number   = noc;
    value_number    = nov;
    value_length_max=0;
    alloc(nor,nov,10);
    initialize_value_length(v_length);

    for (int i=0; i<row_number; i++){//pro kazdy radek
        value_length[i] = v_length[i];
        if(value_length_max<value_length[i])//nastav maxim delku radku
            value_length_max = value_length[i];

        for (int j=0; j<value_length[i]; j++){//nakopiruj vsechny hodnoty
            value_start[i][j]  = v_start[i][j];
            column_start[i][j] = c_start[i][j];
            if(v_start[i][j]==*d[i])diagonal[i]=&value_start[i][j];
        }
    }
}


// ==============================================================================
/**
\brief Zkopiruje vstupni hodnoty do radku r. Radek musi byt dostatecne dlouhy, zde se nekontroluje
\param v - Ridke pole hodnot, ktere se ma zapsat
\param ind - Pole sloupcu odpovidajicich v
\param r - Index radku do ktereho se ma zapisovat
\param l - nova delka radku
 */
void SPARSE::copy_row_in(const double* v, const long* ind,  const long& r, const long& l)
{
    //Predani informaci o radku
    long e = l - value_length[r];//prodlouzeni radku po pricteni hodnot
    value_length[r]=l;
    if(value_length[r]>value_length_max)
        value_length_max=value_length[r];

    value_number += e;

    //Predani hodnot
    for (long i=0;i<l;i++){
        value_start[r][i]   =v[i];
        column_start[r][i]  =ind[i];
        if(column_start[r][i]==r)diagonal[r]=&value_start[r][i];
    }
}


// ==============================================================================
/**
\brief Setridi pole c,r,v tak aby byly hodnoty r u sebe
\param rs - row-start - pointery na pole r kde zacina teoreticky dalsi radek
\param rl - pocet hodnot v kazdem radku
\param rn - pocet radku
 */
void SPARSE::order_values_by_row(long* c, long* r, double* v, long** rs, const long rn)
{
    long** re; //row_end - pole pointeru na posledni potvrzenou hodnotu v kazdem radku
    re = new long* [rn+1];
    long ri=0; //index aktualne posuzovaneho radku
    long i;//Index aktualne posuzovanych poli r,c,v

    for (i=0; i<rn; i++)
        re[i] = rs[i];

    for(i=0,ri=0;;){
        if(r[i]==ri){//na danem miste i je hodnota ze spravneho radku
            i++; re[ri]++;
            if(&r[i] == rs[ri+1]){//jsi na zacatku noveho radku?
                ri++;
                if(ri==rn)//jsi za poslednim radkem?
                    break;
            }
        }else{//na miste i neni hodnota z daneho radku ri
            while(*re[r[i]]==r[i])//dokud je na miste kam chci vlozit aktualni hodnotu spravne cislo
                re[r[i]]++;//posun teor. konec dal
            //prohod hodnoty, koukas na spatne umistene cislo
            swap(c[i],c[c[i]]);
            swap(r[i],r[r[i]]);
            swap(v[i],v[c[i]]);
        }
    }
    delete[] re;
}


// ==============================================================================
//                             SPRAVA PAMETI
// ==============================================================================
// ==============================================================================
/**
    \brief Najde nove dostatecne dlouhe pole, nebo si ho alokuje
    \return hodnoty v a c, jsou to pointery na zacatky volnych poli
    \param l length of required space
    \return  pointing to beginning of new row of value
*/
double* SPARSE::M_get_free_place(const long& l)
{
    bool    found = false;
    double* v; //ukazatel na zacatek nalezeneho radku

    v = M_find_free_place(l,found);
    if (found)
        return v;

    // jestli neni dost mista v poli, alokuj nove pole
    M_alloc_new_field();


    //uz by tam melo byt dost mista
    v = M_find_free_place(l,found);
    if (!found){
        throw("M_get_free_place se zacykloval a nebude dale pokracovat. Je treba prodlouzit delku pole value[>0]");
        /// \warning Muze se zacyklovat, kdyz bude pozadovat delsi radek nez je delka extra_field_length \see extra_field_length
    }
    return (v);
}


// ==============================================================================
/**
    \brief Najde nove dostatecne dlouhe pole, jestli ne, informuje o tom ve found
    \return hodnoty v a c --> jsou to pointery na zacatky volnych poli
    \param l length of required space
    \param v return value, pointing to beginning of new row of value
    \param c return value, pointing to beginning of new column of value
*/
double* SPARSE::M_find_free_place(const long& l, bool& found)
{
    double* v;
    found = false;
    long empt_space_length = M_empty_end[allocated_field] - M_empty_start[allocated_field]+1;

    if(empt_space_length > l){//v alokovanem poli je dost mista
        v = M_empty_start[allocated_field];
        M_update_empty_pointers(l);
        found = true;
    }
    return v;
}


// ==============================================================================
/**
    \brief Zalokuje nove pole, jestli nejni nove misto, tak prealokuje matici
*/
bool SPARSE::M_alloc_new_field()
{

    bool found = false;
    if(allocated_field == (extra_values-1)){
        const int rm = 10; //o kolik bude rozsirena matice
        M_shake_down(rm);
        found = true;
    }else{
        allocated_field++;
        M_allocated_field_length[allocated_field] = extra_field_length;

        value[allocated_field]  = new double [M_allocated_field_length[allocated_field]];
        column[allocated_field] = new long  [M_allocated_field_length[allocated_field]];

        M_empty_start[allocated_field] = value[allocated_field];
        M_empty_end[allocated_field] = &value[allocated_field][M_allocated_field_length[allocated_field]];

        found = true;
    }

    return found;
}


// ==============================================================================
/**
    \brief Setrese pole, vysledek bude jen jedno zalokovane pole.
    \warning Muze zpusobit bad_alloc, protoze alokuje dvakrat tak velkou pamet nez je potreba pro jednu matici
    \warning Dodela spravu pameti, kdyz bude shakedown pouzit funkci LDU nemusi byt prealokovan s plnym poctem radku
    \todo - Zkusit udelat tak aby se neodalokovavaly vsechna pole
    \param rm Rozsireni matice - o kolik bude nulte pole delsi nez je pocet hodnot, v %
    \param new_value Nove pole value, ktere nahradi stare pole value
    \param new_column Nove pole column, ktere nahradi stare pole column
 */
void SPARSE::M_shake_down(const int& rm)
{
    printf("\n\n\n\n\n\n ============================Jsem v Shakedown============================\n\n\n\n\n\n");
    // - Zalokuj nove pole
    long    i,j,ii;
    long     new_field_length   = enlarge_field(value_number,rm);
    double*  new_value          = new double[new_field_length];
    long*    new_column         = new long[new_field_length];
    double** new_value_start    = new double*[row_number];
    long**   new_column_start   = new long*[row_number];
    long*    new_value_length   = new long[row_number];
    double** new_diagonal       = new double*[row_number];

    ii=0;

    // - Zkopiruj do nich hodnoty
    for (i=0; i<row_number; i++){//projed po radcich
        new_value_start[i] = &new_value[ii];
        new_column_start[i] = &new_column[ii];
            new_value_length[i] = value_length[i];

        for (j=0; j<value_length[i]; j++){//projed po sloupcich
            new_value[ii]=value_start[i][j];
            new_column[ii]=column_start[i][j];
            if(new_column[ii]==i){
                new_diagonal[i]=&new_value[ii];
//                dia = true;
            }
            ii++;
        }
    }

    // - Odalokuj puvodne alokovana pole
    dealloc();

    // - Preukazuj ukazatele
    M_aloc          = true;
    value[0]        = new_value;
    column[0]       = new_column;
    value_start     = new_value_start;
    column_start    = new_column_start;
    value_length    = new_value_length;
    diagonal        = new_diagonal;
    value_number    = ii;

    // - Predelej veci na spravu pameti
    allocated_field = 0;
    M_allocated_field_length[0]     = new_field_length;
    M_empty_start[allocated_field]  = value[allocated_field]+ii;
    M_empty_end[allocated_field]    = &value[allocated_field][M_allocated_field_length[0]-1];
}


// ==============================================================================
/**
    \brief Nastavi pointery na zacatky radku.
*/
void SPARSE::initialize_value_length(const long* length)
{
    long i,k;
    k=0;
    for(i=0;i<row_number;i++){
        value_start[i]=&value[0][k];
        column_start[i]=&column[0][k];
        k+=length[i];
    }
    M_empty_start[allocated_field]  = &value[0][k];
}

