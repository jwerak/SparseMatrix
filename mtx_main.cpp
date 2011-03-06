#include "sparse.h"
#include "row.h"
#include "vtr.h"
#include "ctime"

int main()
{
    VTR x;
    MTX* A = new SPARSE;
//    A->loadFromFile("txt_matice/matice.txt");
    A->loadFromFile("txt_matice/matice100x100.dat");
//    A->loadFromFile("txt_matice/test_mtx_25x25.dat");
    //    A->loadFromFile("txt_matice/test_mtx_10x10.dat");
    VTR b(A->get_row_number());
    VTR c(A->get_row_number());//check vektor
//    x.load_from_file("txt_matice/x_vektor.txt");//spocteny vektor matice matice.txt
//    b.load_from_file("txt_matice/b_10.dat");
//    b.load_from_file("txt_matice/b_25.dat");
//    b.print();

    b.ones();
    MTX* LDU = new SPARSE;

    clock_t start, stop;
    double t = 0.0;
    start = clock();

//    try{
        A->solve_Axb((*LDU), x, b);
//    }catch(...){
//        printf("Chyba: Axb \n");
//        throw;
//    }

    stop = clock();
    t = (double) (stop-start)/CLOCKS_PER_SEC;
    t /= 60.0;

    LDU->write_log(t);

//    x.print_to_file("x_vektor.txt");
//    x.print();
//    (*A).print_full_mtx();

    //===============================
    //TESTOVANI SPRAVNEHO X
    //===============================
    c = A->multiply(x);
    //zjisti jestli jsou vsude napsany jednicky
//    printf("\nZacinam overovat spravnost hodnot ve vektoru b:\n");
//    for(long i=0;i<A->get_row_number();i++){
//        if(c.get_element(i) != b.get_element(i))
//            printf("V radku %li vektoru c je hodnota: %lf. \tMa byt:%lf \n",i,c.get_element(i),b.get_element(i));
//    }
    c.print_to_file("c_kontrola.txt");



//    printf("\nMatice pred rozlozenim\n");
//    (*A).print_full_mtx();

//===============================================================
    //ROZKLAD MATICE
//    (*A).make_ldu(*LDU);
//    printf("\nRozlozena matice LDU\n");
//    (*LDU).print_full_mtx();

//===============================================================
//    TESTOVANI MATICE*
//    MTX* L = new SPARSE;
//    MTX* D = new SPARSE;
//    MTX* U = new SPARSE;


//    printf("\nZiskavani casti U radku LDU:\n");
//    for(long i=0;i<9;i++){
//        ROW l = LDU->get_rowU(i);
//        l.print_full();
//    }
//    (*L).get_mtx_L(*LDU);
//    printf("\nSpodni cast rozlozene matice LDU\n");
//    (*L).print_full_mtx();

//    (*D).get_mtx_D(*LDU);
//    printf("\nSpodni cast rozlozene matice LDU\n");
//    (*D).print_full_mtx();

//    (*U).get_mtx_U(*LDU);
//    printf("\nSpodni cast rozlozene matice LDU\n");
//    (*U).print_full_mtx();

//===============================================================
//    //KONTROLA NASOBENI
//    MTX* R = new SPARSE;
//ROZLOZENE
//    (*L).multiply_matrices(*D,*R);
//    (*R).multiply_matrices(*U,*D);
//    printf("\nMatice po roznasobeni\n");
//    (*D).print_full_mtx();
//ORIGINAL
//    (*A).multiply_matrices(*LDU,*R);
//    printf("\nROZNASOBENA MATICE\n");
//    (*R).print_full_mtx();
}

