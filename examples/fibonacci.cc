#include <LinRealRecursiveSeq.h>

#include <cmath>    // sqrt
#include <vector>
#include <fstream>  // ofstream
#include <iostream> // cout
#include <iomanip>  // setprecision

using namespace std;

/*******************************************************************
 * This example compares 92 precomputed Fibonacci numbers:         *
 *                                                                 *
 *     fib(n) = fib(n-1) + fib(n-2)                                *
 *                                                                 *
 * with calculated with LinRealRecursiveSeq class. Both absolute   *
 * and relative errors are computed. Program uses two methods:     *
 * straightforward comparison with the result (double variable)    *
 * and comparison with the result rounded to the nearest long long *
 * number. Average error is printed to the standard output. Errors *
 * for each k are saved in 'fibonacci_no_round.dat' and            *
 * 'fibonacci_round.dat' output files.                             *
 *******************************************************************/

int main(){
    /* First two elements and recurrence relation */
    vector<double> fiboFirstElems, fiboRecurRel;

    fiboFirstElems.push_back(0.0);
    fiboFirstElems.push_back(1.0);

    fiboRecurRel.push_back(1.0);
    fiboRecurRel.push_back(1.0);

    /* Recursive sequence itself */
    LinRealRecursiveSeq fiboSeq(fiboFirstElems,
                                fiboRecurRel);

    /* Variable used to store kth element of the sequence */
    double fiboSeqElem;

    /* Precomputed Fibonacci numbers */
    unsigned long long int fib[]
        = {0,
           1,1,2,3,
           5,8,13,21,
           34,55,89,144,
           233,377,610,987,
           1597,2584,4181,6765,
           10946,17711,28657,46368,
           75025,121393,196418,317811,
           514229,832040,1346269,2178309,
           3524578,5702887,9227465,14930352,
           24157817,39088169,63245986,102334155,
           165580141,267914296,433494437,701408733,
           1134903170,1836311903,2971215073,4807526976,
           7778742049,12586269025,20365011074,32951280099,
           53316291173,86267571272,139583862445,225851433717,
           365435296162,591286729879,956722026041,1548008755920,
           2504730781961,4052739537881,6557470319842,10610209857723,
           17167680177565,27777890035288,44945570212853,72723460248141,
           117669030460994,190392490709135,308061521170129,498454011879264,
           806515533049393,1304969544928657,2111485077978050,3416454622906707,
           5527939700884757,8944394323791464,14472334024676221,23416728348467685,
           37889062373143906,61305790721611591,99194853094755497,160500643816367088,
           259695496911122585,420196140727489673,679891637638612258,1100087778366101931,
           1779979416004714189,2880067194370816120,4660046610375530309,7540113804746346429};

    /* Variables used for calculating and storing average errors and their
     * standard deviations */
    double fibo_NoRound_TmpAbs,          fibo_Round_TmpAbs;
    double fibo_NoRound_TmpRel,          fibo_Round_TmpRel;

    double fibo_NoRound_ErrAbsAv  = 0.0, fibo_Round_ErrAbsAv  = 0.0;
    double fibo_NoRound_ErrRelAv  = 0.0, fibo_Round_ErrRelAv  = 0.0;

    double fibo_NoRound_ErrAbsAv2 = 0.0, fibo_Round_ErrAbsAv2 = 0.0;
    double fibo_NoRound_ErrRelAv2 = 0.0, fibo_Round_ErrRelAv2 = 0.0;

    /* Output files */
    ofstream fiboOutNoRound, fiboOutRound;
    fiboOutNoRound.open("fibonacci_no_round.dat");
    fiboOutRound.open("fibonacci_round.dat");

    cout << "Fibonacci test:" << endl;

    fiboOutNoRound << "#k"
                   << "\t" << "abs. err."
                   << "\t" << "rel. err." << endl;

    fiboOutRound   << "#k"
                   << "\t" << "abs. err."
                   << "\t" << "rel. err." << endl;

    for (unsigned int k = 1; k <= 92; ++k){
        /* kth element */
        fiboSeqElem = fiboSeq.Element(k);

        /***************/
        /* No rounding */
        /***************/

        /* Absolute and relative errors */
        fibo_NoRound_TmpAbs     = fib[k] - fiboSeqElem;
        fibo_NoRound_TmpRel     = 1.0 - fiboSeqElem / fib[k];

        /* Absolute and relative averages of errors */
        fibo_NoRound_ErrAbsAv  += fibo_NoRound_TmpAbs;
        fibo_NoRound_ErrRelAv  += fibo_NoRound_TmpRel;

        /* Absolute and relative averages of squared errors */
        fibo_NoRound_ErrAbsAv2 += fibo_NoRound_TmpAbs * fibo_NoRound_TmpAbs;
        fibo_NoRound_ErrRelAv2 += fibo_NoRound_TmpRel * fibo_NoRound_TmpRel;

        /* Save errors to output */
        fiboOutNoRound << k 
                       << "\t" << fibo_NoRound_TmpAbs
                       << "\t" << fibo_NoRound_TmpRel << endl;

        /************/
        /* Rounding */
        /************/

        /* Absolute and relative errors */
        fibo_Round_TmpAbs     = fib[k] - roundl(fiboSeqElem);
        fibo_Round_TmpRel     = 1.0 - roundl(fiboSeqElem) / static_cast<double>(fib[k]);

        /* Absolute and relative averages of errors */
        fibo_Round_ErrAbsAv  += fibo_Round_TmpAbs;
        fibo_Round_ErrRelAv  += fibo_Round_TmpRel;

        /* Absolute and relative averages of squared errors */
        fibo_Round_ErrAbsAv2 += fibo_Round_TmpAbs * fibo_Round_TmpAbs;
        fibo_Round_ErrRelAv2 += fibo_Round_TmpRel * fibo_Round_TmpRel;

        /* Save errors to output */
        fiboOutRound   << k 
                       << "\t" << fibo_Round_TmpAbs
                       << "\t" << fibo_Round_TmpRel << endl;
    }

    fiboOutNoRound.close();
    fiboOutRound.close();

    /* Divide all averages by number of compared elements */
    fibo_NoRound_ErrAbsAv  /= 91.0;
    fibo_NoRound_ErrAbsAv2 /= 91.0;
    fibo_NoRound_ErrRelAv  /= 91.0;
    fibo_NoRound_ErrRelAv2 /= 91.0;

    fibo_Round_ErrAbsAv    /= 91.0;
    fibo_Round_ErrAbsAv2   /= 91.0;
    fibo_Round_ErrRelAv    /= 91.0;
    fibo_Round_ErrRelAv2   /= 91.0;

    cout << setprecision(4);

    /* Print average errors to the standard output */
    cout << "    absolute error (no rounding): "
         << fibo_NoRound_ErrAbsAv
         << " +/- "
         << sqrt(fibo_NoRound_ErrAbsAv2 - fibo_NoRound_ErrAbsAv * fibo_NoRound_ErrAbsAv) << endl;

    cout << "    relative error (no rounding): "
         << fibo_NoRound_ErrRelAv
         << " +/- "
         << sqrt(fibo_NoRound_ErrRelAv2 - fibo_NoRound_ErrRelAv * fibo_NoRound_ErrRelAv) << endl;

    cout << "    absolute error (rounding):    "
         << fibo_Round_ErrAbsAv
         << " +/- "
         << sqrt(fibo_Round_ErrAbsAv2 - fibo_Round_ErrAbsAv * fibo_Round_ErrAbsAv) << endl;

    cout << "    relative error (rounding):    "
         << fibo_Round_ErrRelAv
         << " +/- "
         << sqrt(fibo_Round_ErrRelAv2 - fibo_Round_ErrRelAv * fibo_Round_ErrRelAv) << endl;

    return 0;
}
