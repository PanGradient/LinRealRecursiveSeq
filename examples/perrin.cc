#include <LinRealRecursiveSeq.h>

#include <cmath>    // sqrt
#include <vector>
#include <fstream>  // ofstream
#include <iostream> // cout
#include <iomanip>  // setprecision

using namespace std;

/*******************************************************************
 * This example compares 153 precomputed Perrin numbers:           *
 *                                                                 *
 *     perrin(n) = perrin(n-2) + perrin(n-3)                       *
 *                                                                 *
 * with calculated with LinRealRecursiveSeq class. Both absolute   *
 * and relative errors are computed. Program uses two methods:     *
 * straightforward comparison with the result (double variable)    *
 * and comparison with the result rounded to the nearest long long *
 * number. Average error is printed to the standard output. Errors *
 * for each k are saved in 'perrinn_no_round.dat' and              *
 * 'perrinn_round.dat' output files.                               *
 *******************************************************************/

int main(){
    /* First three elements and recurrence relation */
    vector<double> perrinFirstElems, perrinRecurRel;

    perrinFirstElems.push_back(3.0);
    perrinFirstElems.push_back(0.0);
    perrinFirstElems.push_back(2.0);

    perrinRecurRel.push_back(1.0);
    perrinRecurRel.push_back(1.0);
    perrinRecurRel.push_back(0.0);

    /* Recursive sequence itself */
    LinRealRecursiveSeq perrinSeq(perrinFirstElems,
                                  perrinRecurRel);

    /* Variable used to store kth element of the sequence */
    double perrinSeqElem;

    /* Precomputed Perrin numbers */
    unsigned long long int perrin[]
        = {3,0,2,3,2,5,5,7,
           10,12,17,22,29,39,
           51,68,90,119,158,209,
           277,367,486,644,853,1130,
           1497,1983,2627,3480,4610,6107,
           8090,10717,14197,18807,24914,33004,
           43721,57918,76725,101639,134643,178364,
           236282,313007,414646,549289,727653,963935,
           1276942,1691588,2240877,2968530,3932465,5209407,
           6900995,9141872,12110402,16042867,21252274,28153269,
           37295141,49405543,65448410,86700684,114853953,152149094,
           201554637,267003047,353703731,468557684,620706778,822261415,
           1089264462,1442968193,1911525877,2532232655,3354494070,4443758532,
           5886726725,7798252602,10330485257,13684979327,18128737859,24015464584,
           31813717186,42144202443,55829181770,73957919629,97973384213,129787101399,
           171931303842,227760485612,301718405241,399691789454,529478890853,701410194695,
           929170680307,1230889085548,1630580875002,2160059765855,2861469960550,3790640640857,
           5021529726405,6652110601407,8812170367262,11673640327812,15464280968669,20485810695074,
           27137921296481,35950091663743,47623731991555,63088012960224,83573823655298,110711744951779,
           146661836615522,194285568607077,257373581567301,340947405222599,451659150174378,598320986789900,
           792606555396977,1049980136964278,1390927542186877,1842586692361255,2440907679151155,3233514234548132,
           4283494371512410,5674421913699287,7517008606060542,9957916285211697,13191430519759829,17474924891272239,
           23149346804971526,30666355411032068,40624271696243765,53815702216003594,71290627107275833,94439973912247359,
           125106329323279427,165730601019523192,219546303235526786,290836930342802619,385276904255049978,510383233578329405,
           676113834597852597,895660137833379383,1186497068176182002,1571773972431231980,2082157206009561385,2758271040607413982,
           3653931178440793365,4840428246616975367,6412202219048207347,8494359425057768732};

    /* Variables used for calculating and storing average errors and their
     * standard deviations */
    double perrin_NoRound_TmpAbs,          perrin_Round_TmpAbs;
    double perrin_NoRound_TmpRel,          perrin_Round_TmpRel;

    double perrin_NoRound_ErrAbsAv  = 0.0, perrin_Round_ErrAbsAv  = 0.0;
    double perrin_NoRound_ErrRelAv  = 0.0, perrin_Round_ErrRelAv  = 0.0;

    double perrin_NoRound_ErrAbsAv2 = 0.0, perrin_Round_ErrAbsAv2 = 0.0;
    double perrin_NoRound_ErrRelAv2 = 0.0, perrin_Round_ErrRelAv2 = 0.0;

    /* Output files */
    ofstream perrinOutNoRound, perrinOutRound;
    perrinOutNoRound.open("perrin_no_round.dat");
    perrinOutRound.open("perrin_round.dat");

    cout << "Perrin test:" << endl;

    perrinOutNoRound << "#k"
                     << "\t" << "abs. err."
                     << "\t" << "rel. err." << endl;

    perrinOutRound   << "#k"
                     << "\t" << "abs. err."
                     << "\t" << "rel. err." << endl;

    for (unsigned int k = 3; k < 156; ++k){
        /* kth element */
        perrinSeqElem = perrinSeq.Element(k);

        /***************/
        /* No rounding */
        /***************/

        /* Absolute and relative errors */
        perrin_NoRound_TmpAbs     = perrin[k] - perrinSeqElem;
        perrin_NoRound_TmpRel     = 1.0 - perrinSeqElem / perrin[k];

        /* Absolute and relative averages of errors */
        perrin_NoRound_ErrAbsAv  += perrin_NoRound_TmpAbs;
        perrin_NoRound_ErrRelAv  += perrin_NoRound_TmpRel;

        /* Absolute and relative averages of squared errors */
        perrin_NoRound_ErrAbsAv2 += perrin_NoRound_TmpAbs * perrin_NoRound_TmpAbs;
        perrin_NoRound_ErrRelAv2 += perrin_NoRound_TmpRel * perrin_NoRound_TmpRel;

        /* Save errors to output */
        perrinOutNoRound << k 
                       << "\t" << perrin_NoRound_TmpAbs
                       << "\t" << perrin_NoRound_TmpRel << endl;

        /************/
        /* Rounding */
        /************/

        /* Absolute and relative errors */
        perrin_Round_TmpAbs     = perrin[k] - roundl(perrinSeqElem);
        perrin_Round_TmpRel     = 1.0 - roundl(perrinSeqElem) / static_cast<double>(perrin[k]);

        /* Absolute and relative averages of errors */
        perrin_Round_ErrAbsAv  += perrin_Round_TmpAbs;
        perrin_Round_ErrRelAv  += perrin_Round_TmpRel;

        /* Absolute and relative averages of squared errors */
        perrin_Round_ErrAbsAv2 += perrin_Round_TmpAbs * perrin_Round_TmpAbs;
        perrin_Round_ErrRelAv2 += perrin_Round_TmpRel * perrin_Round_TmpRel;

        /* Save errors to output */
        perrinOutRound   << k 
                       << "\t" << perrin_Round_TmpAbs
                       << "\t" << perrin_Round_TmpRel << endl;
    }

    perrinOutNoRound.close();
    perrinOutRound.close();

    /* Divide all averages by number of compared elements */
    perrin_NoRound_ErrAbsAv  /= 153.0;
    perrin_NoRound_ErrAbsAv2 /= 153.0;
    perrin_NoRound_ErrRelAv  /= 153.0;
    perrin_NoRound_ErrRelAv2 /= 153.0;

    perrin_Round_ErrAbsAv    /= 153.0;
    perrin_Round_ErrAbsAv2   /= 153.0;
    perrin_Round_ErrRelAv    /= 153.0;
    perrin_Round_ErrRelAv2   /= 153.0;

    cout << setprecision(4);

    /* Print average errors to the standard output */
    cout << "    absolute error (no rounding): "
         << perrin_NoRound_ErrAbsAv
         << " +/- "
         << sqrt(perrin_NoRound_ErrAbsAv2 - perrin_NoRound_ErrAbsAv * perrin_NoRound_ErrAbsAv) << endl;

    cout << "    relative error (no rounding): "
         << perrin_NoRound_ErrRelAv
         << " +/- "
         << sqrt(perrin_NoRound_ErrRelAv2 - perrin_NoRound_ErrRelAv * perrin_NoRound_ErrRelAv) << endl;

    cout << "    absolute error (rounding):    "
         << perrin_Round_ErrAbsAv
         << " +/- "
         << sqrt(perrin_Round_ErrAbsAv2 - perrin_Round_ErrAbsAv * perrin_Round_ErrAbsAv) << endl;

    cout << "    relative error (rounding):    "
         << perrin_Round_ErrRelAv
         << " +/- "
         << sqrt(perrin_Round_ErrRelAv2 - perrin_Round_ErrRelAv * perrin_Round_ErrRelAv) << endl;

    return 0;
}
