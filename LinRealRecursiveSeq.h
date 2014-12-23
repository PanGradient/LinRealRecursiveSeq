#ifndef LIN_REAL_RECURSIVE_SEQ_H
#define LIN_REAL_RECURSIVE_SEQ_H

#include <cstdlib>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*
 * Linear Real Recursive Sequence with constant coefficients class
 * Allows quick (although limited by the machine precision) calculation of kth element of
 * a recursive sequence S defined by recurrence relation:
 *
 *     S(k) = a(0) * S(k-N) + ... a(N-1) * S(k-1)
 *
 * with a(m) being constant coefficients. It uses eigenvalue decomposition of recurrence
 * relation matrix representation:
 *
 *     [ S(k)     ]   [ 0    1    0    ... 0      0      ]   [ S(k-1)   ]
 *     [ S(k-1)   ]   [ 0    0    1    ... 0      0      ]   [ S(k-2)   ]
 *     [   .      ]   [ .    .    .    .   .      .      ]   [   .      ]
 *     [   .      ] = [ .    .    .     .  .      .      ] * [   .      ]
 *     [   .      ]   [ .    .    .      . .      .      ]   [   .      ]
 *     [ S(k-N+2) ]   [ 0    0    0    ... 0      1      ]   [ S(k-N+1) ]
 *     [ S(k-N+1) ]   [ a(0) a(1) a(2) ... a(N-2) a(N-1) ]   [ S(k-N)   ]
 *
 * Decomposing above matrix in a form C * diag( { m(0), ..., m(N-1) } ) * C^(-1) (m(i) being
 * its complex eigenvalues and C being a complex matrix build from eigenvectors) allows for
 * quick calculation of kth element as 0th element of vector:
 *
 *     v = C * diag( { m(0)^k, ..., m(N-1)^k } ) * C^(-1) * v0
 *
 * with v0 begin vector composed of the first elements of the sequence.
 * All calculations are performed by using GSL complex vectors and matrices and GSL BLAS
 * Interface.
 */
class LinRealRecursiveSeq
{
    public:
        /*
         * Default constructor
         * Makes a constant sequence of zeroes.
         */
        LinRealRecursiveSeq();

        /*
         * Copy constructor and operator
         */
        LinRealRecursiveSeq (const LinRealRecursiveSeq& realRecurSeq);
        LinRealRecursiveSeq& operator= (const LinRealRecursiveSeq& realRecurSeq);

        /*
         * Constructor
         * Takes as an input stdlib vector of doubles with first N elements of recursive
         * sequence (firstElements) and stdlib vector of doubles defining the recurrence
         * relation (recurrenceRelation) given by:
         *
         *     S(n) = recurrenceRelation[0] * S(n-N)
         *          + ...
         *          + recurrenceRelation[N-1] * S(n-1)
         *
         * Both have to have the same size.
         */
        LinRealRecursiveSeq (const std::vector<double>& firstElements,
                             const std::vector<double>& recurrenceRelation);

        /*
         * Constructor
         * Behaves the same way as constructor using stdlib vectors, only this one uses arrays of
         * doubles and their size N as an input. Since implementation would be (nearly) the same,
         * it uses constructor delegation allowed by C++11 standard (makes it necessary to use
         * compiler compliant with this standard).
         */
        LinRealRecursiveSeq (const double *firstElements,
                             const double *recurrenceRelation,
                             size_t N) : LinRealRecursiveSeq(std::vector<double>(firstElements,
                                                                                 firstElements + N),
                                                             std::vector<double>(recurrenceRelation,
                                                                                 recurrenceRelation + N)) { }

        /*
         * Destructor
         * Deallocation is performed, where it's necessary.
         */
       ~LinRealRecursiveSeq ();

        /*
         * Calculates kth element of recursive sequence using the decomposition of recurrence
         * relation matrix form M performed by class constructor:
         *     M = C * diag( {m1, ..., mN} ) * C^(-1) =>
         *     => M^k = C * diag( {m1^k, ..., mN^k} ) * C^(-1)
         * and multiplying those matrices by first elements vector. The kth element is the
         * first element of the resulting vector. Because of rounding errors, computed vector
         * might have some small imaginary part, which is neglected.
         */
        double Element(unsigned int k);


    private:
        /*
         * Complex GSL vector containing first N elements of recursive sequence.
         */
        gsl_vector_complex *_first_elems_;

        /*
         * Complex GSL vector containing eigenvalues of recurrence relation matrix
         * representation.
         */
        gsl_vector_complex *_eigenvals_;

        /*
         * Complex GSL transition matrix and its inverse of recurrence relation matrix
         * representation.
         */
        gsl_matrix_complex *_transM_, *_invTransM_;

        /*
         * Complex diagonal matrix and complex vectors used in calculation of kth element of
         * sequence. Each evaluation requires them in intermediate calculations. Because of
         * that, they are allocated by the constructor, rather than allocated and deallocated
         * every time Element() function is called.
         */
        gsl_matrix_complex *_diag_;
        gsl_vector_complex *_elems1_, *_elems2_;
        /*
         * Size value for all vector (size _N_) and matrix (size _N_ x _N_) class members.
         */
        size_t _N_;

}; /* class LinRealRecursiveSeq */

#endif /* LIN_REAL_RECURSIVE_SEQ_H */
