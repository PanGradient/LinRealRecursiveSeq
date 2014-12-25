#include <LinRealRecursiveSeq.h>

#include <stdexcept>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>


/*
 * Default constructor
 * Creates a constant sequence of zeroes.
 */
LinRealRecursiveSeq::
  LinRealRecursiveSeq()
{
    /* Set size to 1 */
    _N_ = 1;

    /* Allocate and set first element and eigenvalue to 0.0 + 0.0I */
    _first_elems_ = gsl_vector_complex_alloc (1);
    _eigenvals_   = gsl_vector_complex_alloc (1);

    gsl_complex z;
    GSL_SET_COMPLEX(&z, 0.0, 0.0);
    gsl_vector_complex_set( _first_elems_, 0, z );
    gsl_vector_complex_set( _eigenvals_, 0, z );

    /* Because N = 1, transition matrix, its inverse and memebers used (in general) in
     * intermediate calculations when Element() function is called won't
     * be use in nth sequence element calculation. Therefore, their pointers
     * points to NULL. */
    _transM_    = NULL;
    _invTransM_ = NULL;

    _diag_      = NULL;
    _elems1_    = NULL;
    _elems2_    = NULL;
}


/*
 * Copy constructor
 */
LinRealRecursiveSeq::
  LinRealRecursiveSeq (const LinRealRecursiveSeq& realRecurSeq)
{
    /* Copy size */
    _N_ = realRecurSeq._N_;

    /* Allocate and copy first elements and eigenvalues. */
    _first_elems_ = gsl_vector_complex_alloc (_N_);
    _eigenvals_   = gsl_vector_complex_alloc (_N_);

    gsl_vector_complex_memcpy ( _first_elems_, realRecurSeq._first_elems_ );
    gsl_vector_complex_memcpy ( _eigenvals_,   realRecurSeq._eigenvals_   );

    /* If the recurrence relation is defined as
     *     S(n) = a*S(n-1)
     * for some a, it's unnecessary to calculate eigenvalues (since
     * there is only one equal to a) and eigenvectors (since there
     * is only one equal to (1.0 + 0.0*I)). Because of this, for N = 1
     * transition matrix, its inverse and memebers used (in general) in
     * intermediate calculations when Element() function is called won't
     * be use in nth sequence element calculation. Therefore, their pointers
     * points to NULL. */
    if (_N_ == 1){
        _transM_    = NULL;
        _invTransM_ = NULL;

        _diag_      = NULL;
        _elems1_    = NULL;
        _elems2_    = NULL;
    }

    /* Otherwise, allocate memory for those members. Off-diagonal elements of
     * diagonal matrix have to be initialized to 0.0 + 0.0I, therefore calloc
     * method is used.*/
    else {
        _transM_    = gsl_matrix_complex_alloc  (_N_, _N_);
        _invTransM_ = gsl_matrix_complex_alloc  (_N_, _N_);

        _diag_      = gsl_matrix_complex_calloc (_N_, _N_);
        _elems1_    = gsl_vector_complex_alloc  (_N_);
        _elems2_    = gsl_vector_complex_alloc  (_N_);

        gsl_matrix_complex_memcpy ( _transM_,    realRecurSeq._transM_    );
        gsl_matrix_complex_memcpy ( _invTransM_, realRecurSeq._invTransM_ );
    }
}


/*
 * Copy operator
 */
LinRealRecursiveSeq&
LinRealRecursiveSeq::
  operator= (const LinRealRecursiveSeq& realRecurSeq)
{
    /* Avoid self-assignment */
    if (this != &realRecurSeq){

        /* Copy size */
        this->_N_ = realRecurSeq._N_;

        /* Reallocate and copy first elements and eigenvalues. */
        gsl_vector_complex_free(this->_first_elems_);
        gsl_vector_complex_free(this->_eigenvals_);

        this->_first_elems_ = gsl_vector_complex_alloc (_N_);
        this->_eigenvals_   = gsl_vector_complex_alloc (_N_);

        gsl_vector_complex_memcpy ( this->_first_elems_, realRecurSeq._first_elems_ );
        gsl_vector_complex_memcpy ( this->_eigenvals_,   realRecurSeq._eigenvals_   );

        /* If not null pointers, deallocate other vector and metrix members. */
        if ( this->_transM_    != NULL )  gsl_matrix_complex_free(this->_transM_);
        if ( this->_invTransM_ != NULL )  gsl_matrix_complex_free(this->_invTransM_);
        if ( this->_diag_      != NULL )  gsl_matrix_complex_free(this->_diag_);
        if ( this->_elems1_    != NULL )  gsl_vector_complex_free(this->_elems1_);
        if ( this->_elems2_    != NULL )  gsl_vector_complex_free(this->_elems2_);

        /* If the recurrence relation is defined as
         *     S(n) = a*S(n-1)
         * for some a, it's unnecessary to calculate eigenvalues (since
         * there is only one equal to a) and eigenvectors (since there
         * is only one equal to (1.0 + 0.0*I)). Because of this, for N = 1
         * transition matrix, its inverse and memebers used (in general) in
         * intermediate calculations when Element() function is called won't
         * be use in nth sequence element calculation. Therefore, their pointers
         * points to NULL. */
        if (_N_ == 1){
            this->_transM_    = NULL;
            this->_invTransM_ = NULL;

            this->_diag_      = NULL;
            this->_elems1_    = NULL;
            this->_elems2_    = NULL;
        }

        /* Otherwise, allocate memory for those members. Off-diagonal elements of
         * diagonal matrix have to be initialized to 0.0 + 0.0I, therefore calloc
         * method is used.*/
        else {
            this->_transM_    = gsl_matrix_complex_alloc  (_N_, _N_);
            this->_invTransM_ = gsl_matrix_complex_alloc  (_N_, _N_);

            this->_diag_      = gsl_matrix_complex_calloc (_N_, _N_);
            this->_elems1_    = gsl_vector_complex_alloc  (_N_);
            this->_elems2_    = gsl_vector_complex_alloc  (_N_);

            gsl_matrix_complex_memcpy ( this->_transM_,    realRecurSeq._transM_    );
            gsl_matrix_complex_memcpy ( this->_invTransM_, realRecurSeq._invTransM_ );
        }
    }

    return *this;
}


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
LinRealRecursiveSeq::
  LinRealRecursiveSeq (const std::vector<double>& firstElements,
                       const std::vector<double>& recurrenceRelation)
{
    /* Local complex variable, used throughout the constructor */
    gsl_complex z;

    /* Assuming nth sequence element is dependent on at least (n-1)th
     * element, otherwise throw lenghth exception.
     * Note that it's not the same as saying its dependent on one of
     * the previous elements. For instance if:
     *     S(n) = 2.0 * S(n-3)
     * requires recurrenceRelation to be { 1.0, 0.0, 0.0 }. */
    if (firstElements.size() == 0 || recurrenceRelation.size() == 0){
        throw std::length_error("declared matrix size is 0");
    }

    if (firstElements.size() != recurrenceRelation.size()){
        throw std::length_error("number of the first elements not equal to number of recurrence relation coefficients");
    }

    _N_ = firstElements.size();

    /* If the recurrence relation is defined as
     *     S(n) = a*S(n-1)
     * for some a, it's unnecessary to calculate eigenvalues (since
     * there is only one equal to a) and eigenvectors (since there
     * is only one equal to (1.0 + 0.0*I)). Because of this, for N = 1
     * transition matrix, its inverse and memebers used (in general) in
     * intermediate calculations when Element() function is called won't
     * be use in nth sequence element calculation. Therefore, their pointers
     * points to NULL. */
    if (_N_ == 1){
        _first_elems_ = gsl_vector_complex_alloc(1);
        _eigenvals_   = gsl_vector_complex_alloc(1);

        GSL_SET_COMPLEX( &z, firstElements[0], 0.0 );
        gsl_vector_complex_set( _first_elems_, 0, z );

        GSL_SET_COMPLEX( &z, recurrenceRelation[0], 0.0 );
        gsl_vector_complex_set( _eigenvals_, 0, z );

        _transM_    = NULL;
        _invTransM_ = NULL;

        _diag_      = NULL;
        _elems1_    = NULL;
        _elems2_    = NULL;
    }

    /* Otherwise, full eigensystem calculation is required. Since matrix
     * expressing the recurrence relation isn't in general symmetric, its
     * eigenvalues may be complex. This is the reason GSL complex vectors
     * and matrices are used. */
    else {
        size_t i;

        /* Members allocation. Off-diagonal elements of diagonal matrix have
         * to be initialized to 0.0 + 0.0I, therefore calloc method is
         * used. */
        _first_elems_ = gsl_vector_complex_alloc  (_N_);
        _eigenvals_   = gsl_vector_complex_alloc  (_N_);
        _transM_      = gsl_matrix_complex_alloc  (_N_, _N_);
        _invTransM_   = gsl_matrix_complex_alloc  (_N_, _N_);

        _diag_        = gsl_matrix_complex_calloc (_N_, _N_);
        _elems1_      = gsl_vector_complex_alloc  (_N_);
        _elems2_      = gsl_vector_complex_alloc  (_N_);

        /* Assigning first N sequence elements. */
        for (i = 0; i < _N_; ++i){
            GSL_SET_COMPLEX ( &z, firstElements[i], 0.0 );
            gsl_vector_complex_set ( _first_elems_, i, z );
        }

        /* Matrix form of recurrence sequence:
         *         [ 0  1  0  ... 0      0      ]
         *         [ 0  0  1  ... 0      0      ]
         *         [ .  .  .  .   .      .      ]
         *     M = [ .  .  .   .  .      .      ]
         *         [ .  .  .    . .      .      ]
         *         [ 0  0  0  ... 0      1      ]
         *         [ r0 r1 r2 ... r(n-2) r(n-1) ]
         * with rk = recurrenceRelation[k] array. Diagonal eigenvalue matrix,
         * transition matrix and its inverse are calculated from decomposition
         *     M = C * diag * C^(-1) =>
         *          => diag   = _eigenvalsM_
         *             C      = _transM_
         *             C^(-1) = _invTransM_
         * This can be obtained by calculating eigenvalues and eigenvectors of M
         * using GSL methods.
         */
        gsl_matrix *M = gsl_matrix_calloc(_N_, _N_);

        for (i = 0; i < _N_-1; ++i){
            gsl_matrix_set ( M, i, i+1, 1.0 );
            gsl_matrix_set ( M, _N_-1, i, recurrenceRelation[i] );
        }
        gsl_matrix_set ( M, _N_-1, _N_-1, recurrenceRelation[_N_-1] );

        /* GSL workspace for eigensystem calculations. */
        gsl_eigen_nonsymmv_workspace * work = gsl_eigen_nonsymmv_alloc(_N_);
        gsl_eigen_nonsymmv ( M, _eigenvals_, _transM_, work );

        /* Workspace and M matrix aren't needed anymore, deallocation performed. */
        gsl_eigen_nonsymmv_free(work);
        gsl_matrix_free(M);

        /* Inverting transition matrix. GSL methods requires previously preformed
         * LU decomposition for inverting matrices. For this a dummy int variable is
         * used and a GSL permutation object. After the matrix inversion, permutation
         * is deallocated. */
        int dummy;
        gsl_permutation *perm = gsl_permutation_alloc(_N_);
        gsl_linalg_complex_LU_decomp ( _transM_, perm, &dummy );
        gsl_linalg_complex_LU_invert ( _transM_, perm, _invTransM_ );
        gsl_permutation_free(perm);
    }
}


/*
 * Destructor
 * Deallocation is performed, where it's necessary.
 */
LinRealRecursiveSeq::
  ~LinRealRecursiveSeq()
{
    gsl_vector_complex_free(_first_elems_);
    gsl_vector_complex_free(_eigenvals_);

    /* Transition matrix might be a null pointer, deallocate only if it isn't. */
    if (_transM_ != NULL){
        gsl_matrix_complex_free(_transM_);
    }

    /* Inverse of transition matrix might be a null pointer, deallocate only if it
     * isn't. */
    if (_invTransM_ != NULL){
        gsl_matrix_complex_free(_invTransM_);
    }

    /* Diagonal matrix might be a null pointer, deallocate only if it
     * isn't. */
    if (_diag_ != NULL){
        gsl_matrix_complex_free(_diag_);
    }

    /* First complex vector, used in Element() function, might be a null pointer,
     * deallocate only if it isn't. */
    if (_elems1_ != NULL){
        gsl_vector_complex_free(_elems1_);
    }

    /* Second complex vector, used in Element() function, might be a null pointer,
     * deallocate only if it isn't. */
    if (_elems2_ != NULL){
        gsl_vector_complex_free(_elems2_);
    }
}


/*
 * Calculates kth element of recursive sequence using the decomposition of recurrence
 * relation matrix form M performed by class constructor:
 *     M = C * diag( {m1, ..., mN} ) * C^(-1) =>
 *     => M^k = C * diag( {m1^k, ..., mN^k} ) * C^(-1)
 * and multiplying those matrices by first elements vector. The kth element is the
 * first element of the resulting vector. Because of rounding errors, computed vector
 * might have some small imaginary part, which is neglected.
 */
double
LinRealRecursiveSeq::
  Element(unsigned int k)
{
    gsl_complex z;

    /* If kth element is part of first sequence elements, calculations are unnecessary,
     * return already stored value. */
    if (k < _N_){
        z = gsl_vector_complex_get(_first_elems_, k);
        return GSL_REAL(z);
    }

    /* If the recurrence relation is defined as
     *     S(n) = a*S(n-1)
     * the kth element is simply a^k * S(0). */
    else if (_N_ == 1){
        z = gsl_vector_complex_get(_eigenvals_, 0);
        z = gsl_complex_pow_real(z, k);
        return GSL_REAL(z) * GSL_REAL( gsl_vector_complex_get(_first_elems_, 0) );
    }

    /* Otherwise... */
    else {
        /* ... fill diagonal matrix with m_i^k, with m_i being stored eigenvalues. */
        for (unsigned int i = 0; i < _N_; ++i){
            z = gsl_vector_complex_get(_eigenvals_, i);
            gsl_matrix_complex_set ( _diag_, i, i,
                                     gsl_complex_pow_real(z, k) );
        }

        /* Complex constants, used below. */
        gsl_complex unit_complex = gsl_complex_rect(1.0, 0.0);
        gsl_complex zero_complex = gsl_complex_rect(0.0, 0.0);

        /* Matrix-vector intermediate products. Multiplication is provided by GSL BLAS
         * Interface member gsl_blas_zgemv(), which computes:
         *     y = a*op(M)*x + b*y
         * with a,b being complex variables, x some complex vector, y complex vector
         * modified by the function, M - a given matrix and op() operator being chosen
         * matrix transformation (none, transposition or conjugate transposition).
         * Since this function requires y = M*x calculation, the chosen parameters are:
         *     a    = 1.0 + 0.0I
         *     b    = 0.0 + 0.0I
         *     op() = none
         * and _elems1_ and _elems2_ are used, so this function doesn't change value of
         * _first_elems_ and doesn't give undefined behaviour. */
        gsl_blas_zgemv( CblasNoTrans,
                        unit_complex, _invTransM_, _first_elems_,
                        zero_complex, _elems1_ );
        gsl_blas_zgemv( CblasNoTrans,
                        unit_complex, _diag_,   _elems1_,
                        zero_complex, _elems2_ );
        gsl_blas_zgemv( CblasNoTrans,
                        unit_complex, _transM_, _elems2_,
                        zero_complex, _elems1_ );

        /* The kth elements is the first element of calculated complex vector. Because
         * of rounding errors, imaginary number can be non-zero, but insignificant,
         * therefore it's neglected. */
        z = gsl_vector_complex_get(_elems1_, 0);

        return GSL_REAL(z);
    }
}
