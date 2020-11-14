//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_SPACE_H_INCLUDED )
#define  KRATOS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <numeric>




// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/matrix_market_interface.h"
#include "utilities/dof_updater.h"
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

template <class TDataType, class TMatrixType, class TVectorType>
class KratosSpace;

typedef std::size_t IndexType;

template <class TDataType>
using TKratosSmpSparseSpace = KratosSpace<TDataType, CsrMatrix<TDataType,IndexType> , SystemVector<TDataType,IndexType>>;

template <class TDataType>
using TKratosSmpDenseSpace  =  KratosSpace<TDataType, DenseMatrix<TDataType>, DenseVector<TDataType>>;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.

/** Detail class definition.
 */
template<class TDataType, class TMatrixType, class TVectorType>
class KratosSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosSpace
    KRATOS_CLASS_POINTER_DEFINITION(KratosSpace);

    typedef TDataType DataType;

    typedef TMatrixType MatrixType;

    typedef TVectorType VectorType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef typename Kratos::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename Kratos::shared_ptr< TVectorType > VectorPointerType;

    typedef DofUpdater< KratosSpace<TDataType,TMatrixType,TVectorType> > DofUpdaterType;
    typedef typename DofUpdaterType::UniquePointer DofUpdaterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    KratosSpace()
    {
    }

    /// Destructor.

    virtual ~KratosSpace()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(new TMatrixType(0, 0));
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(new TVectorType(0));
    }

    /// return size of vector rV

    static IndexType Size(VectorType const& rV)
    {
        return rV.size();
    }

    /// return number of rows of rM

    static IndexType Size1(MatrixType const& rM)
    {
        return rM.size1();
    }

    /// return number of columns of rM

    static IndexType Size2(MatrixType const& rM)
    {
        return rM.size2();
    }

    /// rXi = rMij
	// This version is needed in order to take one column of multi column solve from AMatrix matrix and pass it to an ublas vector
	template<typename TColumnType>
	static void GetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
	{
		if (rX.size() != rM.size1())
			rX.resize(rM.size1(), false);

		for (IndexType i = 0; i < rM.size1(); i++) {
			rX[i] = rM(i, j);
		}
	}

	// This version is needed in order to take one column of multi column solve from AMatrix matrix and pass it to an ublas vector
	template<typename TColumnType>
	static void SetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
	{
		for (IndexType i = 0; i < rM.size1(); i++) {
			rM(i,j) = rX[i];
		}
	}


    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY = rX;
    }

    /// rY = rX

    static void Copy(VectorType const& rX, VectorType& rY)
    {
        rY = rX;
    }

    /// rX * rY

    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
        return rX.inner_prod(rY);
    }


    /// ||rX||2

    static TDataType TwoNorm(VectorType const& rX)
    {
        return std::sqrt(rX.inner_prod(rX));
    }

    static TDataType TwoNorm(const Matrix& rA) // Frobenious norm
    {
        TDataType aux_sum = TDataType();
        #pragma omp parallel for reduction(+:aux_sum)
        for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
            for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                aux_sum += std::pow(rA(i,j),2);
            }
        }
        return std::sqrt(aux_sum);
    }

    static TDataType TwoNorm(const CsrMatrix<TDataType> & rA) // Frobenious norm
    {
        TDataType sum2 = IndexPartition<IndexType>(rA.value_data().size())
            .for_each<SumReduction<TDataType>>( [&](IndexType i){
                return std::pow(rA.value_data()[i],2);
            });
        return std::sqrt(sum2);

        //TODO implement
        //return rA.FrobeniusNorm();
    }

    /**
     * This method computes the Jacobi norm
     * @param rA The matrix to compute the Jacobi norm
     * @return aux_sum: The Jacobi norm
     */
    static TDataType JacobiNorm(const Matrix& rA)
    {
        TDataType aux_sum = TDataType();
        #pragma omp parallel for reduction(+:aux_sum)
        for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
            for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                if (i != j) {
                    aux_sum += std::abs(rA(i,j));
                }
            }
        }
        return aux_sum;
    }

    static TDataType JacobiNorm(const CsrMatrix<TDataType>& rA)
    {
        TDataType aux_sum = TDataType();
        IndexPartition<IndexType>(rA.size1()).for_each( [&](IndexType i){
            IndexType row_begin = rA.index1_data()[i];
            IndexType row_end   = rA.index1_data()[i+1];
            for(IndexType k = row_begin; k < row_end; ++k){
                IndexType col = rA.index2_data()[k];
                if(i != col)
                    AtomicAdd(aux_sum, std::abs(rA.value_data()[k]));
            }  
        });
        return aux_sum;
    }

    static void Mult(const Matrix& rA, VectorType& rX, VectorType& rY)
    {
        axpy_prod(rA, rX, rY, true);
    }

    static void Mult(const CsrMatrix<TDataType>& rA, VectorType& rX, VectorType& rY)
    {
        // rY = rA * rX
        rA.SpMV(rY,rX);
    }

    template< class TOtherMatrixType >
    static void TransposeMult(TOtherMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        // rY = rA^t * rX
		rA.TransposeSpMV(rY,rX);
    } 

    static inline SizeType GraphDegree(IndexType i, TMatrixType& A)
    {
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        return ( std::distance(a_iterator.begin(), a_iterator.end()));
#else
        return ( std::distance(begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
                               end(a_iterator, boost::numeric::ublas::iterator1_tag())));
#endif
    }

    static inline void GraphNeighbors(IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
    {
        neighbors.clear();
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        for (typename MatrixType::iterator2 row_iterator = a_iterator.begin();
                row_iterator != a_iterator.end(); ++row_iterator)
        {
#else
        for (typename MatrixType::iterator2 row_iterator = begin(a_iterator,
                boost::numeric::ublas::iterator1_tag());
                row_iterator != end(a_iterator,
                                    boost::numeric::ublas::iterator1_tag()); ++row_iterator)
        {
#endif
            neighbors.push_back(row_iterator.index2());
        }
    }


    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise

    static void InplaceMult(VectorType& rX, const double A)
    {

        if (A == 1.00){ //do nothing 
        }
        else{
            rX *= A;
        }
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;

    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
        if (A == 1.00)
            rX = rY;
        else if (A == -1.00)
        {
            InplaceMult(rX,-1);
        }
        else
        {
            rX = rY;
            rX *= A;
        }

    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;

    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        if (A == 1.00)
            rX += rY;
        else if (A == -1.00)
            rX -= rY;
        else
            rX.Add(A,rY);
    }

    //********************************************************************

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        Assign(rZ, A, rX); //rZ = A*rX
        UnaliasedAdd(rZ, B, rY); //rZ += B*rY
    }

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        InplaceMult(rY, B);
        UnaliasedAdd(rY, A, rX);
    }


    /// rA[i] * rX

    static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        return inner_prod(row(rA, i), rX);
    }


    static void SetValue(VectorType& rX, IndexType local_i, TDataType value)
    {
        rX[local_i] = value;
    }

    /// rX = A
    static void Set(VectorType& rX, TDataType A)
    {
        rX.SetValue(A);
    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m, n, false);
    }

    static void Resize(MatrixPointerType& pA, SizeType m, SizeType n)
    {
        pA->resize(m, n, false);
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n, false);
    }

    static void Resize(VectorPointerType& pX, SizeType n)
    {
        pX->resize(n, false);
    }

    static void Clear(MatrixPointerType& pA)
    {
        pA->clear();
        pA->resize(0, 0, false);
    }

    static void Clear(VectorPointerType& pX)
    {
        pX->clear();
        pX->resize(0, false);
    }

    template<class TOtherMatrixType>
    inline static void ResizeData(TOtherMatrixType& rA, SizeType m)
    {
        rA.resize(m, false);
        //            std::fill(rA.begin(), rA.end(), TDataType());
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        DataType* vals = rA.value_data().begin();
        #pragma omp parallel for firstprivate(m)
        for(int i=0; i<static_cast<int>(m); ++i)
            vals[i] = TDataType();
#endif
    }

    //TODO: remove
    inline static void ResizeData(CsrMatrix<TDataType>& rA, SizeType m)
    {
        rA.value_data().resize(m);
        std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
    }

    //TODO remove
    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        rX.resize(m, false);
        std::fill(rX.begin(), rX.end(), TDataType());
    }

    template<class TOtherMatrixType>
    inline static void SetToZero(TOtherMatrixType& rA)
    {
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        TDataType* vals = rA.value_data().begin();
        const int size = rA.value_data().end() - rA.value_data().begin();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            vals[i] = TDataType();
#endif
    }

    inline static void SetToZero(CsrMatrix<TDataType>& rA)
    {
        rA.SetValue(TDataType());
    }

    inline static void SetToZero(VectorType& rX)
    {
#ifndef _OPENMP
        std::fill(rX.begin(), rX.end(), TDataType());
#else
        const int size = rX.size();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            rX[i] = TDataType();
#endif
    }

    template<class TOtherMatrixType, class TEquationIdVectorType>
    inline static void AssembleLHS(
        MatrixType& A,
        TOtherMatrixType& LHS_Contribution,
        TEquationIdVectorType& EquationId
    )
    {
        A.Assemble(A,LHS_Contribution,EquationId);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    virtual std::string Info() const
    {
        return "KratosSpace";
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "KratosSpace";
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    //***********************************************************************

    inline static constexpr bool IsDistributed()
    {
        return false;
    }

    //***********************************************************************

    inline static TDataType GetValue(const VectorType& x, IndexType LocalI)
    {
        return x[LocalI];
    }
    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<IndexType>& IndexArray, TDataType* pValues)
    {
        KRATOS_TRY

        for (IndexType i = 0; i < IndexArray.size(); i++)
            pValues[i] = x[IndexArray[i]];

        KRATOS_CATCH("")
    }

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char* pFileName, /*const*/ TOtherMatrixType& rM, const bool Symmetric)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketMatrix(pFileName, rM, Symmetric);
    }

    template< class VectorType >
    static bool WriteMatrixMarketVector(const char* pFileName, const VectorType& rV)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketVector(pFileName, rV);
    }

    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosSpace & operator=(KratosSpace const& rOther);

    /// Copy constructor.
    KratosSpace(KratosSpace const& rOther);


    ///@}

}; // Class KratosSpace



///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    KratosSpace& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream,
// 				    const KratosSpace& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
///@}


} // namespace Kratos.

#endif // KRATOS_SPACE_H_INCLUDED  defined
