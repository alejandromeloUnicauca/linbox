/*
 * Copyright (C) 2017      Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file matrix/matrixdomain/blocked-domain.h
 * @ingroup matrixdomain
 * @ingroup blocked
 * @brief NO DOC
 */

#ifndef __LINBOX_matrix_matrixdomain_blockeddomain_H
#define __LINBOX_matrix_matrixdomain_blockeddomain_H

#include "linbox/matrix/matrixdomain/blas-matrix-domain.h"
#include "linbox/matrix/matrixdomain/blas-matrix-domain.h"

namespace LinBox{

	/**
	 *
	 */
	template <class Field_>
	class BlockedMatrixDomain{
	public:
		typedef Field_                          Field;
		typedef typename Field::Element         Element;

	protected:
		const Field& _F;
		BlasMatrixDomain<Field> BMD;
	
	public:

		//! Constructor of OpenCLDomain.
		BlockedMatrixDomain(const Field& F) : _F(F), BMD(&F){}

		//! Copy constructor
		BlockedMatrixDomain(const BlockedMatrixDomain<Field> & BlockMD) :
			// Init list begins here.
			_F(BlockMD._F),
			BMD(BlockMD.BMD){}

		//! Deconstructor
		~BlockedMatrixDomain(){}

		//! Field accessor
		const Field& field() const{
			return _F;
		}

		/*
		 * Basics operation available matrix respecting BlasMatrix interface
		 */

		//! addition.
		//! C = A+B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& add(Operand1& C, const Operand2& A, const Operand3& B) const{
			for(size_t i = 0; i < C.blockRowdim(); i++){
				for(size_t j = 0; j < C.blockColdim(); j++){
					BMD.add(*(C.refBlock(i,j)),
						*(A.getBlock(i,j)),
						*(B.getBlock(i,j)));
				}
			}
			return C;
		}

		//! addition (in place)
		//! C += B
		template <class Operand1, class Operand3>
		Operand1& addin(Operand1& C, const Operand3& B) const{
			for(size_t i = 0; i < C.blockRowdim(); i++){
				for(size_t j = 0; j < C.blockColdim(); j++){
					BMD.addin(*(C.refBlock(i,j)),
						  *(B.getBlock(i,j)));
				}
			}
			return C;
		}

		//! substraction
		//! C = A-B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& sub(Operand1& C, const Operand2& A, const Operand3& B) const{
			for(size_t i = 0; i < C.blockRowdim(); i++){
				for(size_t j = 0; j < C.blockColdim(); j++){
					BMD.sub(*(C.refBlock(i,j)),
						*(A.getBlock(i,j)),
						*(B.getBlock(i,j)));
				}
			}
			return C;
		}

		//! substraction (in place)
		//! C -= B
		template <class Operand1, class Operand3>
		Operand1& subin(Operand1& C, const Operand3& B) const{
			for(size_t i = 0; i < C.blockRowdim(); i++){
				for(size_t j = 0; j < C.blockColdim(); j++){
					BMD.subin(*(C.refBlock(i,j)),
						  *(B.getBlock(i,j)));
				}
			}
			return C;
		}

		//! copy.
		//! B = A
		template <class Operand1, class Operand2>
		Operand1& copy(Operand1& B, const Operand2& A) const{
			for(size_t i = 0; i < B.blockRowdim(); i++){
				for(size_t j = 0; j < B.blockColdim(); j++){
					BMD.copy(*(B.refBlock(i,j)),
						 *(A.getBlock(i,j)));
				}
			}
			return B;
		}

		//! multiplication.
		//! C = A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladd(C,_F.zero,C,_F.one,A,B);
		}

		//! multiplication with scaling.
		//! C = alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& mul(
			Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			return muladd(C,_F.zero,C,alpha,A,B);
		}

		//! In place multiplication.
		//! A = A*B
		template <class Operand1, class Operand2>
		Operand1& mulin_left(Operand1& A, const Operand2& B) const{
			Operand1* tmp = new Operand1(A);
			mul(A,*tmp,B);
			delete tmp;
			return A;
		}

		template <class Operand1, class Operand2>
		Operand1& mulin(Operand1& A, const Operand2& B) const{
			return mulin_left<Operand1, Operand2>(A,B);
		}
		//! In place multiplication.
		//! B = A*B
		template <class Operand1, class Operand2>
		Operand2& mulin_right(const Operand1& A, Operand2& B) const{
			Operand2* tmp = new Operand2(B);
			mul(B,A,*tmp);
			delete tmp;
			return B;
		}


		//! axpy.
		//! D = A*B + C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.one,C,_F.one,A,B);
		}

		//! axpyin.
		//! C += A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axpyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.one,C,_F.one,A,B);
		}

		//! maxpy.
		//! D = C - A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.one,C,_F.mOne,A,B);
		}

		//! maxpyin.
		//! C -= A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& maxpyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.one,C,_F.mOne,A,B);
		}

		//! axmy.
		//! D= A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmy(
			Operand1& D,
			const Operand2& A,
			const Operand3& B,
			const Operand1& C) const{

			return muladd(D,_F.mOne,C,_F.one,A,B);
		}

		//! axmyin.
		//! C = A*B - C
		template <class Operand1, class Operand2, class Operand3>
		Operand1& axmyin(Operand1& C, const Operand2& A, const Operand3& B) const{
			return muladdin(_F.mOne,C,_F.one,A,B);
		}

		//!  general matrix-matrix multiplication and addition with scaling.
		//! D = beta.C + alpha.A*B
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladd(
			Operand1& D,
			const Element& beta,
			const Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			for(size_t i = 0; i < D.blockRowdim(); i++){
				for(size_t j = 0; j < D.blockColdim(); j++){
					BMD.muladd(*(D.refBlock(i,j)),
						   beta,
						   *(C.getBlock(i,j)),
						   alpha,
						   *(A.getBlock(i,0)),
						   *(B.getBlock(0,j)));

					for(size_t k = 1; k < A.blockColdim(); k++){
						BMD.muladdin(_F.one,
							     *(D.refBlock(i,j)),
							     alpha,
							     *(A.getBlock(i,k)),
							     *(B.getBlock(k,j)));
					}
				}
			}
			return C;

		}

		//! muladdin.
		//! C = beta.C + alpha.A*B.
		template <class Operand1, class Operand2, class Operand3>
		Operand1& muladdin(
			const Element& beta,
			Operand1& C,
			const Element& alpha,
			const Operand2& A,
			const Operand3& B) const{

			Operand1* tmp = new Operand1(C);
			muladd(C,beta,*tmp,alpha,A,B);
			delete tmp;
			return C;
		}

		/*
		 * Advanced operations.
		 */

		/**
		 * In-place LU factorization.
		 */
		template <class Operand1>
		Operand1& LUFactorize(Operand1& A){
			size_t n = A.getBlockColdim() - 1;
			
			for(size_t i = 0; i < n - 1; i++){
				AlignedBlockedSubmatrix<Operand1> M(A,i,i,n,n);
				LUStep(M);
			}

			LQUPMatrix<Field> LU(*(A.refBlock(n,n)));

			return A;
		}

	private:
		template <class Operand1>
		Operand1& LUStep(Operand1& M){
			size_t n = M.getBlockColdim() - 1;

			LQUPMatrix<Field> A(*(M.refBlock(0,0)));

			AlignedBlockedSubmatrix<Operand1> B(M,0,1,1,n);
			AlignedBlockedSubmatrix<Operand1> C(M,1,0,n,1);
			AlignedBlockedSubmatrix<Operand1> D(M,1,1,n,n);

			for(size_t i = 0; i < n; i++){
				A.left_Lsolve(B.refBlock(0,i));
				A.right_Usolve(C.refBlock(i,0));
			}

			maxpyin(D,C,B);

			return M;
		}
	}; // end of class BlockedMatrixDomain
} // end of namespace LinBox

#endif // __LINBOX_matrix_matrixdomain_blockeddomain_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
