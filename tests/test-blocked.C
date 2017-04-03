/* tests/test-blocked.C
 * Copyright (C) 2017 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
 *
 * --------------------------------------------------------
 *
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
 *
 */


/*! @file  tests/test-blocked.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include "linbox/linbox-config.h"

#include "linbox/matrix/blockedmatrix/dense-block-allocator.h"
#include "linbox/matrix/blockedmatrix/blocked-matrix.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"

#include <iostream>
#include <string>
#include <sys/time.h>

using namespace LinBox;

template<class _F, class _M>
void blockedMul(const _F& F, _M& C, _M& A, _M& B){
	typedef typename _M::Block_t Block_t;
	BlasMatrixDomain<_F> BMD(F);

	for(size_t i = 0; i < C.blocksPerColumn(); i++){
		for(size_t j = 0; j < C.blocksPerRow(); j++){
			Block_t& C_ij = C.refBlock(i,j);
			for(size_t k = 0; k < A.blocksPerRow(); k++){
				const Block_t& A_ik = A.getBlock(i,k);
				const Block_t& B_kj = B.getBlock(k,j);

				BMD.axpyin(*C_ij,*A_ik,*B_kj);
			}
		}
	}
}

int main (int argc, char **argv){
	typedef Givaro::Modular<double> Field;
	
	bool pass = true;
	static size_t n = 128;
	static int q = 100003U;
	Field F(q);

	static Argument args[] = {
		{'n', "-n N", "Set dimension of the test matrices to NxN", TYPE_INT, &n},
		END_OF_ARGUMENTS
	};

	parseArguments(argc, argv, args);

	srand((unsigned)time(NULL));

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_UNIMPORTANT);

	commentator().start("Blocked matrix test suite", "BlockedMatrix");

	BlasMatrixDomain<Field> BMD(F);

	BlockedMatrix<Field, DenseBlockAllocator<Field> > A1(F,2,2,n/2,n/2);
	BlockedMatrix<Field, DenseBlockAllocator<Field> > B1(F,2,2,n/2,n/2);
	BlockedMatrix<Field, DenseBlockAllocator<Field> > C1(F,2,2,n/2,n/2);

	BlasMatrix<Field> A2(F,n,n);
	BlasMatrix<Field> B2(F,n,n);
	BlasMatrix<Field> C2(F,n,n);
	BlasMatrix<Field> C_temp(F,n,n);

	Field::RandIter G(F);

	for(size_t i = 0; i < n; i++){
		for(size_t j = 0; j < n; j++){
			A2.setEntry(i,j,G.random());
			B2.setEntry(i,j,G.random());
		}
	}

	A1.copyFromMatrix<BlasMatrix<Field> >(A2);
	B1.copyFromMatrix<BlasMatrix<Field> >(B2);

	BMD.mul(C2,A2,B2);
	blockedMul(F,C1,A1,B1);

	C1.copyToMatrix(C_temp);

	pass = BMD.areEqual(C2,C_temp);

	commentator().stop(MSG_STATUS(pass), " Blocked matrix test suite", (const char*)0);
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
