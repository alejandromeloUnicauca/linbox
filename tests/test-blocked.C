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

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/blockedmatrix/dense-block-allocator.h"
#include "linbox/matrix/blockedmatrix/blocked-matrix.h"

#include <iostream>

using namespace LinBox;

int main (int argc, char **argv){
	bool pass = true;
	static size_t n = 128;
	static int q = 100003U;

	typedef Givaro::Modular<double> Field;
	Field F(q);

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_UNIMPORTANT);

	commentator().start("Blocked matrix test suite", "BlockedMatrix");
	
	BlockedMatrix<Field, DenseBlockAllocator<Field> > A1(F,n,n,n/2,n/2);
	BlockedMatrix<Field, DenseBlockAllocator<Field> > B1(F,n,n,n/2,n/2);
	BlockedMatrix<Field, DenseBlockAllocator<Field> > C1(F,n,n,n/2,n/2);

	BlasMatrix<Field> A2(F,n,n);
	BlasMatrix<Field> B2(F,n,n);
	BlasMatrix<Field> C2(F,n,n);

	Field::Element tmp;
	Field::RandIter G(F);

	for(size_t i = 0; i < n; i++){
		for(size_t j = 0; j < n; j++){
			tmp = G.random();
			A1.setEntry(i,j,tmp);
			A2.setEntry(i,j,tmp);

			tmp = G.random();
			B1.setEntry(i,j,tmp);
			B2.setEntry(i,j,tmp);
		}
	}

	commentator().stop("Blocked matrix test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
