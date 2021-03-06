/* tests/test-block-wiedemann.C
 * evolved by -bds from test-solve.C
 *
 * --------------------------------------------------------
 *
 * Copyright (c) LinBox
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

/*! @file   tests/test-block-wiedemann.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */

// to print out the number of try for BW with sigma basis code
/*  not needed -- this will print when commentator is on (test report file set).
#define _BW_LASVEGAS_COUNT
*/

#include "linbox/linbox-config.h"
#include <iostream>

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/block-wiedemann.h"
#include "linbox/algorithms/coppersmith.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/scalar-matrix.h"

#include "test-common.h"

using namespace LinBox;

/* Tests nonsingular solving for a random right-hand side.
 *
 * Solver - a block Wiedemann solver (see coppersmith.h or block-wiedemann.h).
 * Blackbox - nonsingular matrix (rhs will be a random vector).
 *
 * Checks the result, returning true on success and false on failure
 */
template <class Solver, class Blackbox>
bool testBlockSolver(Solver & S, Blackbox & M, string desc){
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	typedef typename Blackbox::Field Field;
	typedef BlasVector<Field> Vector;
	bool pass = true;
	size_t n = M.coldim();
	Vector b(M.field(),n), x(M.field(),n), y(M.field(),n);
        typename Field::RandIter gen(M.field());
	RandomDenseStream<Field> s (M.field(), gen, n, 1);
	s.next (b);
	VectorDomain<Field> VD (M.field());
	VD.write (report << "Right-hand side: b =  ", b) << endl;

	S.solveNonSingular(x, M, b);

	VD.write (report << desc << " solution:  ", x) << endl;
	M.apply (y, x);
	VD.write ( report << "Output:           ", y) << endl;

	if (!VD.areEqual (y, b)) {
		pass = false;
		report << "ERROR: " << desc << " solution is incorrect" << endl;
	}
	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 9; // blocking + 1 <= n/2 is required.
//	static size_t N = 16;
	static size_t q = 65521U;
	static size_t blocking = 1; // if blocking is 0, default blocksize 8 is used. 
        static size_t seed = time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to n.", TYPE_INT,     &n },
//		{ 'N', "-N N", "Set blocking factor to N.", TYPE_INT,     &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		{ 'b', "-b N", "Set the blocking size", TYPE_INT, &blocking },
                { 's', "-s N", "Set the seed for randomness", TYPE_INT, &seed },
		END_OF_ARGUMENTS
	};
        
        commentator().setMaxDepth(-1);
        commentator().setMaxDetailLevel(-1); 
        
	parseArguments (argc, argv, args);

	typedef Givaro::Modular<double> Field;
	//typedef Givaro::Modular<uint32_t> Field;
	typedef BlasVector<Field> Vector;  

	Field F ( (uint32_t) q);
        Field::RandIter G(F, 0, seed); //random generator over F
        Field::NonZeroRandIter NzG(G); //non-zero random generator over F 

	MatrixDomain<Field> MD(F);

	commentator().start("block wiedemann test suite", "block-wiedemann");
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

        RandomDenseStream<Field, Vector, Field::NonZeroRandIter> s (F, NzG, n, 2);
	Vector d(F,n);
	s.next (d);
	for (size_t i = 0; i < d.size(); ++i) 
		if (F.isZero(d[i])) report << "Warning: zero in vector" << endl;

	// some matrices
	ScalarMatrix<Field> SC (F, n, n, F.one); // identity

	Diagonal <Field> D (d); // random nonsingular diagonal

	SparseMatrix<Field, SparseMatrixFormat::TPL> S (F, n, n);
	for (size_t i = 1; i < n; ++i) S.setEntry(i, i-1, F.one); // subdiag 
	for (size_t i = 0; i < n; ++i) S.setEntry(i, n-1, d[i]); // last col
	S.finalize(); // companion matrix of d

#if 0
// RCS is Yuhasz' Matrix Berlekamp Massey method.
	CoppersmithSolver< MatrixDomain<Field> > RCS(MD,blocking);


	commentator().start("Ident, CoppersmithSolver", "I-Coppersmith");
	pass = pass and testBlockSolver(RCS, SC, "ScalarMatrix(I), Matrix Berlekamp Massey");
	commentator().stop("Ident, CoppersmithSolver");

	commentator().start("Diag, CoppersmithSolver", "D-Coppersmith");
	pass = pass and testBlockSolver(RCS, D, "Diagonal, Matrix Berlekamp Massey");
	commentator().stop("Diag, CoppersmithSolver");

	commentator().start("Companion, CoppersmithSolver", "C-Coppersmith");
	pass = pass and testBlockSolver(RCS, S, "Companion, Matrix Berlekamp Massey");
	commentator().stop("Companion, CoppersmithSolver");
#endif
        
#if 1
// LBWS is Giorgi's block method, SigmaBasis based.

#ifdef __LINBOX_HAVE_OCL
// using this as a test of the opencl matrix domain
	typedef OpenCLMatrixDomain<Field> Context;
// but note, shouldn't need the ifdef.  OpenCLMatrixDomain defaults to BlasMatrixDomain calls.
#else
	typedef BlasMatrixDomain<Field> Context;
#endif
	Context BMD(F);
	BlockWiedemannSolver<Context> LBWS(BMD,blocking,blocking+1);

/*
	commentator().start("Ident, BlockWiedemannSolver", "I-Sigma Basis");
	pass = pass and testBlockSolver(LBWS, SC, "ScalarMatrix(I), Sigma Basis");
	commentator().stop("Ident, BlockWiedemannSolver");
*/

	commentator().start("Diag, BlockWiedemannSolver", "D-Sigma Basis");
	pass = pass and testBlockSolver(LBWS, D, "Diagonal, Sigma Basis");
        commentator().stop(MSG_STATUS (pass), (const char *) 0,"Diagonal, Sigma Basis");

	commentator().start("Companion, BlockWiedemannSolver", "C-Sigma Basis");
	pass = pass and testBlockSolver(LBWS, S, "Companion, Sigma Basis");
        commentator().stop(MSG_STATUS (pass), (const char *) 0,"Companion, Sigma Basis");
#endif
        
        commentator().stop(MSG_STATUS (pass), (const char *) 0,"block wiedemann test suite");
        //std::cout << (pass ? "passed" : "FAILED" ) << std::endl;

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
