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

#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/blockedmatrix/dense-block-allocator.h"
#include "linbox/matrix/blockedmatrix/blocked-matrix.h"
#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"

#include <iostream>
#include <string>
#include <sys/time.h>

using namespace LinBox;

const int maxpretty = 35;

std::string pretty(std::string a) {
	std::string blank;
	blank = a;
	int msgsize= maxpretty - (int)blank.size();
	std::string dot(".");
	for(int i=0;i<msgsize ;++i){
		blank += dot;
	}
	return blank;
}

template <class Field, class Matrix>
Matrix& fillRandom(const Field& F, Matrix& M){
	size_t r = M.rowdim();
	size_t c = M.rowdim();

	typename Field::RandIter G(F);

	for(size_t i = 0; i < r; i++){
		for(size_t j = 0; j < c; j++){
			M.setEntry(i,j,G.random());
		}
	}

	return M;
}

template <class Field>
bool testAdd(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing add"),"testAdd");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);

	MD.add(C1,A,B);
	BMD.add(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAdd");

	return ret;
}

template <class Field>
bool testAddin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing addin"),"testAddin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	B = fillRandom(F,B);
	blockB.copyFromMatrix(B);

	MD.addin(C1,B);
	BMD.addin(blockC,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAddin");

	return ret;
}

template <class Field>
bool testSub(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing sub"),"testSub");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);

	MD.sub(C1,A,B);
	BMD.sub(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testSub");

	return ret;
}

template <class Field>
bool testSubin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing subin"),"testSubin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	B = fillRandom(F,B);
	blockB.copyFromMatrix(B);

	MD.subin(C1,B);
	BMD.subin(blockC,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testSubin");

	return ret;
}

template <class Field>
bool testCopy(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing copy"),"testCopy");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	B = fillRandom(F,B);
	blockB.copyFromMatrix(B);

	MD.copy(C1,B);
	BMD.copy(blockC,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testCopy");

	return ret;
}

template <class Field>
bool testMul(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mul"),"testMul");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);

	MD.mul(C1,A,B);
	BMD.mul(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMul");

	return ret;
}

template <class Field>
bool testMulScale(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulscale"),"testMulScale");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);

	MD.mul(C1,2.0,A,B);
	BMD.mul(blockC,2.0,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulScale");

	return ret;
}

template <class Field>
bool testMulinLeft(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulin_left"),"testMulinLeft");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A1(F,n,n);
	Matrix B(F,n,n);
	Matrix A2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);

	A1 = fillRandom(F,A1);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A1);
	blockB.copyFromMatrix(B);

	MD.mulin_left(A1,B);
	BMD.mulin_left(blockA,blockB);
	
	blockA.copyToMatrix(A2);

	if(!MD.areEqual(A1,A2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulinLeft");

	return ret;
}

template <class Field>
bool testMulin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulin"),"testMulin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A1(F,n,n);
	Matrix B(F,n,n);
	Matrix A2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);

	A1 = fillRandom(F,A1);
	B = fillRandom(F,B);
	blockA.copyFromMatrix(A1);
	blockB.copyFromMatrix(B);

	MD.mulin(A1,B);
	BMD.mulin(blockA,blockB);
	
	blockA.copyToMatrix(A2);

	if(!MD.areEqual(A1,A2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulin");

	return ret;
}

template <class Field>
bool testMulinRight(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing mulin_right"),"testMulinRight");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B1(F,n,n);
	Matrix B2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B1 = fillRandom(F,B1);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B1);

	MD.mulin_right(A,B1);
	BMD.mulin_right(blockA,blockB);
	
	blockB.copyToMatrix(B2);

	if(!MD.areEqual(B1,B2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMulinRight");

	return ret;
}

template <class Field>
bool testAxpy(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axpy"),"testAxpy");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C(F,n,n);
	Matrix D1(F,n,n);
	Matrix D2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);
	BlockMatrix blockD(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C = fillRandom(F,C);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C);

	MD.axpy(D1,A,B,C);
	BMD.axpy(blockD,blockA,blockB,blockC);
	
	blockD.copyToMatrix(D2);

	if(!MD.areEqual(D1,D2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxpy");

	return ret;
}

template <class Field>
bool testAxpyin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axpyin"),"testAxpyin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C1 = fillRandom(F,C1);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C1);

	MD.axpyin(C1,A,B);
	BMD.axpyin(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxpyin");

	return ret;
}

template <class Field>
bool testMaxpy(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing maxpy"),"testMaxpy");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C(F,n,n);
	Matrix D1(F,n,n);
	Matrix D2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);
	BlockMatrix blockD(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C = fillRandom(F,C);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C);

	MD.maxpy(D1,A,B,C);
	BMD.maxpy(blockD,blockA,blockB,blockC);
	
	blockD.copyToMatrix(D2);

	if(!MD.areEqual(D1,D2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMaxpy");

	return ret;
}

template <class Field>
bool testMaxpyin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing maxpyin"),"testMaxpyin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C1 = fillRandom(F,C1);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C1);

	MD.maxpyin(C1,A,B);
	BMD.maxpyin(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMaxpyin");

	return ret;
}

template <class Field>
bool testAxmy(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axmy"),"testAxmy");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C(F,n,n);
	Matrix D1(F,n,n);
	Matrix D2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);
	BlockMatrix blockD(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C = fillRandom(F,C);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C);

	MD.axmy(D1,A,B,C);
	BMD.axmy(blockD,blockA,blockB,blockC);
	
	blockD.copyToMatrix(D2);

	if(!MD.areEqual(D1,D2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxmy");

	return ret;
}

template <class Field>
bool testAxmyin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing axmyin"),"testAxmyin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C1 = fillRandom(F,C1);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C1);

	MD.axmyin(C1,A,B);
	BMD.axmyin(blockC,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testAxmyin");

	return ret;
}

template <class Field>
bool testMuladd(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing muladd"),"testMuladd");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C(F,n,n);
	Matrix D1(F,n,n);
	Matrix D2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);
	BlockMatrix blockD(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C = fillRandom(F,C);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C);

	MD.muladd(D1,2.0,C,3.0,A,B);
	BMD.muladd(blockD,2.0,blockC,3.0,blockA,blockB);
	
	blockD.copyToMatrix(D2);

	if(!MD.areEqual(D1,D2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMuladd");

	return ret;
}

template <class Field>
bool testMuladdin(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing muladdin"),"testMuladdin");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A(F,n,n);
	Matrix B(F,n,n);
	Matrix C1(F,n,n);
	Matrix C2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);
	BlockMatrix blockB(F,2,2,n/2,n/2);
	BlockMatrix blockC(F,2,2,n/2,n/2);

	A = fillRandom(F,A);
	B = fillRandom(F,B);
	C1 = fillRandom(F,C1);
	blockA.copyFromMatrix(A);
	blockB.copyFromMatrix(B);
	blockC.copyFromMatrix(C1);

	MD.muladdin(2.0,C1,3.0,A,B);
	BMD.muladdin(2.0,blockC,3.0,blockA,blockB);
	
	blockC.copyToMatrix(C2);

	if(!MD.areEqual(C1,C2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testMuladdin");

	return ret;
}

template <class Field>
bool testLU(const Field& F, size_t n){
	typedef BlasMatrix<Field> Matrix;
	typedef BlockedMatrix<Field, DenseBlockAllocator<Field> > BlockMatrix;

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth(3);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_NORMAL);
	commentator().start(pretty("Testing LUFactorize"),"testLU");

	bool ret = true;
	BlasMatrixDomain<Field> MD(F);
	BlockedMatrixDomain<Field> BMD(F);

	Matrix A1(F,n,n);
	Matrix A2(F,n,n);

	BlockMatrix blockA(F,2,2,n/2,n/2);

	A1 = fillRandom(F,A1);
	blockA.copyFromMatrix(A1);

	LQUPMatrix<Field> LU(A1);
	BMD.LUFactorize(blockA);
	
	blockA.copyToMatrix(A2);

	if(!MD.areEqual(A1,A2)){ ret = false; }

	commentator().stop(MSG_STATUS(ret), (const char*)0, "testLU");

	return ret;
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

	if(n % 2){ n += 1;}

	srand((unsigned)time(NULL));

	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass(INTERNAL_DESCRIPTION).setMaxDetailLevel(Commentator::LEVEL_UNIMPORTANT);

	commentator().start("Blocked matrix test suite", "BlockedMatrix");

	pass &= testAdd(F,n);
	pass &= testAddin(F,n);

	pass &= testSub(F,n);
	pass &= testSubin(F,n);
	
	pass &= testCopy(F,n);

	pass &= testMul(F,n);
	pass &= testMulScale(F,n);
	pass &= testMulinLeft(F,n);
	pass &= testMulin(F,n);
	pass &= testMulinRight(F,n);

	pass &= testAxpy(F,n);
	pass &= testAxpyin(F,n);

	pass &= testMaxpy(F,n);
	pass &= testMaxpyin(F,n);

	pass &= testAxmy(F,n);
	pass &= testAxmyin(F,n);

	pass &= testMuladd(F,n);
	pass &= testMuladdin(F,n);

	pass &= testLU(F,n);

	commentator().stop(MSG_STATUS(pass), (const char*)0, "BlockedMatrix");
	return (pass ? 0 : -1);
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
