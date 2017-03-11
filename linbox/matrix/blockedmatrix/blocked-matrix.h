/*
 * Copyright (C) 2017 Matthew Wezowicz
 *
 * Written by :
 *               Matthew Wezowicz <mwezz@udel.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** @file matrix/blockedmatrix/blocked-matrix.h
 * @ingroup blockedmatrix
 * A \c BlockedMatrix<_Field,_Allocator>\c represents a matrix as an array of
 * submatrices.
 *
 */

#ifndef __LINBOX_matrix_blockedmatrix_blocked_matrix_H
#define __LINBOX_matrix_blockedmatrix_blocked_matrix_H

#include "linbox/matrix/blockedmatrix/block.h"

namespace LinBox
{
	/**
	 *
	 */
	template<class _Field, class _Allocator<_Field> >
	class BlockedMatrix{
	public:
		typedef _Field                                     Field;
		typedef typename Field::Element                  Element;    //!< Element type
		typedef _Allocator<Field>                      Allocator;
		typedef BlockedMatrix<Field,Allocator>            Self_t;    //!< Self type
		typedef const BlockedMatrix<Field,Allocator> constSelf_t;    //!< Self type
		typedef Block<Field,_Repr>                  Block<_Repr>;

                typedef Self_t                                matrixType;    //!< matrix type
                typedef constSelf_t                      constMatrixType;    //!< matrix type
                typedef Self_t                                  blasType;    //!< blas matrix type
	
	protected:
		size_t         _rows; //!< Number of rows of blocks
		size_t         _cols; //!< Number of columns of blocks
		size_t   _block_rows; //!< Number of rows per block
		size_t   _block_cols; //!< Number of columns per block
		Allocator _allocator;

	public:
		//////////////////
		// CONSTRUCTORS //
		//////////////////

		/** Allocates a new zero block matrix of \f$ M \times N\f$ 
		 * blocks, with each block \f$ m \times n\f$ in dimension.
		 * @param F
		 * @param m rows of blocks
		 * @param n cols of blocks
		 * @param block_m rows per block
		 * @param block_n cols per block
		 *
		 */
		BlockedMatrix(
			const Field& F,
			const size_t m,
			const size_t n,
			const size_t block_m,
			const size_t block_n);

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the overall number of rows in the matrix.
		 * @returns Overall number of rows in matrix
		 */
		size_t rowdim() const;

		/** Get the overall number of columns in the matrix.
		 * @returns Overall number of columns in matrix
		 */
		size_t coldim() const;

		/** Get the number of rows per block in the matrix.
		 * @returns Number of rows per block in matrix
		 */
		size_t blockRowdim() const;

		/** Get the number of columns per block in the matrix.
		 * @returns Number of columns per block in matrix
		 */
		size_t blockColdim() const;

		/** Get the number of blocks per row in the matrix.
		 * @returns Number of blocks per row in matrix
		 */
		size_t blocksPerRow() const;

		/** Get the number of blocks per column in the matrix.
		 * @returns Number of blocks per column in matrix
		 */
		size_t blocksPerColumn() const;


		//////////////////
		//    BLOCKS    //
		//////////////////

		/** Set the block at the (i, j) position to a_ij.
		 * @param i Block Row number, 0...blockRowdim() - 1
		 * @param j Block Column number 0...blockColdim() - 1
		 * @param a_ij Block to set
		 */
		template<class _Repr>
		void setBlock(size_t i, size_t j, const Block<_Repr>& a_ij);

		/** Get a writeable reference to the block in the (i, j) position.
		 * @param i Block Row index of entry
		 * @param j Block Column index of entry
		 * @returns Reference to block
		 */
		template<class _Repr>
		Block<_Repr>& refBlock(size_t i, size_t j);

		/** Get a read-only reference to the block in the (i, j) position.
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Const reference to block
		 */
		template<class _Repr>
		const Block<_Repr>& getBlock(size_t i, size_t j) const;

		/** Copy the (i, j) block into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Block in which to store result
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Reference to x
		 */
		template<class _Repr>
		Block<_Repr>& getBlock(Block<_Repr>& x, size_t i, size_t j) const;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/** Set the entry at the (i, j) position to a_ij.
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry(size_t i, size_t j, const Element& a_ij);

		/** Get a writeable reference to the entry in the (i, j) position.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		Element& refEntry(size_t i, size_t j);

		/** Get a read-only reference to the entry in the (i, j) position.
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		const Element& getEntry(size_t i, size_t j) const;

		/** Copy the (i, j) entry into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @returns Reference to x
		 */
		Element& getEntry(Element& x, size_t i, size_t j) const;

	}; // end of class BlockedMatrix

	/**
	 *
	 */
	template<class _BlockedMatrix>
	class AlignedBlockedSubmatrix{
	public:
		typedef typename _BlockedMatrix::Field                 Field;
		typedef typename Field::Element                      Element;
		typedef typename _BlockedMatrix::Allocator         Allocator;
		typedef AlignedBlockedSubmatrix<_Matrix>              Self_t;
		typedef const AlignedBlockedSubmatrix<_Matrix>   constSelf_t;
		typedef Block<Field,_Repr>                      Block<_Repr>;
		typedef Self_t                                 subMatrixType;
		typedef constSelf_t                       constSubMatrixType;
		typedef typename _BlockedMatrix::Self_t           matrixType;
		typedef typename _BlockedMatrix::constSelf_t constMatrixType;
		typedef matrixType                                  blasType;

	protected:
		size_t         _rows; //!< Number of rows of blocks
		size_t         _cols; //!< Number of columns of blocks
		size_t   _block_rows; //!< Number of rows per block
		size_t   _block_cols; //!< Number of columns per block
		size_t           _r0; //!< Upper left corner block row of AlignedBlockedSubmatrix in \p _Mat
		size_t           _c0; //!< Upper left corner block column of AlignedBlockedSubmatrix in \p _Mat
		_BlockedMatrix& _Mat; //!< Parent BlockedMatrix
	
	public:	
		//////////////////
		// CONSTRUCTORS //
		//////////////////

		/** Constructor from an existing @ref BlasMatrix and dimensions.
		 * @param M Pointer to @ref BlockedMatrix of which to construct submatrix
		 * @param rowbeg Starting block row
		 * @param colbeg Starting block column
		 * @param m rows of blocks 
		 * @param n cols of blocks
		 */
		AlignedBlockedSubmatrix(
			constMatrixType& M,
			const size_t rowbeg,
			const size_t colbeg,
			const size_t m,
			const size_t n);

		AlignedBlockedSubmatrix(
			matrixType& M,
			const size_t rowbeg,
			const size_t colbeg,
			const size_t m,
			const size_t n);

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the overall number of rows in the submatrix.
		 * @returns Overall number of rows in submatrix
		 */
		size_t rowdim() const;

		/** Get the overall number of columns in the submatrix.
		 * @returns Overall number of columns in submatrix
		 */
		size_t coldim() const;

		/** Get the number of rows per block in the submatrix.
		 * @returns Number of rows per block in submatrix
		 */
		size_t blockRowdim() const;

		/** Get the number of columns per block in the submatrix.
		 * @returns Number of columns per block in submatrix
		 */
		size_t blockColdim() const;

		/** Get the number of blocks per row in the submatrix.
		 * @returns Number of blocks per row in submatrix
		 */
		size_t blocksPerRow() const;

		/** Get the number of blocks per column in the submatrix.
		 * @returns Number of blocks per column in submatrix
		 */
		size_t blocksPerColumn() const;

		//////////////////
		//    BLOCKS    //
		//////////////////

		/** Set the block at the (i, j) position to a_ij.
		 * @param i Block Row number, 0...blockRowdim() - 1
		 * @param j Block Column number 0...blockColdim() - 1
		 * @param a_ij Block to set
		 */
		template<class _Repr>
		void setBlock(size_t i, size_t j, const Block<_Repr>& a_ij);

		/** Get a writeable reference to the block in the (i, j) position.
		 * @param i Block Row index of entry
		 * @param j Block Column index of entry
		 * @returns Reference to block
		 */
		template<class _Repr>
		Block<_Repr>& refBlock(size_t i, size_t j);

		/** Get a read-only reference to the block in the (i, j) position.
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Const reference to block
		 */
		template<class _Repr>
		const Block<_Repr>& getBlock(size_t i, size_t j) const;

		/** Copy the (i, j) block into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Block in which to store result
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Reference to x
		 */
		template<class _Repr>
		Block<_Repr>& getBlock(Block<_Repr>& x, size_t i, size_t j) const;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/** Set the entry at the (i, j) position to a_ij.
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry(size_t i, size_t j, const Element& a_ij);

		/** Get a writeable reference to the entry in the (i, j) position.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		Element& refEntry(size_t i, size_t j);

		/** Get a read-only reference to the entry in the (i, j) position.
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		const Element& getEntry(size_t i, size_t j) const;

		/** Copy the (i, j) entry into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @returns Reference to x
		 */
		Element& getEntry(Element& x, size_t i, size_t j) const;

	}; // end of class AlignedBlockedSubmatrix
} // end of namespace LinBox

#include "blocked-matrix.inl"
#include "blocked-submatrix.inl"

#endif // __LINBOX_matrix_blockedmatrix_blocked_matrix_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
