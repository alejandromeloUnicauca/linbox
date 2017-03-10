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

/*! @file matrix/blockedmatrix/blocked-matrix.h
 * @ingroup blockedmatrix
 * A \c BlockedMatrix<\c _Field > represents a matrix as an array of submatrices.
 *
 */

#ifndef __LINBOX_matrix_blockedmatrix_blocked_matrix_H
#define __LINBOX_matrix_blockedmatrix_blocked_matrix_H

namespace LinBox // Forward declarations.
{
	template<class _Field>
	class Block;
} // end of namespace LinBox

namespace LinBox
{
	template<class _Field, class _Allocator<_Field> >
	class BlockedMatrix{
	public:
		typedef _Field                                     Field;
		typedef typename Field::Element                  Element;    //!< Element type
		typedef _Allocator<Field>                      Allocator;
		typedef BlockedMatrix<Field,Allocator>            Self_t;    //!< Self typeype
		typedef const BlockedMatrix<Field,Allocator> constSelf_t;    //!< Self typeype
		typedef Block<Field>                               Block;

                typedef Self_t                                matrixType;    //!< matrix type
                typedef constSelf_t                      constMatrixType;    //!< matrix type
                typedef Self_t                                  blasType;    //!< blas matrix type
	
	private:
		size_t        _row;
		size_t        _col;
		size_t _block_rows;
		size_t _block_cols;
		Allocator     _rep;

	public:
		BlockedMatrix(const Field &F,
			      const size_t &m,
			      const size_t &n,
			      const size_t &block_m,
			      const size_t &block_n);

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the number of rows in the matrix.
		 * @returns Number of rows in matrix
		 */
		size_t rowdim() const;

		/** Get the number of columns in the matrix.
		 * @returns Number of columns in matrix
		 */
		size_t coldim() const;

		/** Get the number of block rows in the matrix.
		 * @returns Number of block rows in matrix
		 */
		size_t blockRowdim() const;

		/** Get the number of block columns in the matrix.
		 * @returns Number of block columns in matrix
		 */
		size_t blockColdim() const;


		//////////////////
		//    BLOCKS    //
		//////////////////

		/** Set the block at the (i, j) position to a_ij.
		 * @param i Block Row number, 0...blockRowdim() - 1
		 * @param j Block Column number 0...blockColdim() - 1
		 * @param a_ij Block to set
		 */
		void setBlock(size_t i, size_t j, const Block &a_ij);

		/** Get a writeable reference to the block in the (i, j) position.
		 * @param i Block Row index of entry
		 * @param j Block Column index of entry
		 * @returns Reference to block
		 */
		Block& refBlock(size_t i, size_t j);

		/** Get a read-only reference to the block in the (i, j) position.
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Const reference to block
		 */
		const Block& getBlock(size_t i, size_t j) const;

		/** Copy the (i, j) block into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Block in which to store result
		 * @param i Block Row index
		 * @param j Block Column index
		 * @returns Reference to x
		 */
		Block &getBlock(Block &x, size_t i, size_t j) const;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/** Set the entry at the (i, j) position to a_ij.
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry(size_t i, size_t j, const Element &a_ij);

		/** Get a writeable reference to the entry in the (i, j) position.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		Element &refEntry(size_t i, size_t j);

		/** Get a read-only reference to the entry in the (i, j) position.
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		const Element &getEntry(size_t i, size_t j) const;

		/** Copy the (i, j) entry into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @returns Reference to x
		 */
		Element &getEntry(Element &x, size_t i, size_t j) const;


	}; // end of class BlockedMatrix
} // end of namespace LinBox

#endif // __LINBOX_densematrix_blas_matrix_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
