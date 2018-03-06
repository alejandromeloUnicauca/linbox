/* Copyright (C) 2007 LinBox
 * Written by bds and zw
 *
 * author: B. David Saunders and Zhendong Wan
 * parallelized for BOINC computing by Bryan Youse
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_cra_mpi_H
#define __LINBOX_cra_mpi_H

#define MPICH_IGNORE_CXX_SEEK //BB: ???
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-cra2.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/util/mpicpp.h"
/*
template <typename T > class chooseMPItype;
template <> struct chooseMPItype<unsigned int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED;};
template <> struct chooseMPItype<unsigned long long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG_LONG;};
template <> struct chooseMPItype<unsigned long int>{ static constexpr MPI_Datatype val = MPI_UNSIGNED_LONG;};
#include <gmp++/gmp++.h>
#include <string>
*/
#include <unordered_set>

namespace LinBox
{

	template<class CRABase>
	struct MPIChineseRemainder  {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
	protected:
		CRABase Builder_;
		Communicator* _commPtr;
		unsigned int _numprocs;

	public:
		template<class Param>
		MPIChineseRemainder(const Param& b, Communicator *c) :
			Builder_(b), _commPtr(c), _numprocs(c->size())
		{}

		/** \brief The CRA loop.
		 *
		 * termination condition.
		 *
		 * \param Iteration  Function object of two arguments, \c
		 * Iteration(r, p), given prime \c p it outputs residue(s) \c
		 * r.  This loop may be parallelized.  \p Iteration must be
		 * reentrant, thread safe.  For example, \p Iteration may be
		 * returning the coefficients of the minimal polynomial of a
		 * matrix \c mod \p p.
		 @warning  we won't detect bad primes.
		 *
		 * \param primeg  RandIter object for generating primes.
		 * \param[out] res an integer
		 */
		template<class Function, class PrimeIterator>
		Integer & operator() (Integer& res, Function& Iteration, PrimeIterator& primeg)
		{
			//  defer to standard CRA loop if no parallel usage is desired
			if(_commPtr == 0 || _commPtr->size() == 1) {
				ChineseRemainder< CRABase > sequential(Builder_);
				return sequential(res, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();

			//  parent process
			if(process == 0 ){
				//  create an array to store primes
				int primes[procs - 1];
				DomainElement r;
				//  send each child process a new prime to work on
				for(int i=1; i<procs; i++){
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					primes[i - 1] = *primeg;
					_commPtr->send(primes[i - 1], i);
				}
				bool first_time = true;
				int poison_pills_left = procs - 1;
				//  loop until all execution is complete
				while( poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive sub-answers from child procs
					_commPtr->recv(r, MPI_ANY_SOURCE);
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
					//  assimilate results
					if(first_time){
						Builder_.initialize(D, r);
						first_time = false;
					}
					else
						Builder_.progress( D, r );
					//  queue a new prime if applicable
					if(! Builder_.terminated()){
						++primeg;
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
					//  send the prime or poison pill
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // end while
				return Builder_.result(res);
			}  // end if(parent process)
			//  child processes
			else{
				int pp;
				while(true){
					//  receive the prime to work on, stop
					//  if signaled a zero
					_commPtr->recv(pp, 0);
					if(pp == 0)
						break;
					Domain D(pp);
					DomainElement r; D.init(r);
					Iteration(r, D);
					//Comm->buffer_attach(rr);
					// send the results
					_commPtr->send(r, 0);
				}
				return res;
			}
		}

#if 0
		template<class V, class F, class P>
		V & operator() (V& res, F& it, P&primeg){ return res; }
#endif
		template<class Vect, class Function, class PrimeIterator>
		Vect & operator() (Vect& res, Function& Iteration, PrimeIterator& primeg)
		{
			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				ChineseRemainder< CRABase > sequential(Builder_);
				return sequential(res, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();
// 			std::vector<DomainElement> r;
			typename Rebind<Vect, Domain>::other r;

			//  parent propcess
			if(process == 0){
				int primes[procs - 1];
				Domain D(*primeg);
				//  for each slave process...
				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
					//  send the prime to a slave process
					_commPtr->send(primes[i - 1], i);
				}
				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
				while(poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive the beginnin and end of a vector in heapspace
					_commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);
					//  determine which process sent answer
					//  and give them a new prime
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
					Builder_.progress(D, r);
					//  if still working, queue a prime
					if(! Builder_.terminated()){
						++primeg;
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
					//  send the prime or poison
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // while
				return Builder_.result(res);
			}
			//  child process
			else{
				int pp;
				//  get a prime, compute, send back start and end
				//  of heap addresses
				while(true){
					_commPtr->recv(pp, 0);
					if(pp == 0)
						break;
					Domain D(pp);
					Iteration(r, D);
					_commPtr->send(r.begin(), r.end(), 0, 0);
				}
				return res;
			}
		}
	};
        
        
        
        
	template<class RatCRABase>
	struct MPIratChineseRemainder  {
		typedef typename RatCRABase::Domain	Domain;
		typedef typename RatCRABase::DomainElement	DomainElement;
	protected:
		RatCRABase Builder_;
		Communicator* _commPtr;
		unsigned int _numprocs;
                
	public:
		template<class Param>
		MPIratChineseRemainder(const Param& b, Communicator *c) :
			Builder_(b), _commPtr(c), _numprocs(c->size())
		{}

		template<class Function, class PrimeIterator>
		Integer & operator() (Integer& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
		{
			//  defer to standard CRA loop if no parallel usage is desired
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();

			//  parent process
			if(process == 0 ){
				//  create an array to store primes
				int primes[procs - 1];
				DomainElement r;
				//  send each child process a new prime to work on
				for(int i=1; i<procs; i++){
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					primes[i - 1] = *primeg;
					_commPtr->send(primes[i - 1], i);
				}
				bool first_time = true;
				int poison_pills_left = procs - 1;
				//  loop until all execution is complete
				while( poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive sub-answers from child procs
					_commPtr->recv(r, MPI_ANY_SOURCE);
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
					//  assimilate results
					if(first_time){
						Builder_.initialize( D, Iteration(r, D) );
						first_time = false;
					}
					else
						Builder_.progress( D, Iteration(r, D) );
					//  queue a new prime if applicable
					if(! Builder_.terminated()){
						++primeg;
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}
					//  send the prime or poison pill
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // end while
				return Builder_.result(num,den);
			}  // end if(parent process)
			//  child processes
			else{
				int pp;
				while(true){
					//  receive the prime to work on, stop
					//  if signaled a zero
					_commPtr->recv(pp, 0);
					if(pp == 0)
						break;
					Domain D(pp);
					DomainElement r; D.init(r);
					Iteration(r, D);
					//Comm->buffer_attach(rr);
					// send the results
					_commPtr->send(r, 0);
				}
				return num;
			}
		}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<BEGIN<<<<ROI<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		template<class Function, class PrimeIterator>
		BlasVector<Givaro::ZRing<Integer> > & operator() ( BlasVector<Givaro::ZRing<Integer> > & num, Integer& den, Function& Iteration, PrimeIterator& primeg)
		{

#if 1

			//  if there is no communicator or if there is only one process,
			//  then proceed normally (without parallel)
			if(_commPtr == 0 || _commPtr->size() == 1) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();
			//typename Rebind<BlasVector< Givaro::ZRing<Integer> > , Domain>::other r;
Domain D(*primeg);
BlasVector<Domain> r(D);
Timer chrono;

			//  parent propcess
			if(process == 0){
  std::unordered_set<int> prime_sent;
				int primes[procs - 1];
				//Domain D(*primeg);
				//  for each slave process...
				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) || prime_sent.find(*primeg) != prime_sent.end()) ++primeg;	//while(Builder_.noncoprime(*primeg)) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
prime_sent.insert(*primeg);
					//  send the prime to a slave process
					_commPtr->send(primes[i - 1], i);
				}  
				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
				while(poison_pills_left > 0 ){
					int idle_process = 0;
					//  receive the beginnin and end of a vector in heapspace
					_commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0); 

					//  determine which process sent answer
					//  and give them a new prime
					idle_process = (_commPtr->get_stat()).MPI_SOURCE;
					Domain D(primes[idle_process - 1]);
//chrono.start();//////////////////////////////////////////////////////////////////
					Builder_.progress(D, r);
//chrono.stop();///////////////////////////////////////////////////////////////////
//std::cout<<"Builder_.progress(D, r) in the manager process used CPU time (seconds): " <<chrono.usertime()<<std::endl;
                                        
					//  if still working, queue a prime
					if(! Builder_.terminated()){
						++primeg; while(Builder_.noncoprime(*primeg) || prime_sent.find(*primeg) != prime_sent.end()) ++primeg;
prime_sent.insert(*primeg);
						primes[idle_process - 1] = *primeg;
					}
					//  otherwise, queue a poison pill
					else{
						primes[idle_process - 1] = 0;
						poison_pills_left--;
					}


					//  send the prime or poison
					_commPtr->send(primes[idle_process - 1], idle_process);
				}  // while

				return Builder_.result(num,den);
                                
			}
			//  child process
			else{
				int pp;
				//  get a prime, compute, send back start and end
				//  of heap addresses
				while(true){
					_commPtr->recv(pp, 0);

					if(pp == 0)
						break;
					Domain D(pp);
//chrono.start();///////////////////////////////////////////////////////////////
					Iteration(r, D);
//chrono.stop();////////////////////////////////////////////////////////////////
//std::cout << "Iteration(r,D) in the worker process used CPU time (seconds):  " << chrono.usertime() << std::endl;
					_commPtr->send(r.begin(), r.end(), 0, 0); 
				}
			}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#else

			if(_commPtr == 0) {
				RationalRemainder< RatCRABase > sequential(Builder_);
				return sequential(num, den, Iteration, primeg);
			}

			int procs = _commPtr->size();
			int process = _commPtr->rank();
			//typename Rebind<BlasVector< Givaro::ZRing<Integer> > , Domain>::other r;
Domain D(*primeg);
Domain D2(*primeg);

BlasVector<Domain> r(D);
BlasVector<Domain> r2(D);


			//  Manager propcess
			if(process == 0){

int primes[procs - 1];
int primes2[procs - 1];

				//  send the prime to a slave process using Linear-array scatter
				for (int i=1; i<procs; i++) {
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
					_commPtr->send(primes[i - 1], process+1); 
				}
				for (int i=1; i<procs; i++) {
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes2[i - 1] = *primeg;
					_commPtr->send(primes2[i - 1], process+1); 
				}

				Builder_.initialize( D, Iteration(r, D) );
				int poison_pills_left = procs - 1;
std::vector<BlasVector<Domain>> R;
std::vector<BlasVector<Domain>> R2;

while(poison_pills_left > 0 ){

	for  (int i=1; i<procs; i++){
		_commPtr->recv(r.begin(), r.end(), i, 0); 
		R.push_back(r);
	} 

	for  (int i=1; i<procs; i++){
		_commPtr->recv(r.begin(), r.end(), i, 0); 
		R2.push_back(r);

	} 
 
	for  (int i=1; i<procs; i++){
		Domain D(primes[i - 1]);		
		Builder_.progress(D,R[i - 1]);

	} R.clear();

	for  (int i=1; i<procs; i++){
		Domain D(primes2[i - 1]);
			
		Builder_.progress(D,R2[i - 1]);
	} R2.clear();



					//  if still working, queue a prime
			if(! Builder_.terminated()){

				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes[i - 1] = *primeg;
					_commPtr->send(primes[i - 1], process+1); 
				}

				for(int i=1; i<procs; i++){
					//  generate a new prime
					++primeg; while(Builder_.noncoprime(*primeg) ) ++primeg;
					//  fix the array of currently sent primes
					primes2[i - 1] = *primeg;
					_commPtr->send(primes2[i - 1], process+1); 
				}

			}
					//  otherwise, queue a poison pill
			else{
				for(int i=1; i<procs; i++){
					 poison_pills_left--;
					//  fix the array of currently sent primes
					primes[i - 1] = 0;
					_commPtr->send(primes[i - 1], process+1); 
				}

				for(int i=1; i<procs; i++){
					 poison_pills_left--;
					//  fix the array of currently sent primes
					primes2[i - 1] = 0;
					_commPtr->send(primes2[i - 1], process+1); 
				}
			}

} 

return Builder_.result(num,den);

			}else{

				int p[procs-1], pp[procs-1];
				//  get a prime, compute, send back start and end
				//  of heap addresses

				while(true){

if(process<procs-1){
 _commPtr->recv(p[0], process-1);

for (int i=1; i<procs-1; i++) MPI_Sendrecv(&p[i-1], 1, MPI_INT, process+1, 0, &p[i], 1, MPI_INT, process-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
_commPtr->send(p[procs-2], process+1);

}else{
for (int j=0; j<procs-1; j++) _commPtr->recv(p[j], process-1);

}

if(process<procs-1){
 _commPtr->recv(pp[0], process-1);

for (int i=1; i<procs-1; i++) MPI_Sendrecv(&pp[i-1], 1, MPI_INT, process+1, 0, &pp[i], 1, MPI_INT, process-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

_commPtr->send(pp[procs-2], process+1);

}else{

for (int j=0; j<procs-1; j++) _commPtr->recv(pp[j], process-1);

}


 					if(p[process-1] == 0)
						break;
 					if(pp[process-1] == 0)
						break;

					Domain D(p[process-1]);
					Iteration(r, D);
					Domain D2(pp[process-1]);
					Iteration(r2, D2);

					_commPtr->send(r.begin(), r.end(), 0, 0); 
					_commPtr->send(r2.begin(), r2.end(), 0, 0); 

				}


			}//if(process == 0)//


#endif

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
		}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>END>>>>ROI>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



	};



}

#undef MPICH_IGNORE_CXX_SEEK
#endif // __LINBOX_cra_mpi_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
