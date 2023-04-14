#include <memory.h>
#include <stdio.h>
#include <assert.h>
#include "BigInt.h"

//------------------------------------------------------
// memory related functions
//------------------------------------------------------
BigInt::BigInt( int size){
	this->size=size; data=new uint64_t[size];
	memset(data,0,sizeof(uint64_t)*size);
}
BigInt::BigInt( const uint64_t* data, int size ){
	this->size=size; this->data=new uint64_t[size];
	memcpy(this->data,data,sizeof(uint64_t)*size);
}
BigInt::BigInt( const BigInt& Z){
	size=Z.size; data=new uint64_t[size];
	memcpy(data,Z.data,sizeof(uint64_t)*size);
}
BigInt::~BigInt(){delete[] data;}
void BigInt::resize( int new_size, bool copy ){
	if(size==new_size) return; 
	uint64_t* new_data = new uint64_t[new_size];
	memset(new_data,0,sizeof(uint64_t)*new_size);
	if(copy)
		memcpy(new_data,data,sizeof(uint64_t)*((new_size<size)?new_size:size));
	delete[] data; data=new_data; size=new_size;
}
void BigInt::operator=( const BigInt& Z){
	resize(Z.size,false);
	memcpy(data,Z.data,sizeof(uint64_t)*size);
}
//------------------------------------------------------
// accessors
//------------------------------------------------------
void BigInt::print( const char* name ) const{
	printf("%s:",name);
	for(int i=size-1;i>=1;i--)
		printf("%" PRIu64 "*beta^%d+",data[i],i);
	printf("%" PRIu64 "\n",data[0]);
}
//------------------------------------------------------
// comparison and arithmetics
//------------------------------------------------------
bool BigInt::operator>=( const BigInt& B) const{ 
	const BigInt& A = *this;
	int Asize = A.size, Bsize=B.size;
	while(Asize!=Bsize){
		if(Asize>Bsize) if(A.data[Asize-1]>0) return true; else Asize--;
		if(Asize<Bsize) if(B.data[Bsize-1]>0) return false; else Bsize--;
	}
	for(int i=Asize-1;i>=0;i--){
		if(A.data[i]>B.data[i]) return true;
		if(A.data[i]<B.data[i]) return false;
	}
	return true;
}

uint64_t BigInt::mod(uint64_t a ) const{
	uint64_t r=0;
	for(int i=size-1;i>=0;i--)
		r=mod_asm(r,data[i]%a,a);
	return r;
}
void BigInt::add( const BigInt& B, BigInt& C ) const{
	const BigInt& A=*this;
	int Csize = (A.size>B.size)?A.size:B.size;
	C.resize(Csize,false); 
	int i=0; bool c=0;
	for(;(i<A.size)&&(i<B.size);i++){
		C[i]=A[i]+B[i]; bool c_next=C[i]<A[i];
		C[i]=C[i]+c; 
		c=(c&&(C[i]==0))||c_next;
	}
	for(;i<Csize;i++){
		if(i<A.size) C[i]=A[i]+c;
		if(i<B.size) C[i]=B[i]+c;
		c=c&&(C[i]==0);
	}
	if(c){
		C.resize(Csize+1,true);
		C[Csize]=1;
	}
}
void BigInt::sub( const BigInt& B, BigInt& C ) const{
	const BigInt& A=*this;
	assert(A>=B); int Csize=A.size; C.resize(Csize,false);
	int i=0; bool b=0;
	for(;i<(B.size>Csize?Csize:B.size);i++){
		C[i]=A[i]-B[i];
		bool b_next=(C[i]>A[i])||(b&&(C[i]==0));
		C[i]=C[i]-b; b=b_next;
	}
	for(;i<Csize;i++){
		bool b_next = b&&(A[i]==0);
		C[i]=A[i]-b;
		b=b_next;
	}
	while( (Csize>0)&&(C[Csize-1]==0))
		Csize--;
	C.resize(Csize,true);
}

void BigInt::mul( uint64_t b, BigInt& C ) const{
	BigInt Cl(size  ); const BigInt& A=*this;
	BigInt Ch(size+1); Ch[0]=0;
	for(int i=0;i<size;i++)
		mul_asm(A[i],b,&Ch[i+1],&Cl[i]);
	Cl.add(Ch,C);
}

void BigInt::div( uint64_t b, BigInt& Q, uint64_t& r ) const{
	Q.resize(size); r=0;
	const BigInt& A = *this;
	for(int i=size-1;i>=0;i--)
		Q[i]=div_asm(r,A[i],b,&r);
}

double BigInt::to_real()const{
	double out=0; double beta=(1ULL<<63); beta*=2;
	for(int i=size-1;i>=0;i--)
		out=out*beta+data[i];
	return out;
}