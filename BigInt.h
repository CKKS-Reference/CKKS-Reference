#pragma once
#include <stdint.h>
#include <inttypes.h>

extern "C"{
	void     mul_asm( uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo );
	uint64_t div_asm( uint64_t hi, uint64_t lo, uint64_t a, uint64_t* r );
	uint64_t mod_asm( uint64_t hi, uint64_t lo, uint64_t a );
}

class BigInt{
public:
	//------------------------------------------------------
	// memory related functions
	//------------------------------------------------------
	BigInt( int size=1);
	BigInt( const uint64_t* data, int size );
	BigInt( const BigInt& );
	~BigInt();
	void resize( int size, bool copy=false );
	void operator=( const BigInt& );
	//------------------------------------------------------
	// accessors
	//------------------------------------------------------
	uint64_t  operator[](int i)const{ return data[i]; }
	uint64_t& operator[](int i)     { return data[i]; }
	void print( const char* name ) const;
	//------------------------------------------------------
	// comparison and arithmetics
	//------------------------------------------------------
	bool operator>=( const BigInt& ) const;
	uint64_t mod(uint64_t a ) const;
	void add( const BigInt& B, BigInt& C ) const; // C=A+B
	void sub( const BigInt& B, BigInt& C ) const; // C=A-B
	void mul( uint64_t      b, BigInt& C ) const; // C=b*A
	void div( uint64_t      b, BigInt& Q, uint64_t& r ) const; // A=bQ+r
	double to_real() const;
protected:
	uint64_t* data;
	int       size;
};