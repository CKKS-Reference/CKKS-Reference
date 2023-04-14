#include <stdlib.h>
#include <set>
#include <unordered_set>
#include <map>
#include "NumberTheory.h"

extern "C"{
	void     mul_asm( uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo );
	uint64_t div_asm( uint64_t hi, uint64_t lo, uint64_t a, uint64_t* r );
	uint64_t mod_asm( uint64_t hi, uint64_t lo, uint64_t a );
}

//
uint64_t mul_mod( uint64_t a, uint64_t b, uint64_t q){
	uint64_t hi,lo; mul_asm(a,b,&hi,&lo);
	return mod_asm(hi,lo,q);
}

//
uint64_t power_mod( uint64_t a, uint64_t n, uint64_t q ){
	uint64_t m=a%q, out=1;
	while(n>0){
		if(n%2==1) out=mul_mod(out,m,q);
		m=mul_mod(m,m,q);
		n>>=1;
	}
	return out;
}

//
uint64_t inv_mod( uint64_t a, uint64_t q ){ return power_mod(a,q-2,q);}
//
#if RAND_MAX/256 >= 0xFFFFFFFFFFFFFF
  #define LOOP_COUNT 1
#elif RAND_MAX/256 >= 0xFFFFFF
  #define LOOP_COUNT 2
#elif RAND_MAX/256 >= 0x3FFFF
  #define LOOP_COUNT 3
#elif RAND_MAX/256 >= 0x1FF
  #define LOOP_COUNT 4
#else
  #define LOOP_COUNT 5
#endif

uint64_t rand_uint64(){
	uint64_t r = 0;
	for (int i=LOOP_COUNT; i > 0; i--) {
		r = r*(RAND_MAX + (uint64_t)1) + rand();
	}
	return r;
}

//
bool isprime(uint64_t n){
	if(n<=1) return false;
	if(n<=3) return true;
	if(n%2==0 || n%3==0) return false;
	for(uint64_t i=5;i*i<=n;i+=6)
		if((n%i==0) || (n%(i+2)==0))
			return false;
	return true;
}
//
void findPrimeFactors( uint64_t n, std::unordered_set<uint64_t>& s ){
	while(n%2==0){ s.insert(2); n=n/2;}
	while(n%3==0){ s.insert(3); n=n/3;}
	for(uint64_t i=5;i*i<=n;i=i+6){
		while(n% i   ==0){s.insert(i  );n=n/ i   ;}
		while(n%(i+2)==0){s.insert(i+2);n=n/(i+2);}
	}
	if(n>1)
		s.insert(n);
}
//
uint64_t findPrimitive(uint64_t n){
	// MODIFIED1022
	static std::map<uint64_t, uint64_t> table;
	static bool table_initialized= false;
	if(table_initialized==false)
	{
		uint64_t q[60]={ 1125899915231233ULL,1125899927027713ULL,1125899930959873ULL,1125899935940609ULL,1125899941445633ULL,1125899963990017ULL,
						 1125899976572929ULL,1125899978145793ULL,1125899989155841ULL,1125899996233729ULL,1125899996495873ULL,1125899998068737ULL,
						 1125899998855169ULL,1125900015108097ULL,1125900015370241ULL,1125900027166721ULL,1125900043419649ULL,1125900047351809ULL,
						 1125900052856833ULL,1125900053905409ULL,1125900055478273ULL,1125900063866881ULL,1125900070944769ULL,1125900076711937ULL,
						 1125900078284801ULL,1125900079595521ULL,1125900083003393ULL,1125900093751297ULL,1125900094013441ULL,1125900097159169ULL,
			             1152921504616808449ULL,1152921504618381313ULL,1152921504622575617ULL,1152921504634109953ULL,1152921504643809281ULL,
						 1152921504650100737ULL,1152921504663994369ULL,1152921504666615809ULL,1152921504672645121ULL,1152921504687849473ULL,
						 1152921504690208769ULL,1152921504690995201ULL,1152921504694927361ULL,1152921504695451649ULL,1152921504695713793ULL,
						 1152921504697286657ULL,1152921504717733889ULL,1152921504722190337ULL,1152921504730316801ULL,1152921504730841089ULL,
						 1152921504738181121ULL,1152921504741851137ULL,1152921504746569729ULL,1152921504756006913ULL,1152921504760201217ULL,
						 1152921504760987649ULL,1152921504767016961ULL,1152921504770424833ULL,1152921504783794177ULL,1152921504786677761ULL};
		uint64_t psi[60]={5,5,5,3,3,5,7,13,37,11,3,3,3,7,3,3,19,13,10,3,3,11,23,3,3,7,3,5,6,3,11,5,3,5,3,3,7,3,11,3,3
			,3,3,13,3,3,3,11,3,11,3,10,7,10,3,3,43,3,3,17};
		for(int i=0;i<60;i++)
			table[q[i]]=psi[i];
		table_initialized = true;
	}

	//
	if(table.find(n) != table.end())
		return table[n];
	else{
		if(isprime(n)==false)
			return -1;
		std::unordered_set<uint64_t> s;
		findPrimeFactors(n-1,s);
		for(uint64_t r=2;r<=n-1;r++){
			bool flag=false;
			for(auto it=s.begin();it!=s.end();it++){
				if(power_mod(r,(n-1)/(*it),n)==1){
					flag=true;
					break;
				}
			}
			if(flag==false) {
				table[n] = r;
				return r;
			}
		}
		return 0;
	}
}
