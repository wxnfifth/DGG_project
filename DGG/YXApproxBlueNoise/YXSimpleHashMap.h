#ifndef __YX_SIMPLE_HASH_MAP_H__
#define __YX_SIMPLE_HASH_MAP_H__

#ifndef UINT32_DEFINED
#define UINT32_DEFINED
typedef unsigned int uint32;
#endif

// default hash function
template <class KEY_TYPE>
struct YXHashFunction{
	uint32 operator()(KEY_TYPE & key){
		uint32 hash = 0;
		uint32 * p = (uint32 *) (&key);
		int k = sizeof(KEY_TYPE) / 4;
		for(int i = 0; i < k; ++i){
			hash = hash * 9997 + p[i]; // signature. similar like Rabin-Karp
		}
		return hash;
	}
};

template <class KEY_TYPE, class VALUE_TYPE, class HASH_FUNCTION = YXHashFunction<KEY_TYPE> >
class YXSimpleHashMap{
	int bucketSize;
	std::vector < std::map <KEY_TYPE, VALUE_TYPE> > bucket;
	HASH_FUNCTION hash_function;
public:
	void init(int __bucketSize){
		bucketSize = __bucketSize;
		bucket.clear();
		bucket.resize(bucketSize);
	}
	void set(KEY_TYPE key, VALUE_TYPE value){
		uint32 hash = hash_function(key) % bucketSize;
		bucket[hash][key] = value;
	}
	VALUE_TYPE * get(KEY_TYPE key){
		uint32 hash = hash_function(key) % bucketSize;
		std::map <KEY_TYPE, VALUE_TYPE>::iterator itr = bucket[hash].find(key);
		if(itr == bucket[hash].end()) return NULL;
		return &(itr->second);
	}
	void clear(){
		bucket.clear();
	}
};

#endif // __YX_SIMPLE_HASH_MAP_H__