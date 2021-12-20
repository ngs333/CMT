#ifndef MEMORY_AUX
#define MEMORY_AUX


/***
This file contains several classes used as alternative memory management aids for the search trees.
In particular, for the DistanceInterval snd PQNode objects.
***/

/*
  MemoryFSA:
  - A fixed size array based memory class; the client needs to know ahead of time of an upper bound on the memory to be used
  -   during the entire run of the program.
  - It overrides the new and delete constructors for classes that inherit this (See,for example, DistanceInterval).
  - Data is allocated in a static std::vector, so the functions reflct this. e.g. static data is not deleted.
*/

template <typename T>
class MemoryFSA {
protected:
	static std::vector<T> data;
	static unsigned int used;
public:
	static void initialize(unsigned int n) {
		data.reserve(n);
		used = 0;
	}
	//static void free() {used = 0;}
	static void reset() {
		used = 0;
	}
	static void* operator new(size_t sz) {
		if (used + 1 > data.capacity())
			error("DI::used += size >= maxSize ");
		T* ed = data.data();
		ed += used;
		used += 1;
		return ed;
	}
	// Overloading CLass specific new operator
	static void* operator new[](size_t sz) {
		unsigned int size = sz - sizeof(int*);
		size = size / sizeof(T);
		if (used + size > data.capacity())
			error("DI::used += size >= maxSize ");
		T* ed = data.data();
		ed += used;
		used += size;
		return ed;
	}
		// Overloading class specific delete operators:
		//Note that static memory should not be deleted
		static void operator delete(void* m) {}
		static void operator delete[](void* m) {/*delete[] m;*/}
};

/***
	A fixed size vector providing minimal size overhead;
	-Its may be used with MemoryFSA class, which provides fixed size static memory arrays.
	- Overhead is 16 bytes, vs 24 bytes for the std::vector in 64bit Windows 10).
	- Only the interfaces used by the CMT are included, and all of these except for function reset are
	  the same as the std::vector. In this way it will be easy to replace it with std::vector when  desired.
*/
template<typename T>
class FSVector {
private:
	unsigned int mSize;
	unsigned int ctr;
	T* data;
public:
	/*
	Constructor does not do much. You must call function reserve().
	*/
	FSVector() {
		mSize = 0;
		ctr = 0;
		data = nullptr;
	}
	~FSVector() {
		clear();
	}
	T& back() { return data[ctr - 1]; }
	T& operator [](unsigned int n) { return data[n]; }
	void push_back(const T& nd) {
		if (ctr >= mSize) {
			error("VectorFM:push_back ctr >= size");
		}
		data[ctr] = nd;
		++ctr;
	}
	unsigned int size() { return ctr; }
	void reserve(unsigned int sz) {
		data = new T[sz];
		mSize = sz;
		ctr = 0;
	}
	void reset() {ctr = 0;}
	
	//
	void clear() {
		delete[] data; //The delete[] operator of class T
		data = nullptr;
		mSize = 0;
		ctr = 0;
	}
	
};

#endif