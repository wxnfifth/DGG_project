#ifndef LINEAR_QUEUE_H
#define LINEAR_QUEUE_H

//#define OUTPUTINFO

#include <vector>
#include "MyList.h"
#include "ListGroup.h"
#include "geodesic_memory.h"
#include <cstdio>

template<typename T>
class Linear_Queue
{
public:

	struct BinItem
	{
		T element;
		unsigned index;
	};

	Linear_Queue();
	Linear_Queue(int w, int h);
	~Linear_Queue();

public:
	void setParameters(double _binWidth, unsigned _binNum);

	void setPhase(unsigned _phase);

	void* push(const T &elem, double key);
	T top();
	T pop();
	void remove(void *ptr);

	int size();
	bool empty();

private:
	// 	std::vector< mylist<BinItem*> > bins;
	// 	std::vector<unsigned> binSizes;
	ListGroup<BinItem> bins;

	//geodesic::MemoryAllocator<BinItem> m_memory_bin;

	unsigned totalSize;

	unsigned firstBin;				//first bin that is not empty;
	unsigned lastBin;				//last bin that is not empty;
	double binWidth;

	unsigned phase;

	BinItem binItem;
	T tmpElem;
};

template <typename T>
Linear_Queue<T>::Linear_Queue()
{
	phase = 0;
	totalSize = 0;
}

template <typename T>
Linear_Queue<T>::Linear_Queue(int w, int h) : 
bins(w, h)
{
	phase = 0;
	totalSize = 0;
}

template <typename T>
Linear_Queue<T>::~Linear_Queue()
{

}

template <typename T>
void Linear_Queue<T>::setParameters(double _binWidth, unsigned _binNum)
{
	binWidth = _binWidth;

	bins.resize(_binNum);

	firstBin = -1;
	lastBin = -1;
}

template <typename T>
void Linear_Queue<T>::setPhase(unsigned _phase)
{
	phase = _phase;
}

template <typename T>
void* Linear_Queue<T>::push(const T &elem, double key)
{
	binItem.element = elem;
	binItem.index = (unsigned)(key / binWidth);

	void *rtPtr = bins.push_back(binItem, binItem.index);

	++totalSize;

	//update first & lastBin (if it need to be updated)
	if (lastBin == -1 || binItem.index >= lastBin) lastBin = binItem.index + 1;
	if (firstBin == -1 || binItem.index < firstBin) firstBin = binItem.index;

	return rtPtr;
}

template <typename T>
T Linear_Queue<T>::top()
{
	while (firstBin< lastBin && bins.empty(firstBin))
		++firstBin;

	void *headOfBinDummy = bins.begin(firstBin);
	return ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element.element;
}

template <typename T>
T Linear_Queue<T>::pop()
{
	while (firstBin< lastBin && bins.empty(firstBin))
		++firstBin;

	void *headOfBinDummy = bins.begin(firstBin);
	binItem = ((ListGroup<BinItem>::ListItem*)headOfBinDummy)->element;

	--totalSize;

	tmpElem = binItem.element;
	bins.erase(headOfBinDummy, binItem.index);

	return tmpElem;
}

template <typename T>
void Linear_Queue<T>::remove(void *ptr)
{
	binItem = ((ListGroup<BinItem>::ListItem*)ptr)->element;
	unsigned index = binItem.index;

	bins.erase(ptr, index);
	--totalSize;
}

template <typename T>
int Linear_Queue<T>::size()
{
	return totalSize;
}

template <typename T>
bool Linear_Queue<T>::empty()
{
	return totalSize == 0;
}

#endif
