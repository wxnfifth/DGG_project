#include "PriorityQueue.h"

using namespace MMP_HalfEdge;

bool operator==(const Window& w1, const Window& w2)
{
	return EQUALZERO(w1.b1-w2.b1) && EQUALZERO(w1.b2-w2.b2)
		&& EQUALZERO(w1.d1-w2.d1) && EQUALZERO(w1.d2-w2.d2)
		&& (w1.Edge_idx == w2.Edge_idx) && EQUALZERO(w1.sigma-w2.sigma)
		&& (w1.winType == w2.winType);
}

Priority_Queue::Priority_Queue()
{
	minHeap.clear();

	PriorityQueueItem qi;
	qi.key = -1;
	minHeap.push_back(qi);
}

Priority_Queue::~Priority_Queue()
{
	minHeap.clear();
}

void Priority_Queue::push(const Window& elem, const double& key, const list<Window>::iterator& ptr)
{
	PriorityQueueItem qi;
	qi.value = elem;
	qi.key = key;
	qi.ptr = ptr;

	minHeap.push_back(qi);
	int cur_idx = minHeap.size() - 1;
	int p_idx = cur_idx / 2;

	while(p_idx >= 0 && minHeap[p_idx].key > qi.key)
	{
		minHeap[cur_idx] = minHeap[p_idx];
		minHeap[cur_idx].ptr->idx_in_pq = cur_idx;

		cur_idx = p_idx;
		p_idx = cur_idx / 2;
	}
	minHeap[cur_idx] = qi;
	minHeap[cur_idx].ptr->idx_in_pq = cur_idx;
}

Window Priority_Queue::top()
{
	return minHeap[1].value;
}

double Priority_Queue::topKey()
{
	return minHeap[1].key;
}

Window Priority_Queue::pop()
{
	PriorityQueueItem rt = minHeap[1];
	rt.ptr->idx_in_pq = -1;

	int last_idx = minHeap.size() - 1;
	minHeap[1] = minHeap[last_idx];
	minHeap[1].ptr->idx_in_pq = 1;
	PriorityQueueItem qi = minHeap[1];
	minHeap.pop_back();

	if(minHeap.size() == 1)
		return rt.value;

	unsigned int cur_idx = 1;
	unsigned int lc = 2*cur_idx;
	unsigned int rc = lc+1;
	unsigned int minidx;
	if(rc < minHeap.size())										//有右叶子
		minidx = minHeap[lc].key < minHeap[rc].key ? lc : rc;
	else if(lc < minHeap.size())									//只有左叶子
		minidx = lc;
	else
		minidx = cur_idx;

	while(cur_idx != minidx && minHeap[minidx].key < qi.key)
	{
		minHeap[cur_idx] = minHeap[minidx];
		minHeap[cur_idx].ptr->idx_in_pq = cur_idx;

		cur_idx = minidx;
		lc = 2*cur_idx;
		rc = lc+1;
		if(rc < minHeap.size())
			minidx = minHeap[lc].key < minHeap[rc].key ? lc : rc;
		else if(lc < minHeap.size())	
			minidx = lc;
		else
			minidx = cur_idx;
	}
	minHeap[cur_idx] = qi;
	minHeap[cur_idx].ptr->idx_in_pq = cur_idx;

	return rt.value;
}

unsigned int Priority_Queue::find(const Window& elem)
{
	unsigned int idx = 1;
	unsigned int end_idx = minHeap.size();
	for(;idx < end_idx;++idx)
	{
		if(minHeap[idx].value == elem)
			break;
	}
	return idx;
}

void Priority_Queue::remove(const Window& elem)
{
	int rm_idx = find(elem);
	if(rm_idx == minHeap.size())								//找不到
		return;

	remove(rm_idx);
}

void Priority_Queue::remove(unsigned int rm_idx)
{
	int last_idx = minHeap.size() - 1;
	minHeap[rm_idx].ptr->idx_in_pq = -1;

	minHeap[rm_idx] = minHeap[last_idx];
	minHeap[rm_idx].ptr->idx_in_pq = rm_idx;

	PriorityQueueItem qi = minHeap[rm_idx];
	minHeap.pop_back();

	if(rm_idx == minHeap.size())
		return;

	unsigned int cur_idx = rm_idx;
	unsigned int lc = 2*cur_idx;
	unsigned int rc = lc+1;
	unsigned int p_idx = rm_idx/2;

	if(minHeap[cur_idx].key < minHeap[p_idx].key)			//如果最后一个(已移至rm_idx的位置)比其当前父亲还小，则应上滤
	{
		while(p_idx >= 0 && minHeap[p_idx].key > qi.key)
		{
			minHeap[cur_idx] = minHeap[p_idx];
			minHeap[cur_idx].ptr->idx_in_pq = cur_idx;
			cur_idx = p_idx;
			p_idx = cur_idx / 2;
		}
		minHeap[cur_idx] = qi;
		minHeap[cur_idx].ptr->idx_in_pq = cur_idx;
	}
	else
	{
		unsigned int minidx;
		if(rc < minHeap.size())										//有右叶子
			minidx = minHeap[lc].key < minHeap[rc].key ? lc : rc;
		else if(lc < minHeap.size())									//只有左叶子
			minidx = lc;
		else
			minidx = cur_idx;
		
		while(cur_idx != minidx && minHeap[minidx].key < qi.key)
		{
			minHeap[cur_idx] = minHeap[minidx];
			minHeap[cur_idx].ptr->idx_in_pq = cur_idx;
			cur_idx = minidx;
			lc = 2*cur_idx;
			rc = lc+1;
			if(rc < minHeap.size())
				minidx = minHeap[lc].key < minHeap[rc].key ? lc : rc;
			else if(lc < minHeap.size())	
				minidx = lc;
			else
				minidx = cur_idx;
		}
		minHeap[cur_idx] = qi;
		minHeap[cur_idx].ptr->idx_in_pq = cur_idx;
	}
}

unsigned int Priority_Queue::size()
{
	return minHeap.size();
}

bool Priority_Queue::empty()
{
	return minHeap.size() == 1;
}

void Priority_Queue::clear()
{
	minHeap.resize(1);
}