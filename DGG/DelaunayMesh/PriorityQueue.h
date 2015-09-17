#ifndef PRIORITYQUEUE_H1
#define PRIORITYQUEUE_H1

#include <vector>
#include <list>
#include <bitset>
#include "Geometry.h"

#define MAXSOURCENUM 10

using namespace std;

enum Window_Type
{POINTWINDOW, LINEWINDOW};

namespace MMP_HalfEdge
{
	//���ڽṹ
	struct Window
	{
		//�������ڵı�
		unsigned int Edge_idx;
		//��Ԫ���ʾ�Ĳ�غ���
		double b1, b2;
		double d1, d2;
		double sigma;
		//��ʾαԴ�����ĸ���
		bool tao;
		//���ڵ�����
		Window_Type winType;
		//���������ȶ����е��±�
		int idx_in_pq;

		unsigned int srcId;
		unsigned pseuSrcId;
	};

	class Priority_Queue
	{
	public:
		Priority_Queue();
		~Priority_Queue();

	public:
		struct PriorityQueueItem
		{
			Window value;
			double key;

			list<Window>::iterator ptr;
		};
		void push(const Window& elem, const double& key, const list<Window>::iterator& ptr);
		Window top();
		double topKey();
		Window pop();

		unsigned int find(const Window& elem);
		void remove(const Window& elem);
		void remove(unsigned int rm_idx);

		unsigned int size();
		bool empty();

		void clear();

	private:
		vector<PriorityQueueItem> minHeap;
	};

}

#endif
