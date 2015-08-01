#ifndef NAVMESH_H
#define NAVMESH_H

#include <vector>
#include <map>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

class dlist;
struct dnode{
	dnode():next(NULL),pre(NULL),_dlist(NULL){}
	dnode *next;
	dnode *pre;
	dlist *_dlist;
};

class dlist{
public:
	dlist():size(0){
		head.next = &tail;
		tail.pre = &head;
	}
	size_t Size(){
		return size;
	}

	bool Empty(){
		return size == 0;
	}

	dnode *Begin(){
		if(Empty())return &tail;
		return head.next;
	}

	dnode *End(){
		return &tail;
	}

	void Push(dnode *n){
		if(n->_dlist || n->next || n->pre) return;
		tail.pre->next = n;
		n->pre = tail.pre;
		tail.pre = n;
		n->_dlist = this;
		n->next = &tail;
		++size;
	}

	void Remove(dnode *n)
	{
		if(n->_dlist != this || (!n->pre && !n->next))
			return;
		n->pre->next = n->next;
		n->next->pre = n->pre;
		n->pre = n->next = NULL;
		n->_dlist = NULL;
		--size;
	}

	dnode* Pop(){
		if(Empty())
			return NULL;
		else
		{
			dnode *n = head.next;
			Remove(n);
			return n;
		}
	}

private:
	dnode head;
	dnode tail;
	size_t size;
};

struct heapele
{
	int index;//index in minheap;
};

struct mapnode : public heapele,public dnode
{
	mapnode *parent;
	double G;//从初始点到当前点的开销
	double H;//从当前点到目标点的估计开销
	double F;
	float  x;
	float  y;
	int    i;
	double v[3];
	int links[3];
	int    z;
};

struct minheap{
private:
	int size;
	int max_size;
	typedef bool (*less)(struct heapele*l,struct heapele*r);
	typedef void (*clear)(struct heapele*); 
	less cmp_function;
	clear clear_function;
	heapele **elements;

private:
	void swap(int idx1,int idx2)
	{
		heapele *ele = elements[idx1];
		elements[idx1] = elements[idx2];
		elements[idx2] = ele;
		elements[idx1]->index = idx1;
		elements[idx2]->index = idx2;
	}

	int parent(int idx){
		return idx/2;
	}

	int left(int idx){
		return idx*2;
	}

	inline int right(int idx){
		return idx*2+1;
	}

	void up(int idx){
		int p = parent(idx);
		while(p)
		{
			assert(elements[idx]);
			assert(elements[p]);
			if(cmp_function(elements[idx],elements[p]))
			{
				swap(idx,p);
				idx = p;
				p = parent(idx);
			}
			else
				break;
		}
	}

	void down(int idx){

		int l = left(idx);
		int r = right(idx);
		int min = idx;
		if(l <= size)
		{
			assert(elements[l]);
			assert(elements[idx]);
		}

		if(l <= size && cmp_function(elements[l],elements[idx]))
			min = l;

		if(r <= size)
		{
			assert(elements[r]);
			assert(elements[min]);
		}

		if(r <= size && cmp_function(elements[r],elements[min]))
			min = r;

		if(min != idx)
		{
			swap(idx,min);
			down(min);
		}		
	}


public:
	minheap(int size,less cmp_func,clear clear_func):size(0),max_size(size),cmp_function(cmp_func),clear_function(clear_func){
		elements = (heapele **)calloc(size,sizeof(*elements));
	}

	~minheap(){
		free(elements);
	}

	void change(heapele *e)
	{
		int idx = e->index;
		down(idx);
		if(idx == e->index)
			up(idx);
	}

	void insert(heapele *e)
	{
		if(e->index)
			return change(e);
		if(size >= max_size-1)
		{
			//expand the heap
			unsigned int new_size = max_size*2;
			heapele** tmp = (heapele**)calloc(new_size,sizeof(*tmp));
			if(!tmp) return;
			memcpy(tmp,elements,max_size*sizeof(*tmp));
			free(elements);
			elements = tmp;
			max_size = new_size;
		}	
		++size;
		elements[size] = e;
		e->index = size;
		up(e->index);	
	}

	/*void remove(heapele *e)
	{
		heapele **back = (heapele**)calloc(1,sizeof(*back)*(size-1));
		int i = 1;
		int c = 1;
		for( ; i <= size;++i)
		{
			elements[i]->index = 0;
			if(elements[i] != e)
				back[c++] = elements[i];
		}
		memset(&(elements[0]),0,sizeof(*back)*max_size);
		size = 0;
		i = 1;
		for(; i < c; ++i) insert(back[i]);
		free(back);
	}*/

	heapele* minheap_min()
	{
		if(!size)	return NULL;
		return elements[1];
	}

	//return the min element and remove it
	heapele* popmin()
	{
		if(size)
		{
			heapele *e = elements[1];
			swap(1,size);
			elements[size] = NULL;
			--size;
			down(1);
			e->index = 0;
			return e;
		}
		return NULL;
	}

	void Clear()
	{
		for(int i = 1; i <= size; ++i){
			elements[i]->index = 0;
			clear_function(elements[i]);
		}
		size = 0;
	}

	int Size()
	{
		return size;
	}
};

struct Navmesh
{
	float* verts;
	int nverts;
	unsigned short* tris;
	int ntris;
	mapnode* defaultMap;
	dlist* close_list;
	minheap* open_list;
};

// Creates navmesh from a polygon.
struct Navmesh* navmeshCreateEx(std::vector<int> const& trisList, std::vector<std::pair<float, float> > const& pointList);

// Find nearest triangle
int navmeshFindNearestTri(struct Navmesh* nav, const float* pos, float* nearest);

// Find path
int navmeshFindPath(struct Navmesh* nav, const float* start, const float* end, unsigned short* path, const int maxpath);

// optimal Find path
int navmeshFindBestPath(struct Navmesh* nav, const float* start, const float* end, unsigned short* path);

// Find tight rope path
int navmeshStringPull(struct Navmesh* nav, const float* start, const float* end,
					  const unsigned short* path, const int npath,
					  float* pts, const int maxpts);

// Deletes navmesh.
void navmeshDelete(struct Navmesh* nav);


#define MAX_CORNERS 64
#define AGENT_MAX_TRAIL 64
#define AGENT_MAX_PATH 128

#define VEL_HIST_SIZE 6

struct NavmeshAgent
{
	float pos[2];
	float target[2];
	float oldpos[2];
	float corners[MAX_CORNERS];
	int ncorners;
	unsigned short path[AGENT_MAX_PATH];
	int npath;
};
    
struct MemPool
{
    unsigned char* buf;
    unsigned int cap;
    unsigned int size;
};

//对外查找接口
int FindPath(struct Navmesh* nav, struct NavmeshAgent* agent);
void agentInit(struct NavmeshAgent* agent);

#ifdef __cplusplus
}
#endif
		

#endif