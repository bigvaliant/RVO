#include "navmesh.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
// HACK!!
#define inline static
#endif

inline float lerp(float a, float b, float t) { return a + (b-a)*t; }
inline float mini(int a, int b) { return a < b ? a : b; }
inline float maxi(int a, int b) { return a > b ? a : b; }
inline float minf(float a, float b) { return a < b ? a : b; }
inline float maxf(float a, float b) { return a > b ? a : b; }
inline float clamp(float a, float mn, float mx) { return a < mn ? mn : (a > mx ? mx : a); }
inline float sqr(float x) { return x*x; }
inline void vscale(float* v, const float* a, const float s) { v[0] = a[0]*s; v[1] = a[1]*s; }
inline void vset(float* a, const float x, const float y) { a[0]=x; a[1]=y; }
inline void vcpy(float* a, const float* b) { a[0]=b[0]; a[1]=b[1]; }
inline float vdot(const float* a, const float* b) { return a[0]*b[0] + a[1]*b[1]; }
inline float vperp(const float* a, const float* b) { return a[0]*b[1] - a[1]*b[0]; }
inline void vsub(float* v, const float* a, const float* b) { v[0] = a[0]-b[0]; v[1] = a[1]-b[1]; }
inline void vadd(float* v, const float* a, const float* b) { v[0] = a[0]+b[0]; v[1] = a[1]+b[1]; }

bool _less(heapele*l,heapele*r) { return ((mapnode*)l)->F < ((mapnode*)r)->F; }
void _clear(heapele*e){ ((mapnode*)e)->F = ((mapnode*)e)->G = ((mapnode*)e)->H = ((mapnode*)e)->z = 0; ((mapnode*)e)->parent = NULL; }

inline float vdistsqr(const float* a, const float* b)
{
	const float dx = b[0]-a[0];
	const float dy = b[1]-a[1];
	return dx*dx + dy*dy;
}

inline void vlerp(float* v, const float* a, const float* b, const float t)
{
	v[0] = a[0] + (b[0]-a[0])*t;
	v[1] = a[1] + (b[1]-a[1])*t;
}


inline int next(int i, int n) { return i+1 < n ? i+1 : 0; }

inline float triarea(const float* a, const float* b, const float* c)
{
	return (b[0]*a[1] - a[0]*b[1]) + (c[0]*b[1] - b[0]*c[1]) + (a[0]*c[1] - c[0]*a[1]);
}

inline int left(const float* a, const float* b, const float* c)
{
	const float EPS = 0.00001f;
	return triarea(a,b,c) < EPS;
}

inline int pntri(const float* a, const float* b, const float* c, const float* pt)
{
	return left(a, b, pt) && left(b, c, pt) && left(c, a, pt);
}


inline float distpt(const float* a, const float* b)
{
	const float dx = b[0] - a[0];
	const float dy = b[1] - a[1];
	return dx*dx + dy*dy;
}

void* poolAlloc( void* userData, unsigned int size )
{
	struct MemPool* pool = (struct MemPool*)userData;
	if (pool->size + size < pool->cap)
	{
		unsigned char* ptr = pool->buf + pool->size;
		pool->size += size;
		return ptr;
	}
	return 0;
}

bool requestLink(mapnode *tmp, float* va, float* vb, float* vc, float x1, float y1, float x2, float y2, int i)
{
	if (va[0] == x1 && va[1] == y1) {
		if (vb[0] == x2 && vb[1] == y2) {
			return true;
		} else if (vc[0] == x2 && vc[1] == y2) {
			return true;
		}
	} else if (vb[0] == x1 && vb[1] == y1) {
		if (va[0] == x2 && va[1] == y2) {
			return true;
		} else if (vc[0] == x2 && vc[1] == y2) {
			return true;
		}
	} else if (vc[0] == x1 && vc[1] == y1) {
		if (va[0] == x2 && va[1] == y2) {
			return true;
		} else if (vb[0] == x2 && vb[1] == y2) {
			return true;
		}
	}
	
	return false;
}

struct Edge
{
	unsigned short vert[2];
	unsigned short polyEdge[2];
	unsigned short poly[2];
};

static int buildadj(unsigned short* tris, const int ntris, const int nverts)
{
	int maxEdgeCount = ntris*3;
	unsigned short* firstEdge = 0;
	unsigned short* nextEdge = 0;
	struct Edge* edges = 0;
	int edgeCount = 0;
	int i,j;
	
	firstEdge = (unsigned short*)malloc(sizeof(unsigned short)*(nverts + maxEdgeCount));
	if (!firstEdge)
		goto cleanup;
	nextEdge = firstEdge + nverts;
	
	edges = (struct Edge*)malloc(sizeof(struct Edge)*maxEdgeCount);
	if (!edges)
		goto cleanup;
	
	for (i = 0; i < nverts; i++)
		firstEdge[i] = 0xffff;
	
	for (i = 0; i < ntris; ++i)
	{
		unsigned short* t = &tris[i*6];
		for (j = 0; j < 3; ++j)
		{
			unsigned short v0 = t[j];
			unsigned short v1 = t[(j+1) % 3];
			if (v0 < v1)
			{
				struct Edge* edge = &edges[edgeCount];
				edge->vert[0] = v0;
				edge->vert[1] = v1;
				edge->poly[0] = (unsigned short)i;
				edge->polyEdge[0] = (unsigned short)j;
				edge->poly[1] = (unsigned short)i;
				edge->polyEdge[1] = 0;
				// Insert edge
				nextEdge[edgeCount] = firstEdge[v0];
				firstEdge[v0] = (unsigned short)edgeCount;
				edgeCount++;
			}
		}
	}
	
	for (i = 0; i < ntris; ++i)
	{
		unsigned short* t = &tris[i*6];
		for (j = 0; j < 3; ++j)
		{
			unsigned short v0 = t[j];
			unsigned short v1 = t[(j+1) % 3];
			if (v0 > v1)
			{
				unsigned short e;
				for (e = firstEdge[v1]; e != 0xffff; e = nextEdge[e])
				{
					struct Edge* edge = &edges[e];
					if (edge->vert[1] == v0 && edge->poly[0] == edge->poly[1])
					{
						edge->poly[1] = (unsigned short)i;
						edge->polyEdge[1] = (unsigned short)j;
						break;
					}
				}
			}
		}
	}
	
	// Store adjacency

	for (i = 0; i < ntris; ++i)
	{
		unsigned short* t = &tris[i*6];
		t[3] = 0xffff;
		t[4] = 0xffff;
		t[5] = 0xffff;
	}

	for (i = 0; i < edgeCount; ++i)
	{
		struct Edge* e = &edges[i];
		if (e->poly[0] != e->poly[1])
		{
			unsigned short* t0 = &tris[e->poly[0]*6];
			unsigned short* t1 = &tris[e->poly[1]*6];
			t0[3+e->polyEdge[0]] = e->poly[1];
			t1[3+e->polyEdge[1]] = e->poly[0];
		}
	}

	free(firstEdge);
	free(edges);
	return 1;
	
cleanup:
	if (firstEdge)
		free(firstEdge);
	if (edges)
		free(edges);
	return 0;
}

static int buildmapnode(struct Navmesh* nav)
{
	if (!nav) return 0;

	//分配map_node内存
	nav->defaultMap = (mapnode*)calloc(nav->ntris*2,sizeof(*nav->defaultMap));
	if (!nav->defaultMap) return 0;

	for (int i = 0; i < nav->ntris; ++i)
	{
		float cx1, cy1, cx2, cy2, cx3, cy3;
		unsigned short* tri = &nav->tris[i*6];
		float* va = &nav->verts[tri[0]*2];
		float* vb = &nav->verts[tri[1]*2];
		float* vc = &nav->verts[tri[2]*2];
		mapnode *tmp = &nav->defaultMap[i];
		tmp->i  = i;
		
		tmp->links[0] = -1;
		tmp->links[1] = -1;
		tmp->links[2] = -1;

		//重心
		tmp->x  = (va[0] + vb[0] + vc[0])/3;
		tmp->y  = (va[1] + vb[1] + vc[1])/3;
		
		//中心点坐标
		cx1 = (va[0] + vb[0])/2;
		cy1 = (va[1] + vb[1])/2;
		cx2 = (vc[0] + vb[0])/2;
		cy2 = (vc[1] + vb[1])/2;
		cx3 = (vc[0] + va[0])/2;
		cy3 = (vc[1] + va[1])/2;

		//三条边中心点距离
		tmp->v[0] = sqrt((cx1 - cx2)*(cx1 - cx2) + (cy1 - cy2)*(cy1 - cy2));
		tmp->v[1] = sqrt((cx2 - cx3)*(cx2 - cx3) + (cy2 - cy3)*(cy2 - cy3));
		tmp->v[2] = sqrt((cx1 - cx3)*(cx1 - cx3) + (cy1 - cy3)*(cy1 - cy3));

		//邻接三角形
		for (int i = 0; i < 3; ++i)
		{
			const unsigned short nei = tri[3+i];
			if (nei != 0xffff)
			{
				unsigned short* tri_nei = &nav->tris[nei*6];
				float* na = &nav->verts[tri_nei[0]*2];
				float* nb = &nav->verts[tri_nei[1]*2];
				float* nc = &nav->verts[tri_nei[2]*2];

				//mapnode *nei_node = &nav->map[nei];	
				if (requestLink(tmp, na, nb, nc, va[0], va[1], vb[0], vb[1], nei)) { tmp->links[0] = nei; continue; }
				if (requestLink(tmp, na, nb, nc, vb[0], vb[1], vc[0], vc[1], nei)) { tmp->links[1] = nei; continue; }
				if (requestLink(tmp, na, nb, nc, vc[0], vc[1], va[0], va[1], nei)) { tmp->links[2] = nei; continue; }
			}
		}
	}

	return 1;
}


static float closestPtPtSeg(const float* pt, const float* sp, const float* sq)
{
	float dir[2],diff[3],t,d;
	vsub(dir,sq,sp);
	vsub(diff,pt,sp);
	t = vdot(diff,dir);
	if (t <= 0.0f) return 0;
	d = vdot(dir,dir);
	if (t >= d) return 1;
	return t/d;
}

// navmeshCreateEx
struct Navmesh* navmeshCreateEx(std::vector<int> const& trisList, std::vector<std::pair<float, float> > const& pointList)
{
	int i=0,j=0,z=0;
	struct Navmesh* nav = 0;
	int npts = pointList.size();
	int ntrisL = trisList.size();
	
	nav = (struct Navmesh*)malloc(sizeof(struct Navmesh));
	if (!nav)
		goto cleanup;
	memset(nav, 0, sizeof(struct Navmesh));

	nav->verts = (float*)malloc(sizeof(float)*npts*2);
	if (!nav->verts)
		goto cleanup;
	for(std::vector<std::pair<float, float> >::size_type ix = 0; ix != npts; ++ix)
	{
		nav->verts[i] = pointList[ix].first;
		nav->verts[i+1] = pointList[ix].second;
		i=i+2;
	}
	nav->nverts = npts;

	nav->tris = (unsigned short*)malloc(sizeof(unsigned short)*(ntrisL*2));
	if (!nav->tris)
		goto cleanup;
	for(std::vector<int>::size_type jx = 0; jx != ntrisL; ++jx)
	{
		nav->tris[j] = (unsigned short)trisList[jx];
		++z;
		if(z==3) 
		{
			j=j+4;
			z=0;
		}
		else	
		{
			++j;
		}
	}
	nav->ntris = ntrisL/3;

	if (nav->ntris < 0) nav->ntris = -nav->ntris;
	if (!nav->ntris)
		goto cleanup;

	if (!buildadj(nav->tris, nav->ntris, nav->nverts))
		goto cleanup;

	if (!buildmapnode(nav))
		goto cleanup;
	
	nav->close_list = new dlist();
	nav->open_list = new minheap(8192,_less,_clear);
		
	return nav;
	
cleanup:
	if (nav)
	{
		if (nav->verts)
			free(nav->verts);
		if (nav->tris)
			free(nav->tris);
		if (nav->defaultMap)
			free(nav->defaultMap);
		if (nav->close_list)
			free(nav->close_list);
		if (nav->open_list)
			free(nav->open_list);
		free(nav);
	}

	return 0;
}


// Find nearest triangle
int navmeshFindNearestTri(struct Navmesh* nav, const float* pos, float* nearest)
{
	int i, j, besti = -1;
	float t, d, p[2], bestd = FLT_MAX;

	if (nearest)
		vcpy(nearest, pos);

	for (i = 0; i < nav->ntris; ++i)
	{
		const unsigned short* tri = &nav->tris[i*6];
		const float* va = &nav->verts[tri[0]*2];
		const float* vb = &nav->verts[tri[1]*2];
		const float* vc = &nav->verts[tri[2]*2];
		if (pntri(va,vb,vc,pos))
		{	
			return i;
		}
	}

	for (i = 0; i < nav->ntris; ++i)
	{
		const unsigned short* tri = &nav->tris[i*6];
		for (j = 0; j < 3; ++j)
		{
			const float* va = &nav->verts[tri[j]*2];
			const float* vb = &nav->verts[tri[(j+1)%3]*2];
			if (tri[3+j] != 0xffff)
				continue;
			t = closestPtPtSeg(pos,va,vb);
			vlerp(p,va,vb,t);
			d = distpt(p,pos);
			if (d < bestd)
			{
				if (nearest)
					vcpy(nearest, p);
				bestd = d;
				besti = i;
			}
		}
	}
	
	return besti;
}

int navmeshFindPath(struct Navmesh* nav, const float* start, const float* end, unsigned short* path, const int maxpath)
{
// 根据地图大小适当调整宏
#define MAX_STACK 128*8
#define MAX_PARENT 128*8
	int i, starti, endi, stack[MAX_STACK], nstack;
	unsigned short parent[MAX_PARENT];
	
	starti = navmeshFindNearestTri(nav, start, NULL);
	endi = navmeshFindNearestTri(nav, end, NULL);
	if (starti == -1 || endi == -1)
		return 0;
		
	if (starti == endi)
	{
		path[0] = (unsigned short)starti;
		return 1;
	}

	memset(parent, 0xff, sizeof(unsigned short)*MAX_PARENT);
		
	nstack = 0;
	stack[nstack++] = endi;
	parent[endi] = endi;

	while (nstack)
	{
		unsigned short* tri;
		unsigned short cur;
		
		// Pop front.
		cur = stack[0];
		nstack--;
		for (i = 0; i < nstack; ++i)
			stack[i] = stack[i+1];

		if (cur == starti)
		{
			// Trace and store back.
			int npath = 0;
			for (;;)
			{
				path[npath++] = cur;
				if (npath >= maxpath) break;
				if (parent[cur] == cur) break;
				cur = parent[cur];
			}
			return npath;
		}
		
		tri = &nav->tris[cur*6];
		for (i = 0; i < 3; ++i)
		{
			const unsigned short nei = tri[3+i];
			if (nei == 0xffff) continue;
			if (parent[nei] != 0xffff) continue;
			parent[nei] = cur;
			if (nstack < MAX_STACK)
				stack[nstack++] = nei;
		}
	}
	
	return 0;
}

void reset(minheap* open_list, dlist* close_list)
{
	mapnode *n = NULL;
	while(n = (mapnode*)close_list->Pop()){
		n->G = n->H = n->F = n->z = 0;
		n->parent = NULL;
	}
	open_list->Clear();
}

//计算到达相临节点需要的代价
double cost_2_neighbor(mapnode *from,int i)
{
	int zz = i + from->z;
	if (zz == 0 || zz == 1)
		return from->v[0];
	else if(zz == 2)
		return from->v[2];
	else if(zz == 3)
		return from->v[1];
	else
	{
		assert(0);
		return 0.0f;
	}
}

// H估值计算
double cost_2_goal(mapnode *from,const float* to)
{
	float delta_x = from->x - to[0];
	float delta_y = from->y - to[1];
	return sqrt(delta_x*delta_x + delta_y*delta_y);
}

//记录路径从上一个节点进入该节点的边（如果从终点开始寻路即为穿出边）
void setZ(mapnode *node, int idx)
{	
	node->z = 0;
	if (node->links[0] == idx)
		node->z = 0;
	else if(node->links[1] == idx)
		node->z = 1;
	else if(node->links[2] == idx)
		node->z = 2;
	else
		assert(0);
}

int navmeshFindBestPath(struct Navmesh* nav, const float* start, const float* end, unsigned short* path)
{
	int starti, endi, tpath[AGENT_MAX_PATH];
	starti = navmeshFindNearestTri(nav, start, NULL);
	endi = navmeshFindNearestTri(nav, end, NULL);
	if (starti == -1 || endi == -1)
		return 0;
		
	if (starti == endi)
	{
		path[0] = (unsigned short)starti;
		return 1;
	}
	
	mapnode *from = &nav->defaultMap[starti];
	mapnode *to = &nav->defaultMap[endi];
		
	nav->open_list->insert(to);

	mapnode *current_node = NULL;
	for(;;)
	{
		current_node = (mapnode*)nav->open_list->popmin();
		if(!current_node){ 
			reset(nav->open_list, nav->close_list);
			return 0;
		}

		if(current_node->i == from->i){
			int npath = 0;
			while(current_node)
			{
				path[npath++] = current_node->i;
				mapnode *t = current_node;
				current_node = current_node->parent;
				t->parent = NULL;
				t->F = t->G = t->H = t->z = 0;
				t->index = 0;
			}

			reset(nav->open_list, nav->close_list);
			return npath;
		}

		//current插入到close表
		nav->close_list->Push(current_node);

		for(int i=0; i<3; ++i)
		{
			int nei = current_node->links[i];
			if (nei == -1) continue;

			mapnode *neighbor = &nav->defaultMap[nei];
			if(neighbor->pre || neighbor->next){
				continue;//在close表中,不做处理
			}

			if(neighbor->index)//在openlist中
			{
				double new_G = current_node->G + cost_2_neighbor(current_node,i);
				if(new_G < neighbor->G)
				{
					//经过当前neighbor路径更佳,更新路径
					neighbor->G = new_G;
					neighbor->F = neighbor->G + neighbor->H;
					neighbor->parent = current_node;
					nav->open_list->change(neighbor);
					setZ(neighbor, current_node->i);
				}
				continue;
			}
			neighbor->parent = current_node;
			neighbor->G = current_node->G + cost_2_neighbor(current_node,i);
			neighbor->H = cost_2_goal(neighbor,start);
			neighbor->F = neighbor->G + neighbor->H;
			nav->open_list->insert(neighbor);
			setZ(neighbor, current_node->i);
		}
	}
	
	return 0;
}

inline int vequal(const float* a, const float* b)
{
	static const float eq = 0.001f*0.001f;
	return distpt(a, b) < eq;
}

inline int pushPoint(float* pts, int npts, const int maxpts, const float* pt)
{
	if (npts >= maxpts) return npts;
	if (vequal(pt, &pts[(npts-1)*2])) return npts;
	vcpy(&pts[npts*2], pt);
	return npts+1;
}

static int stringPull(const float* portals, int nportals, float* pts, const int maxpts)
{
	int npts = 0, i;
	float portalApex[2], portalLeft[2], portalRight[2];
	int apexIndex = 0, leftIndex = 0, rightIndex = 0;
	vcpy(portalApex, &portals[0]);
	vcpy(portalLeft, &portals[0]);
	vcpy(portalRight, &portals[2]);
	
	// Add start point.
	vcpy(&pts[npts*2], portalApex);
	npts++;
	
	for (i = 1; i < nportals && npts < maxpts; ++i)
	{
		const float* left = &portals[i*4+0];
		const float* right = &portals[i*4+2];
		
		// Update right vertex.
        if (triarea(portalApex, portalRight, right) <= 0.0f)
		{
			if (vequal(portalApex, portalRight) || triarea(portalApex, portalLeft, right) > 0.0f)
			{
				// Tighten the funnel.
				vcpy(portalRight, right);
				rightIndex = i;
			}
			else
			{
				// Right over left, insert left to path and restart scan from portal left point.
				npts = pushPoint(pts, npts, maxpts, portalLeft);
				// Make current left the new apex.
				vcpy(portalApex, portalLeft);
				apexIndex = leftIndex;
				// Reset portal
				vcpy(portalLeft, portalApex);
				vcpy(portalRight, portalApex);
				leftIndex = apexIndex;
				rightIndex = apexIndex;
				// Restart scan
				i = apexIndex;
				continue;
			}
		}
		
		// Update left vertex.
        if (triarea(portalApex, portalLeft, left) >= 0.0f)
		{
			if (vequal(portalApex, portalLeft) || triarea(portalApex, portalRight, left) < 0.0f)
			{
				// Tighten the funnel.
				vcpy(portalLeft, left);
				leftIndex = i;
			}
			else
			{
				// Left over right, insert right to path and restart scan from portal right point.
				npts = pushPoint(pts, npts, maxpts, portalRight);
				// Make current right the new apex.
				vcpy(portalApex, portalRight);
				apexIndex = rightIndex;
				// Reset portal
				vcpy(portalLeft, portalApex);
				vcpy(portalRight, portalApex);
				leftIndex = apexIndex;
				rightIndex = apexIndex;
				// Restart scan
				i = apexIndex;
				continue;
			}
		}
	}
	// Append last point to path.
	npts = pushPoint(pts, npts, maxpts, &portals[(nportals-1)*4+0]);
	
	return npts;
}

static void getPortalPoints(struct Navmesh* nav, const unsigned short a, const unsigned short b,
							float* left, float* right)
{
	const unsigned short* ta = &nav->tris[a*6];
	int i;
	for (i = 0; i < 3; ++i)
	{
		if (ta[3+i] == b)
		{
			const float* va = &nav->verts[ta[i]*2];
			const float* vb = &nav->verts[ta[(i+1)%3]*2];
			vcpy(right, va);
			vcpy(left, vb);
		}
	}
}

int navmeshStringPull(struct Navmesh* nav, const float* start, const float* end,
					  const unsigned short* path, const int npath,
					  float* pts, const int maxpts)
{
#define MAX_PORTALS 128
	float portals[MAX_PORTALS*4];
	int nportals = 0, i;
	
	// Start portal
	vcpy(&portals[nportals*4+0], start);
	vcpy(&portals[nportals*4+2], start);
	nportals++;
	// Portal between navmesh polygons
	for (i = 0; i < npath-1; ++i)
	{
		getPortalPoints(nav, path[i], path[i+1], &portals[nportals*4+0], &portals[nportals*4+2]);
		nportals++;
	}
	// End portal
	vcpy(&portals[nportals*4+0], end);
	vcpy(&portals[nportals*4+2], end);
	nportals++;
	
	return stringPull(portals, nportals, pts, maxpts);
}

// Deletes navmesh.
void navmeshDelete(struct Navmesh* nav)
{
	if (!nav)
		return;
	if (nav->verts)
		free(nav->verts);
	if (nav->tris)
		free(nav->tris);
	if (nav->defaultMap)
		free(nav->defaultMap);
	if (nav->close_list)
		free(nav->close_list);
	if (nav->open_list)
		free(nav->open_list);
	free(nav);
}

int posValid(const float* p)
{
	return p[0] < FLT_MAX && p[1] < FLT_MAX;
}

void agentInit(struct NavmeshAgent* agent)
{
	vset(agent->pos, FLT_MAX, FLT_MAX);
	vset(agent->oldpos, FLT_MAX, FLT_MAX);
	vset(agent->target, FLT_MAX, FLT_MAX);
	agent->npath = 0;
	agent->ncorners = 0;
}

void agentFindPath(struct NavmeshAgent* agent, struct Navmesh* nav)
{
	agent->npath = 0;

	if (posValid(agent->target))
		agent->npath = navmeshFindBestPath(nav, agent->pos, agent->target, agent->path);
	    //agent->npath = navmeshFindPath(nav, agent->pos, agent->target, agent->path, AGENT_MAX_PATH);
}


int FindPath(struct Navmesh* nav, struct NavmeshAgent* agent)
{	

	//获得联通的三角形
	agentFindPath(agent, nav);

	//拉线获得路径
	agent->ncorners = 0;
	if (agent->npath)
		agent->ncorners = navmeshStringPull(nav, agent->pos, agent->target, agent->path, agent->npath, agent->corners, MAX_CORNERS);
	return agent->ncorners;
}