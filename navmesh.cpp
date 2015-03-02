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

void parsefloat(std::vector<float> const& f, std::vector<std::pair<float, float>> &vector_point,std::vector<int> &vec_idx)
{
	std::map<int, std::pair<float, float>> map_idx_point;
	std::map<std::pair<float, float>, int> map_point_idx;
	std::map<int, int> map_idx1_idx2;
	std::list<int> list_idx;

	std::map<std::pair<float, float>, int>::iterator itr1;
	std::map<int, std::pair<float, float>>::iterator itr2;
	std::map<int, int>::iterator itr3;
	std::list<int>::iterator itr4;
	
	int n = f.size();
	for (int i=0;i<n/2;i++)
	{
		std::pair<float, float> value;
		value = std::make_pair(f[i*2],f[i*2+1]);
		vec_idx.push_back(i);
		map_idx_point.insert(std::make_pair(i, value));
		itr1 = map_point_idx.find(value);
		map_point_idx.insert(std::make_pair(value, i));
		if(itr1 != map_point_idx.end())
		{		
			map_idx1_idx2.insert(std::make_pair(i,itr1->second));
		}
	}
	
	itr1 = map_point_idx.begin();
	for(; itr1 != map_point_idx.end(); ++itr1 )
	{
		itr2 = map_idx_point.begin();
		for(; itr2 != map_idx_point.end(); ++itr2 )
		{
			if(itr2->second.first == itr1->first.first && itr2->second.second == itr1->first.second && itr1->second != itr2->first)
			{
				for(std::vector<int>::size_type ix = 0; ix != vec_idx.size(); ++ix)
				{
					if(vec_idx[ix] == itr2->first)
					{
						vec_idx[ix] = itr1->second;
					}
				}
			}
		}
	}
	
	itr1 = map_point_idx.begin();
	for(; itr1 != map_point_idx.end(); ++itr1 )
	{
		list_idx.push_back(itr1->second);
	}
	list_idx.sort();
	
	itr4 = list_idx.begin();
	for(; itr4 != list_idx.end(); ++itr4 )
	{
		itr2 = map_idx_point.find(*itr4);
		vector_point.push_back(itr2->second);
	}

	for(std::vector<int>::size_type ix = 0; ix != vec_idx.size(); ++ix)
	{
		itr3 = map_idx1_idx2.find(vec_idx[ix]);
		if(itr3 != map_idx1_idx2.end())
		{
			vec_idx[ix] = itr3->second;
		}
	}

	for(std::vector<std::pair<float, float>>::size_type ix = 0; ix != vector_point.size(); ++ix)
	{
		itr1 = map_point_idx.find(vector_point[ix]);
		for(std::vector<int>::size_type jx = 0; jx != vec_idx.size(); ++jx)
		{
			if(vec_idx[jx]==itr1->second)
			{
				vec_idx[jx]=ix;
			}
		}
	}

	//for(std::vector<int>::size_type jx = 0; jx != vec_idx.size(); ++jx)
	//{
	//	printf("%f, %f\n", vector_point[vec_idx[jx]].first, vector_point[vec_idx[jx]].second);
	//}
}

// navmeshCreateEx
struct Navmesh* navmeshCreateEx(std::vector<float> const& f)
{

	std::vector<int> trisList;
	std::vector<std::pair<float, float>> pointList;

	parsefloat(f, pointList, trisList);

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
	for(std::vector<std::pair<float, float>>::size_type ix = 0; ix != npts; ++ix)
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

	//for(int k=0;k<(npts*2);k++)
	//{
	//	printf("k:%d, verts:%f\n", k, nav->verts[k]);
	//}

	//for(int k=0;k<ntrisL*2;k++)
	//{
	//	if(nav->tris[k] != 52685)
	//		printf("k:%d, tris:%d\n", k, nav->tris[k]);
	//}

	if (!buildadj(nav->tris, nav->ntris, nav->nverts))
		goto cleanup;

	return nav;
	
cleanup:
	if (nav)
	{
		if (nav->verts)
			free(nav->verts);
		if (nav->tris)
			free(nav->tris);
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
			//printf("%f:%f %f:%f %f:%f\n", va[0],va[1], vb[0],vb[1], vc[0],vc[1]);
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

	//const unsigned short* tri;
	//const float* va;
	//const float* vb;
	//const float* vc;
	
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
				
				/*tri = &nav->tris[cur*6];
				va = &nav->verts[tri[0]*2];
				vb = &nav->verts[tri[1]*2];
				vc = &nav->verts[tri[2]*2];
				printf("%f:%f %f:%f %f:%f\n", va[0],va[1], vb[0],vb[1], vc[0],vc[1]);*/

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
	agent->nvisited = 0;
	agent->ncorners = 0;
}

void agentFindPath(struct NavmeshAgent* agent, struct Navmesh* nav)
{
	agent->npath = 0;

	if (posValid(agent->target))
		agent->npath = navmeshFindPath(nav, agent->pos, agent->target, agent->path, AGENT_MAX_PATH);
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
