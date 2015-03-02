#ifndef NAVMESH_H
#define NAVMESH_H

#include <vector>
#include <map>
#include <list>

#ifdef __cplusplus
extern "C" {
#endif

struct Navmesh
{
	float* verts;
	int nverts;
	unsigned short* tris;
	int ntris;
};

// Creates navmesh from a polygon.
struct Navmesh* navmeshCreateEx(std::vector<float> const& f);

// Find nearest triangle
int navmeshFindNearestTri(struct Navmesh* nav, const float* pos, float* nearest);

// Find path
int navmeshFindPath(struct Navmesh* nav, const float* start, const float* end, unsigned short* path, const int maxpath);

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
	unsigned short visited[AGENT_MAX_PATH];
	int nvisited;
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
