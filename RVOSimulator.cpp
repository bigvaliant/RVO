/*
 * RVOSimulator.cpp
 * RVO2 Library
 *
 * Copyright (c) 2008-2010 University of North Carolina at Chapel Hill.
 * All rights reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and non-profit purposes, without
 * fee, and without a written agreement is hereby granted, provided that the
 * above copyright notice, this paragraph, and the following four paragraphs
 * appear in all copies.
 *
 * Permission to incorporate this software into commercial products may be
 * obtained by contacting the Office of Technology Development at the University
 * of North Carolina at Chapel Hill <otd@unc.edu>.
 *
 * This software program and documentation are copyrighted by the University of
 * North Carolina at Chapel Hill. The software program and documentation are
 * supplied "as is," without any accompanying services from the University of
 * North Carolina at Chapel Hill or the authors. The University of North
 * Carolina at Chapel Hill and the authors do not warrant that the operation of
 * the program will be uninterrupted or error-free. The end-user understands
 * that the program was developed for research purposes and is advised not to
 * rely exclusively on the program for any reason.
 *
 * IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL OR THE
 * AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
 * CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
 * SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF NORTH CAROLINA AT
 * CHAPEL HILL OR THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 * THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL AND THE AUTHORS SPECIFICALLY
 * DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY
 * STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON
 * AN "AS IS" BASIS, AND THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL AND THE
 * AUTHORS HAVE NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 * ENHANCEMENTS, OR MODIFICATIONS.
 *
 * Please send all bug reports to <geom@cs.unc.edu>.
 *
 * The authors may be contacted via:
 *
 * Jur van den Berg, Stephen J. Guy, Jamie Snape, Ming C. Lin, Dinesh Manocha
 * Dept. of Computer Science
 * 201 S. Columbia St.
 * Frederick P. Brooks, Jr. Computer Science Bldg.
 * Chapel Hill, N.C. 27599-3175
 * United States of America
 *
 * <http://gamma.cs.unc.edu/RVO2/>
 */

#include "RVOSimulator.h"

#include "Agent.h"
#include "KdTree.h"
#include "Obstacle.h"
#include "navmesh.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace RVO {
	RVOSimulator::RVOSimulator() : defaultAgent_(NULL), globalTime_(0.0f), kdTree_(NULL), timeStep_(0.0f), count_(0), obstcount_(0), isBuildObstFlag_(false)
	{
		kdTree_ = new KdTree(this);
		preStepTime_ = 0;
		luaThreadIndex_ = 0; 
		instanceId_ = 0;

		sample_ = 0;
	}

	RVOSimulator::RVOSimulator(float timeStep, float neighborDist, size_t maxNeighbors, float timeHorizon, float timeHorizonObst, float radius, float maxSpeed, const Vector2 &velocity) : defaultAgent_(NULL), globalTime_(0.0f), kdTree_(NULL), timeStep_(timeStep), count_(0), obstcount_(0), isBuildObstFlag_(false)
	{
		preStepTime_ = 0;
		luaThreadIndex_ = 0;

		sample_ = 0;

		kdTree_ = new KdTree(this);
		defaultAgent_ = new Agent(this);

		defaultAgent_->maxNeighbors_ = maxNeighbors;
		defaultAgent_->maxSpeed_ = maxSpeed;
		defaultAgent_->neighborDist_ = neighborDist;
		defaultAgent_->radius_ = radius;
		defaultAgent_->timeHorizon_ = timeHorizon;
		defaultAgent_->timeHorizonObst_ = timeHorizonObst;
		defaultAgent_->velocity_ = velocity;
	}

	RVOSimulator::~RVOSimulator()
	{
		if (defaultAgent_ != NULL) {
			delete defaultAgent_;
		}

		for (size_t i = 0; i < agents_.size(); ++i) {
			delete agents_[i];
		}
		
		//for (size_t i = 0; i < obstacles_.size(); ++i) {
		//	delete obstacles_[i];
		//}

		navmeshDelete(nav_);

		delete sample_;

		portNoMap_.clear();
		delete kdTree_;
	}

	void RVOSimulator::setNavTargetPos(Agent *agent, const float target[])
	{
		if (agent) {
			vcpy(agent->navAgent_.target, target);
		}
	}

	int RVOSimulator::findPath(Agent *agent)
	{
		if (agent) {
			agent->isTrailOverFlag_ = false;
			agent->trailList_.clear();
			FindPath(nav_, &agent->navAgent_);
			for(int j=0; j<agent->navAgent_.ncorners*2;) {
				if (agent->navAgent_.corners[j] != agent->navAgent_.pos[0] || agent->navAgent_.corners[j+1] != agent->navAgent_.pos[1]) {
					agent->trailList_.push_back(RVO::Vector2(agent->navAgent_.corners[j], agent->navAgent_.corners[j+1]));
				}
				j=j+2;
			}
			return 1;
		}

		return NULL;
	}

	size_t RVOSimulator::addAgent(const Vector2 &position)
	{
		if (defaultAgent_ == NULL) {
			return RVO_ERROR;
		}

		Agent *agent = new Agent(this);

		agent->position_ = position;
		agent->maxNeighbors_ = defaultAgent_->maxNeighbors_;
		agent->maxSpeed_ = defaultAgent_->maxSpeed_;
		agent->neighborDist_ = defaultAgent_->neighborDist_;
		agent->radius_ = defaultAgent_->radius_;
		agent->timeHorizon_ = defaultAgent_->timeHorizon_;
		agent->timeHorizonObst_ = defaultAgent_->timeHorizonObst_;
		agent->velocity_ = defaultAgent_->velocity_;

		agentInit(&agent->navAgent_);
		float start[] = {position.x(), position.y()};
		vcpy(agent->navAgent_.pos, start);

		agent->id_ = ++count_;

		agents_.push_back(agent);

		return agent->id_;
	}

	Agent *RVOSimulator::addAgentEx(const Vector2 &position, float radius, float maxSpeed)
	{
		if (defaultAgent_ == NULL) {
			return NULL;
		}

		Agent *agent = new Agent(this);

		agent->position_ = position;
		agent->maxNeighbors_ = defaultAgent_->maxNeighbors_;
		agent->maxSpeed_ = maxSpeed;
		agent->neighborDist_ = radius + 10;
		agent->radius_ = radius;
		agent->timeHorizon_ = defaultAgent_->timeHorizon_;
		agent->timeHorizonObst_ = defaultAgent_->timeHorizonObst_;
		agent->velocity_ = defaultAgent_->velocity_;

		agentInit(&agent->navAgent_);
		float start[] = {position.x(), position.y()};
		vcpy(agent->navAgent_.pos, start);

		agent->id_ = ++count_;

		agents_.push_back(agent);

		return agent;
	}

	size_t RVOSimulator::addAgent(const Vector2 &position, float neighborDist, size_t maxNeighbors, float timeHorizon, float timeHorizonObst, float radius, float maxSpeed, const Vector2 &velocity)
	{
		Agent *agent = new Agent(this);

		agent->position_ = position;
		agent->maxNeighbors_ = maxNeighbors;
		agent->maxSpeed_ = maxSpeed;
		agent->neighborDist_ = neighborDist;
		agent->radius_ = radius;
		agent->timeHorizon_ = timeHorizon;
		agent->timeHorizonObst_ = timeHorizonObst;
		agent->velocity_ = velocity;

		agentInit(&agent->navAgent_);
		float start[] = {position.x(), position.y()};
		vcpy(agent->navAgent_.pos, start);

		agent->id_ = ++count_;

		agents_.push_back(agent);

		return agent->id_;
	}

	void RVOSimulator::delAgent(Agent *agent)
	{
		if (agent) {
			std::vector<Agent *>::iterator iter = agents_.begin();
			for (; iter!=agents_.end(); iter++) {
				if (*iter == agent) {
					agents_.erase(iter);break;
				}
			}

			iter = kdTree_->agents_.begin();
			for (; iter!=kdTree_->agents_.end(); iter++) {
				if (*iter == agent) {
					kdTree_->agents_.erase(iter);break;
				}
			}
			delete agent;
		}
	}

	Obstacle *RVOSimulator::addObstacle(const std::vector<Vector2> &vertices)
	{
		if (vertices.size() < 2) {
			return NULL;
		}

		Obstacle *firstObstacle = NULL;
		const size_t obstacleNo = obstacles_.size();

		for (size_t i = 0; i < vertices.size(); ++i) {
			Obstacle *obstacle = new Obstacle();
			obstacle->point_ = vertices[i];

			if (i != 0) {
				obstacle->prevObstacle_ = obstacles_.back();
				obstacle->prevObstacle_->nextObstacle_ = obstacle;
			}
			else {
				firstObstacle = obstacle;
			}

			if (i == vertices.size() - 1) {
				obstacle->nextObstacle_ = obstacles_[obstacleNo];
				obstacle->nextObstacle_->prevObstacle_ = obstacle;
			}

			obstacle->unitDir_ = normalize(vertices[(i == vertices.size() - 1 ? 0 : i + 1)] - vertices[i]);

			if (vertices.size() == 2) {
				obstacle->isConvex_ = true;
			}
			else {
				obstacle->isConvex_ = (leftOf(vertices[(i == 0 ? vertices.size() - 1 : i - 1)], vertices[i], vertices[(i == vertices.size() - 1 ? 0 : i + 1)]) >= 0.0f);
			}

			obstacle->id_ = obstcount_++;

			obstacles_.push_back(obstacle);
		}

		return firstObstacle;
	}
		
	void RVOSimulator::delObstacle(Obstacle *obstacle)
	{
		while(obstacle){
			std::vector<Obstacle *>::iterator iter = obstacles_.begin();
			for (; iter!=obstacles_.end(); iter++) {
				if (*iter == obstacle) {
					obstacles_.erase(iter);break;
				}
			}

			if (obstacle->prevObstacle_) {
				obstacle->prevObstacle_->nextObstacle_ = NULL;
			}
			
			if (obstacle->nextObstacle_) {
				obstacle = obstacle->nextObstacle_;
				delete obstacle->prevObstacle_;
				obstacle->prevObstacle_ = NULL;
			}
			else {
				delete obstacle;
				obstacle = NULL;
			}
		}
	}

	void RVOSimulator::doStep()
	{
		kdTree_->buildAgentTree();

		if (isBuildObstFlag_) {
			kdTree_->buildObstacleTree();
			isBuildObstFlag_ = false;
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < static_cast<int>(agents_.size()); ++i) {
			agents_[i]->computeNeighbors();
			agents_[i]->computeNewVelocity();
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < static_cast<int>(agents_.size()); ++i) {
			agents_[i]->update();
		}

		globalTime_ += timeStep_;
	}

	void RVOSimulator::doSingleStep(size_t agentNo)
	{
		if (isBuildObstFlag_) {
			kdTree_->buildObstacleTree();
			isBuildObstFlag_ = false;
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif

		agents_[agentNo]->computeNeighborsObst();
		agents_[agentNo]->computeNewVelocity();

#ifdef _OPENMP
#pragma omp parallel for
#endif
		agents_[agentNo]->update();
	}
	
	size_t RVOSimulator::getAgentAgentNeighbor(size_t agentNo, size_t neighborNo) const
	{
		return agents_[agentNo]->agentNeighbors_[neighborNo].second->id_;
	}

	size_t RVOSimulator::getAgentMaxNeighbors(size_t agentNo) const
	{
		return agents_[agentNo]->maxNeighbors_;
	}

	float RVOSimulator::getAgentMaxSpeed(size_t agentNo) const
	{
		return agents_[agentNo]->maxSpeed_;
	}

	float RVOSimulator::getAgentNeighborDist(size_t agentNo) const
	{
		return agents_[agentNo]->neighborDist_;
	}

	size_t RVOSimulator::getAgentNumAgentNeighbors(size_t agentNo) const
	{
		return agents_[agentNo]->agentNeighbors_.size();
	}

	size_t RVOSimulator::getAgentNumObstacleNeighbors(size_t agentNo) const
	{
		return agents_[agentNo]->obstacleNeighbors_.size();
	}

	size_t RVOSimulator::getAgentNumORCALines(size_t agentNo) const
	{
		return agents_[agentNo]->orcaLines_.size();
	}

	size_t RVOSimulator::getAgentObstacleNeighbor(size_t agentNo, size_t neighborNo) const
	{
		return agents_[agentNo]->obstacleNeighbors_[neighborNo].second->id_;
	}

	size_t RVOSimulator::getAgentId(size_t agentNo) const
	{
		return agents_[agentNo]->id_;
	}

	const Line &RVOSimulator::getAgentORCALine(size_t agentNo, size_t lineNo) const
	{
		return agents_[agentNo]->orcaLines_[lineNo];
	}

	const Vector2 &RVOSimulator::getAgentPosition(size_t agentNo) const
	{
		return agents_[agentNo]->position_;
	}

	const Vector2 &RVOSimulator::getAgentTrail(size_t agentNo) const
	{
		return agents_[agentNo]->trailList_[0];
	}

	const Vector2 &RVOSimulator::getAgentPrefVelocity(size_t agentNo) const
	{
		return agents_[agentNo]->prefVelocity_;
	}

	float RVOSimulator::getAgentRadius(size_t agentNo) const
	{
		return agents_[agentNo]->radius_;
	}

	float RVOSimulator::getAgentTimeHorizon(size_t agentNo) const
	{
		return agents_[agentNo]->timeHorizon_;
	}

	float RVOSimulator::getAgentTimeHorizonObst(size_t agentNo) const
	{
		return agents_[agentNo]->timeHorizonObst_;
	}

	const Vector2 &RVOSimulator::getAgentVelocity(size_t agentNo) const
	{
		return agents_[agentNo]->velocity_;
	}

	float RVOSimulator::getGlobalTime() const
	{
		return globalTime_;
	}

	size_t RVOSimulator::getNumAgents() const
	{
		return agents_.size();
	}

	size_t RVOSimulator::getNumTrails(size_t i) const
	{
		return agents_[i]->trailList_.size();
	}

	bool RVOSimulator::getTrailOverFlag(size_t i) const
	{
		return agents_[i]->isTrailOverFlag_;
	}

	void RVOSimulator::setTrailOverFlag(size_t i)
	{
		agents_[i]->isTrailOverFlag_ = false;
	}

	void RVOSimulator::setBuildObstFlag(bool flag)
	{
		isBuildObstFlag_ = flag;
	}

	size_t RVOSimulator::getNumObstacleVertices() const
	{
		return obstacles_.size();
	}

	const Vector2 &RVOSimulator::getObstacleVertex(size_t vertexNo) const
	{
		return obstacles_[vertexNo]->point_;
	}

	size_t RVOSimulator::getNextObstacleVertexNo(size_t vertexNo) const
	{
		return obstacles_[vertexNo]->nextObstacle_->id_;
	}

	size_t RVOSimulator::getPrevObstacleVertexNo(size_t vertexNo) const
	{
		return obstacles_[vertexNo]->prevObstacle_->id_;
	}

	float RVOSimulator::getTimeStep() const
	{
		return timeStep_;
	}

	void RVOSimulator::setObstacles(const std::vector<Obstacle *> obstacles)
	{
		size_t size = obstacles.size();
		obstacles_.resize(size);
		for (size_t i = 0; i < size; ++i) {
			obstacles_[i] = obstacles[i];
		}
		obstcount_ = size;
	}

	void RVOSimulator::processObstacles()
	{
		kdTree_->buildObstacleTree();
	}

	bool RVOSimulator::queryVisibility(const Vector2 &point1, const Vector2 &point2, float radius) const
	{
		return kdTree_->queryVisibility(point1, point2, radius);
	}

	void RVOSimulator::setAgentDefaults(float neighborDist, size_t maxNeighbors, float timeHorizon, float timeHorizonObst, float radius, float maxSpeed, const Vector2 &velocity)
	{
		if (defaultAgent_ == NULL) {
			defaultAgent_ = new Agent(this);
		}

		defaultAgent_->maxNeighbors_ = maxNeighbors;
		defaultAgent_->maxSpeed_ = maxSpeed;
		defaultAgent_->neighborDist_ = neighborDist;
		defaultAgent_->radius_ = radius;
		defaultAgent_->timeHorizon_ = timeHorizon;
		defaultAgent_->timeHorizonObst_ = timeHorizonObst;
		defaultAgent_->velocity_ = velocity;
	}
	
	struct Navmesh* RVOSimulator::createNavmesh(std::vector<int> const& trisList, std::vector<std::pair<float, float> > const& pointList)
	{
		nav_ = navmeshCreateEx(trisList, pointList);
		if (nav_) {
			return nav_;
		}
		else {
			return NULL;
		}
	}

	void RVOSimulator::setAgentMaxNeighbors(size_t agentNo, size_t maxNeighbors)
	{
		agents_[agentNo]->maxNeighbors_ = maxNeighbors;
	}

	void RVOSimulator::setAgentMaxSpeed(size_t agentNo, float maxSpeed)
	{
		agents_[agentNo]->maxSpeed_ = maxSpeed;
	}

	void RVOSimulator::setAgentNeighborDist(size_t agentNo, float neighborDist)
	{
		agents_[agentNo]->neighborDist_ = neighborDist;
	}

	void RVOSimulator::setAgentPosition(size_t agentNo, const Vector2 &position)
	{
		agents_[agentNo]->position_ = position;

		float start[] = {position.x(), position.y()};
		vcpy(agents_[agentNo]->navAgent_.pos, start);
		agents_[agentNo]->stop();
	}

	void RVOSimulator::setAgentPosition(Agent* agent, const Vector2 &position)
	{
		agent->position_ = position;

		float start[] = {position.x(), position.y()};
		vcpy(agent->navAgent_.pos, start);
		agent->stop();
	}

	void RVOSimulator::setAgentPrefVelocity(size_t agentNo, const Vector2 &prefVelocity)
	{
		agents_[agentNo]->prefVelocity_ = prefVelocity;
	}

	void RVOSimulator::setAgentRadius(size_t agentNo, float radius)
	{
		agents_[agentNo]->radius_ = radius;
	}

	void RVOSimulator::setAgentTimeHorizon(size_t agentNo, float timeHorizon)
	{
		agents_[agentNo]->timeHorizon_ = timeHorizon;
	}

	void RVOSimulator::setAgentTimeHorizonObst(size_t agentNo, float timeHorizonObst)
	{
		agents_[agentNo]->timeHorizonObst_ = timeHorizonObst;
	}
	
	void RVOSimulator::setAgentVelocity(size_t agentNo, const Vector2 &velocity)
	{
		agents_[agentNo]->velocity_ = velocity;
	}
	
	void RVOSimulator::getAgentAABB(size_t agentNo, b2AABB& aabb)
	{
		agents_[agentNo]->getAABB(aabb);
	}

	void RVOSimulator::getAgentPolygoShape(size_t agentNo, b2PolygonShape& poly)
	{
		agents_[agentNo]->getPolygoShape(poly);
	}

	void RVOSimulator::setTimeStep(float timeStep)
	{
		timeStep_ = timeStep;
	}

	int RVOSimulator::getAgentType(size_t agentNo)
	{
		return agents_[agentNo]->getAgentType();
	}

	void RVOSimulator::addPortNo(int port_no)
	{
		portNoMap_[port_no] = 1;
	}

	void RVOSimulator::delPortNo(int port_no)
	{
		portNoMap_.erase(port_no);
	}

	Agent* RVOSimulator::getAgentById(size_t id)
	{
		std::vector<Agent *>::iterator iter = agents_.begin();
		for (; iter!=agents_.end(); iter++) {
			if ((*iter)->id_ == id) {
				return (*iter);
			}
		}

		return NULL;
	}

	bool RVOSimulator::createNavigation(InputGeom* geom, BuildContext* ctx)
	{
		sample_ = new Sample_SoloMesh();
		if (sample_ && geom && ctx) 
		{
			sample_->setContext(ctx);
			sample_->handleMeshChanged(geom);
		}else {
			return false;
		}

		if (!sample_->handleBuild())
			return false;

		return true;
	}
}
