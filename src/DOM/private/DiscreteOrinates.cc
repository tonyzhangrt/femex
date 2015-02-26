/*
 * DiscreteOrinates.cpp
 *
 *  Created on: Feb 21, 2015
 *      Author: lurker
 */

#include "DiscreteOrinates.h"

namespace Core {

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle) noexcept {
	nAngle = static_cast<int32_t>(*mxGetPr(_nAngle));
	initAngle = 0.;
	Ray.resize(nAngle);

}

DiscreteOrinates::DiscreteOrinates(MatlabPtr _nAngle, MatlabPtr _initAngle) noexcept {
	nAngle = static_cast<int32_t>(*mxGetPr(_nAngle));
	initAngle = static_cast<int32_t>(*mxGetPr(_initAngle));
	Ray.resize(nAngle);
}

DiscreteOrinates::~DiscreteOrinates() {

	int32_t tmp_i, tmp_j;
	if (nAngle != 0) {
		for (int32_t i = 0 ; i < nAngle; i++){
			tmp_i = Ray[i].size();
			for (int32_t j = 0; j < tmp_i; j++){
				Ray[i][j].clear();
			}
			Ray[i].clear();
		}
		Ray.clear();
	}
#ifdef DEBUG
	mexPrintf("Discrete Orinates Method detached.\n");
#endif
}

void DiscreteOrinates::RayInt(MatlabPtr nodes, MatlabPtr elems,
		MatlabPtr neighbors){


	auto numberofnodes        = mxGetN(nodes);
	auto numberofelems        = mxGetN(elems);
	auto numberofnodesperelem = mxGetM(elems);
//	auto numberofedges        = mxGetN(edges);
//	auto numberofnodesperedge = mxGetM(edges);


	auto pnodes               = mxGetPr(nodes);
	auto pelems               = (int32_t*)mxGetPr(elems);
	auto pneighbors           = (int32_t*)mxGetPr(neighbors);
//	auto pedges               = (int32_t*)mxGetPr(edges);

	mwSize vertex, e_vl, e_vr;;
	Real_t theta , x1, x2, y1, y2,  a, b;
	std::unordered_set<int32_t> visited;

	Real_t ret, t, eta;
	Real_t q_x1, q_y1, q_x2, q_y2, q_x3, q_y3;
	Real_t q_eta, q_t;
	int32_t index;

	std::vector<Real_t> tmp;

	bool intersect;
	std::queue<int32_t> q_elem;
	std::unordered_set<int32_t> s_elem;

	Raylet tmp_ray;

	/*
	 * resize Ray to hold raylets
	 */
	for (int32_t i = 0; i < nAngle; i++) {
		Ray[i].resize(numberofnodes);
	}


	/*
	 * main part
	 */
	for (int32_t i = 0; i < nAngle; i++) {
		std::cout << "the " << i << "th run" << std::endl;
		theta = initAngle + 2 * i * M_PIl / nAngle;
		visited.clear();
		for (int32_t j = 0; j < numberofelems; j++) {
			for (int32_t k = 0; k < numberofnodesperelem; k++){
				/*
				 * index of vertex in nodes.
				 */
				vertex = pelems[numberofnodesperelem * j + k] - 1;
				if (visited.find(vertex) == visited.end()){
					visited.insert(vertex);
					/*
					 * calculate the ray intersects with the boundary.
					 */
					a = pnodes[2 * vertex    ];
					b = pnodes[2 * vertex + 1];

					// clean up set.
					s_elem.clear();
					// push current element
					q_elem.push(j);
					s_elem.insert(j);
					/*
					 * check the jth.
					 */
					index = j;

					q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
					q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

					q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
					q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

					q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
					q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];

					RayTrace(tmp, intersect, q_t, q_eta,
							q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
							a, b, theta);

					RayTrim(tmp, a, b);

					if (tmp.size()) {
							tmp_ray.elem = index;
						if (tmp.size() == 2) {
							tmp_ray.first[0] = a;
							tmp_ray.first[1] = b;
							tmp_ray.second[0] = tmp[0];
							tmp_ray.second[1] = tmp[1];
						}
						else {
							tmp_ray.first[0] = tmp[0];
							tmp_ray.first[1] = tmp[1];
							tmp_ray.second[0] = tmp[2];
							tmp_ray.second[1] = tmp[3];
						}

						/*
						 * remove duplicates, since sorted, only test first two.
						 */

						if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
							fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){

							Ray[i][vertex].push_back(tmp_ray);
						}
					}
#ifdef DEBUG
					if (tmp.size()){
						if (tmp.size() == 2){
							std::cout << index << " ["<<a<< "," << b << "] -> ["<< tmp[0] << "," << tmp[1] <<"]" << std::endl;
						}
						else {
							std::cout << index << " ["<<tmp[0]<< "," << tmp[1] << "] -> ["<< tmp[2] << "," << tmp[3] <<"]" << std::endl;
						}

					}
#endif
					while (!q_elem.empty()){
						auto top = q_elem.front();
						q_elem.pop();
						// loop over neighbors.
						// fail safe strategy, if it is possible to intersect, must push.
						for (int32_t q_elem_i = 0; q_elem_i < 3; q_elem_i ++){

							index = pneighbors[3 * top + q_elem_i] - 1;

							// test if it is going to push them in.
							// first : valid index of element
							if (index > -1 && s_elem.find(index) == s_elem.end()) {

								/*
								 * calculate the open angle limit. And find the suitable one.
								 *
								 * Target function is max(min(abs(langle - angle), abs(rangle - angle)))
								 *
								 * if both gives 0, then push both into queue.
								 */

								q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
								q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

								q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
								q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

								q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
								q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];


								/*
								 * [a , b] 's angle for the element,
								 */

								RayTrace(tmp, intersect, q_t, q_eta,
										q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
										a, b, theta);

								if (intersect){
									q_elem.push(index);
									s_elem.insert(index);
									/*
									 * calculate integral over this element
									 *
									 * first step: calculate the intersection,
									 * second step: integral
									 */
									RayTrim(tmp, a, b);

									if (tmp.size()) {
											tmp_ray.elem = index;
										if (tmp.size() == 2) {
											tmp_ray.first[0] = a;
											tmp_ray.first[1] = b;
											tmp_ray.second[0] = tmp[0];
											tmp_ray.second[1] = tmp[1];
										}
										else {
											tmp_ray.first[0] = tmp[0];
											tmp_ray.first[1] = tmp[1];
											tmp_ray.second[0] = tmp[2];
											tmp_ray.second[1] = tmp[3];
										}

										/*
										 * remove duplicates, since sorted, only test first two.
										 */
										if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
											fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){

											Ray[i][vertex].push_back(tmp_ray);
										}
									}
#ifdef DEBUG
									if (tmp.size()){
										if (tmp.size() == 2){
											std::cout << index << " ["<<a<< "," << b << "] -> ["<< tmp[0] << "," << tmp[1] <<"]" << std::endl;
										}
										else {
											std::cout << index << " ["<<tmp[0]<< "," << tmp[1] << "] -> ["<< tmp[2] << "," << tmp[3] <<"]" << std::endl;
										}

									}
#endif
								}
							}
						}
					}
				}
				Ray[i][vertex].shrink_to_fit();
			}
		}
	}
//	RayShow();
}

void DiscreteOrinates::RayTrace(std::vector<Real_t>& tmp, bool& intersect, Real_t& q_t, Real_t& q_eta,
			Real_t& q_x1, Real_t& q_y1, Real_t& q_x2,
			Real_t& q_y2, Real_t& q_x3, Real_t& q_y3,
			Real_t& a, Real_t& b, Real_t& theta){

	intersect = false;
	tmp.clear();
	// 1 - 2 intersects.
	if (fabs(INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x1 - q_x2) + MEX_EPS < fabs(q_x1 - a) + fabs(q_x2 - a) ||
					fabs(q_y1 - q_y2) + MEX_EPS < fabs(q_y1 - b) + fabs(q_y2 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);

		q_eta = INTERSECT_DET(a, b, q_x2, q_y2, theta)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x1 + (1- q_eta) * q_x2);
			tmp.push_back(q_eta * q_y1 + (1 -q_eta) * q_y2);
		}
	}
	// 2 - 3
	if (fabs(INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x2 - q_x3) + MEX_EPS < fabs(q_x2 - a) + fabs(q_x3 - a) ||
					fabs(q_y2 - q_y3) + MEX_EPS < fabs(q_y2 - b) + fabs(q_y3 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		q_eta = INTERSECT_DET(a, b, q_x3, q_y3, theta)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x2 + (1 - q_eta) * q_x3);
			tmp.push_back(q_eta * q_y2 + (1 - q_eta) * q_y3);
		}
	}
	// 3 - 1
	if (fabs(INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta)) < SCALE * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)) < SCALE * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x3 - q_x1) + MEX_EPS < fabs(q_x3 - a) + fabs(q_x1 - a) ||
					fabs(q_y3 - q_y1) + MEX_EPS < fabs(q_y3 - b) + fabs(q_y1 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);

		q_eta = INTERSECT_DET(a, b, q_x1, q_y1, theta)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x3 + (1 - q_eta) * q_x1);
			tmp.push_back(q_eta * q_y3 + (1 - q_eta) * q_y1);
		}
	}
}

void DiscreteOrinates::RayTrim(std::vector<Real_t>& tmp, Real_t &a, Real_t &b){
	for (int32_t tmp_i = 0; tmp_i < tmp.size()/2; tmp_i ++){

		if (fabs(a - tmp[2 * tmp_i]) + fabs(b - tmp[2 * tmp_i + 1]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2 * tmp_i, tmp.begin() + 2 * tmp_i + 2);
			tmp_i --;
		}
	}

	if (tmp.size() == 6) {

	/*
	 * remove duplicate coordinates
	 */
	// 2 - 3 duplicates
		if (fabs(tmp[2] - tmp[4]) + fabs(tmp[3] - tmp[5]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2, tmp.begin() + 4);
		}
	// 3 - 1 duplicates or 1 - 2 duplicates
		else {
			tmp.erase(tmp.begin(), tmp.begin() + 2);
		}
	}

	if (tmp.size() == 4) {
	/*
	 * remove duplicates
	 */
		if (fabs(tmp[0] - tmp[2]) + fabs(tmp[1] - tmp[3]) < MEX_EPS) {
			tmp.clear();
		}
		else {
			// reorder
			if (pow(tmp[0] - a, 2) + pow(tmp[1] - b, 2) > pow(tmp[2] - a, 2) + pow(tmp[3] - b, 2)){
				swap(tmp[0], tmp[2]);
				swap(tmp[1], tmp[3]);
			}
		}
	}
}

/*
 * depreciated later.
 */
void DiscreteOrinates::RayShow(){

	int32_t tmp_i, tmp_j;
//	size_t tmp_total = 0;
	if (nAngle != 0) {
		for (int32_t i = 0 ; i < nAngle; i++){
			tmp_i = Ray[i].size();
			for (int32_t j = 0; j < tmp_i; j++){
				tmp_j = Ray[i][j].size();
//				tmp_total += tmp_j * 36;

				for (int32_t k = 0; k < tmp_j; k++) {
					std::cout << i << "th Angle, "
							<< j << "th node, "
							<< k << "th raylet: passes through "
							<< Ray[i][j][k].elem << ", starting from "
							<< Ray[i][j][k].first << " --> "
							<< Ray[i][j][k].second << std::endl;
				}
			}
		}
	}
//	std::cout << tmp_total << std::endl;
}



/*
 * Source Iteration.
 *
 */

void DiscreteOrinates::SourceIteration_init(MatlabPtr Fcn,
		MatlabPtr Sigma_t_Fcn, MatlabPtr Sigma_s_Fcn,
		MatlabPtr nodes, MatlabPtr elems){
//	Source.resize(nAngle);

	/*
	 * consider source as f(x).
	 */
	auto numberofnodes = mxGetN(nodes);
	auto interp        = mxGetPr(Fcn);
	auto interp_t      = mxGetPr(Sigma_t_Fcn);
	auto interp_s      = mxGetPr(Sigma_s_Fcn);

	auto pnodes        = mxGetPr(nodes);
	auto pelems        = (int32_t*)mxGetPr(elems);
	auto numberofnodesperelem = mxGetM(elems);

	mxAssert(mxGetN(Fcn) == numberofnodes,
			"DiscreteOrinates::SourceIteration_init::Dimension does not match.\n");

	Source.resize(numberofnodes);
	Sigma_t.resize(numberofnodes);
	Sigma_s.resize(numberofnodes);

	memcpy(&Source[0], interp, sizeof(Real_t) * numberofnodes);
	memcpy(&Sigma_t[0], interp_t, sizeof(Real_t) * numberofnodes);
	memcpy(&Sigma_s[0], interp_s, sizeof(Real_t) * numberofnodes);

	RHS.resize(numberofnodes);
	Average.resize(numberofnodes);

	Output.resize(nAngle);


	for (int32_t s_i; s_i < nAngle; s_i ++) {
		Output[s_i].resize(numberofnodes);
	}
}



void DiscreteOrinates::SourceIteration_iter(MatlabPtr nodes, MatlabPtr elems){
	// kernel as 1/2/pi constant

	auto pnodes        = mxGetPr(nodes);
	auto pelems        = (int32_t *)mxGetPr(elems);
	auto numberofnodesperelem = mxGetM(elems);

	auto numberofnodes = Source.size();

	for (int32_t s_j = 0; s_j < numberofnodes; s_j++){
		Average[s_j] = 0.;
		for (int32_t s_i = 0; s_i < nAngle; s_i++) {
			Average[s_j] += Output[s_i][s_j];
		}
		RHS[s_j] = Sigma_s[s_j] * Average[s_j]/nAngle;
		RHS[s_j] += Source[s_j];
	}

	mwSize vertex_1, vertex_2, vertex_3;
	Real_t x1, y1, x2, y2, x3, y3, det, lambda, eta, length, accum_s, accum_v;

	Real_t lv, rv, ls, rs;

	for (int32_t s_i = 0; s_i < nAngle; s_i++) {

		for (int32_t s_j = 0; s_j < numberofnodes; s_j++){

			accum_s = 0.;
			accum_v = 0.;

			if (Ray[s_i][s_j].size()){

				for (auto it : Ray[s_i][s_j]){
					vertex_1 = pelems[it.elem * numberofnodesperelem ] - 1;
					vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
					vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

					x1 = pnodes[2 * vertex_1];
					y1 = pnodes[2 * vertex_1 + 1];
					x2 = pnodes[2 * vertex_2];
					y2 = pnodes[2 * vertex_2 + 1];
					x3 = pnodes[2 * vertex_3];
					y3 = pnodes[2 * vertex_3 + 1];

					/*
					 * first node
					 */
					det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);

					eta = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
					eta /= det;

					lambda = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
					lambda /= det;

					lv = lambda * RHS[vertex_1] + eta * RHS[vertex_2] +
							(1 - lambda - eta) * RHS[vertex_3];
					ls = lambda * Sigma_t[vertex_1] + eta * Sigma_t[vertex_2] +
							(1 - lambda - eta) * Sigma_t[vertex_3];

					/*
					 * second node
					 */
					eta = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
					eta /= det;

					lambda = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
					lambda /= det;

					rv = lambda * RHS[vertex_1] + eta * RHS[vertex_2] +
							(1 - lambda - eta) * RHS[vertex_3];
					rs = lambda * Sigma_t[vertex_1] + eta * Sigma_t[vertex_2] +
							(1 - lambda - eta) * Sigma_t[vertex_3];

					/*
					 * length
					 */
					length = sqrt(pow(it.first[0] - it.second[0], 2) + pow(it.first[1] - it.second[1], 2));

					accum_v += exp(-accum_s) * lv * length/2.0;

					accum_s += (rs + ls) * length/ 2.0;

					accum_v += exp(-accum_s) * rv * length/2.0;
				}
			}
			else{
				accum_v = 0.;
			}
			Output[s_i][s_j] = accum_v;
		}
	}
}


void DiscreteOrinates::SourceIteration_accl(){

}


void DiscreteOrinates::SourceIteration_output(){
	// move to mexplus
}



} /* namespace Core */


using namespace Core;
template class mexplus::Session<DiscreteOrinates>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<DiscreteOrinates>::create(new DiscreteOrinates(const_cast<mxArray*>(input.get(0)))));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<DiscreteOrinates>::destroy(input.get(0));
}

MEX_DEFINE(rayint) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->RayInt(CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]));
}

MEX_DEFINE(rayshow) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->RayShow();
}

MEX_DEFINE(si_init) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_init(CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]),CAST(prhs[4]), CAST(prhs[5]));
}

MEX_DEFINE(si_iter) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 0);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	DOM->SourceIteration_iter(CAST(prhs[1]), CAST(prhs[2]));
}

MEX_DEFINE(si_output) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 1);

	DiscreteOrinates* DOM = Session<DiscreteOrinates>::get(input.get(0));

	plhs[0] = mxCreateNumericMatrix(1, DOM->Average.size(), mxDOUBLE_CLASS, mxREAL);
	memcpy(mxGetPr(plhs[0]), &(DOM->Average[0]), DOM->Average.size()*sizeof(Real_t));
}

}

MEX_DISPATCH
