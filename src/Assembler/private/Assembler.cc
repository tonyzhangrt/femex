/*
 * Assembler.cc
 *
 *  Created on: Oct 9, 2014
 *      Author: lurker
 */
#include "Assembler.h"

Assembler::Assembler(){


}

Assembler::~Assembler() {
#ifdef DEBUG
	mexPrintf("Assembler detached\n");
#endif
}


/*
 * Extract information of basis on Reference Triangle
 */
void Assembler::Reference(MatlabPtr &F, MatlabPtr &DX, MatlabPtr &DY,
		MatlabPtr Points, MatlabPtr QPoints){

	auto _numberofpoints = mxGetN(Points);
	auto _numberofqpoints = mxGetN(QPoints);


	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	auto Vander = mxCreateNumericMatrix(_numberofpoints,_numberofpoints,mxDOUBLE_CLASS, mxREAL);
	auto VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderY = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	int deg = round((sqrt(8*_numberofpoints + 1) - 3)/2);

	if ((deg + 1) * (deg + 2)/2 != _numberofpoints){
		mexErrMsgTxt("Invalid length of input nodes\n");
	}

	// extract all nodes from promoted nodes.
	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderF_ptr = mxGetPr(VanderF);
	auto VanderX_ptr = mxGetPr(VanderX);
	auto VanderY_ptr = mxGetPr(VanderY);

	auto nodes_ptr   = mxGetPr(Points);
	auto qnodes_ptr  = mxGetPr(QPoints);

	// basis on all nodes

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*Vander_ptr++ = pow(nodes_ptr[2*col], i - j)*pow(nodes_ptr[2*col + 1], j);
			}
		}
	}



	// basis on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++){
				*VanderF_ptr++ = pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j);
			}
		}
	}

	// partial x on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == i){
					*VanderX_ptr++ = 0.;
				}
				else{
					*VanderX_ptr++ = (i - j) * pow(qnodes_ptr[2*col], i - j - 1)*pow(qnodes_ptr[2*col + 1], j);
				}
			}
		}
	}

	// partial y on qnodes
	for (size_t col = 0; col < _numberofqpoints; col++){
		for (size_t i = 0; i < deg + 1; i++){
			for (size_t j = 0; j < i + 1; j++) {
				if (j == 0){
					*VanderY_ptr++ = 0.;
				}
				else{
					*VanderY_ptr++ = j* pow(qnodes_ptr[2*col], i - j)*pow(qnodes_ptr[2*col + 1], j - 1);
				}
			}
		}
	}



	mxArray* RHS_f[] = {Vander, VanderF};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	mxArray* RHS_y[] = {Vander, VanderY};
	mexCallMATLAB(1, &DY, 2, RHS_y, "mldivide");


	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);
	mxDestroyArray(VanderY);


}

/*
 * 1D reference
 */
void Assembler::Reference(MatlabPtr& F, MatlabPtr& DX, MatlabPtr Degree, MatlabPtr QPoints){
	// Reference segment [-1 , 1],
	// with (Degree - 1) equal-spaced points in between.

	auto _numberofpoints = static_cast<int32_t>(*Matlab_Cast<Real_t>(Degree) + 1);

	vector<Real_t> Points(_numberofpoints);
	Points[0] = -1.;
	Points[1] = 1. ;
	for (size_t i = 2; i < _numberofpoints; i++) {
		Points[i] = (-1.) + 2.0*(i - 1)/static_cast<Real_t>(_numberofpoints - 1);
	}

	// 1D vector-row-majored
	auto _numberofqpoints = mxGetN(QPoints);
	auto qnodes_ptr       = mxGetPr(QPoints);

	/*
	 *  Temporary Arrays, will be destroyed later.
	 */
	auto Vander = mxCreateNumericMatrix(_numberofpoints, _numberofpoints ,mxDOUBLE_CLASS, mxREAL);
	auto VanderF = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	auto VanderX = mxCreateNumericMatrix(_numberofpoints, _numberofqpoints, mxDOUBLE_CLASS, mxREAL);


	// extract all nodes from promoted nodes.
	auto Vander_ptr  = mxGetPr(Vander);
	auto VanderF_ptr = mxGetPr(VanderF);
	auto VanderX_ptr = mxGetPr(VanderX);

	for (size_t col = 0; col < _numberofpoints; col++){
		for (size_t i = 0; i < _numberofpoints; i++){
			*Vander_ptr++ = pow(Points[col], i);
		}
	}

	for (size_t col = 0; col < _numberofqpoints; col++) {
		for (size_t i = 0; i < _numberofpoints; i++) {
			*VanderF_ptr++ = pow(qnodes_ptr[col], i);
		}
	}

	for (size_t col = 0; col < _numberofqpoints; col++) {
		for (size_t i = 0 ; i < _numberofpoints ; i++) {
			if (i == 0) {
				*VanderX_ptr++ = 0.;
			}
			else{
				*VanderX_ptr++ = i * pow(qnodes_ptr[col], i - 1);
			}
		}
	}

	mxArray* RHS_f[] = {Vander, VanderF};
	mexCallMATLAB(1, &F, 2, RHS_f, "mldivide");

	mxArray* RHS_x[] = {Vander, VanderX};
	mexCallMATLAB(1, &DX, 2, RHS_x, "mldivide");

	Points.clear();
	mxDestroyArray(Vander);
	mxDestroyArray(VanderF);
	mxDestroyArray(VanderX);

}

void Assembler::AssembleBC(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes,MatlabPtr eRobin,
		MatlabPtr Ref, MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(eRobin);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofedges          = mxGetN(eRobin);
	auto numberofnodesperedge   = mxGetM(eRobin);
	auto numberofqnodes         = mxGetN(Ref);

	mwSize vertex_1, vertex_2;
	Real_t length;

	/*
	 * More codes but there is only one judgment at beginning.
	 * Performance does not drop.
	 */
	if (mxGetNumberOfElements(Fcn) == numberofedges*numberofqnodes){

		for (size_t i = 0; i < numberofedges; i++) {


			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;


			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				for (size_t k = 0; k < numberofnodesperedge; k++) {
					*pI++ = pedge_ptr[i*numberofnodesperedge + j];
					*pJ++ = pedge_ptr[i*numberofnodesperedge + k];
					*pV = 0.;

					for (size_t l = 0; l < numberofqnodes; l++) {
						*pV = *pV + Interp[i*numberofqnodes + l] *
								reference[j + l*numberofnodesperedge] *
								reference[k + l*numberofnodesperedge] *
								weights[l];
					}
					*pV *= length/2.0; pV++;
				}
			}
		}
	}
	else {
		for (size_t i = 0; i < numberofedges; i++) {


			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;


			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				for (size_t k = 0; k < numberofnodesperedge; k++) {
					*pI++ = pedge_ptr[i*numberofnodesperedge + j];
					*pJ++ = pedge_ptr[i*numberofnodesperedge + k];
					*pV = 0.;

					for (size_t l = 0; l < numberofqnodes; l++) {
						*pV = *pV +
								reference[j + l*numberofnodesperedge] *
								reference[k + l*numberofnodesperedge] *
								weights[l];
					}
					*pV *= *(Interp)*length/2.0; pV++;
				}
			}
		}
	}
}


void Assembler::AssembleBC(Real_t*& pNeumann, MatlabPtr Nodes,
		MatlabPtr QNodes, MatlabPtr eNeumann,
		MatlabPtr Ref, MatlabPtr Weights,  MatlabPtr Fcn) {
	// Fcn needs to be same size as eNeumann * qnodes

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(eNeumann);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofedges          = mxGetN(eNeumann);
	auto numberofnodesperedge   = mxGetM(eNeumann);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2;
	Real_t length, tmp;

	if (mxGetNumberOfElements(Fcn)  == numberofedges*numberofqnodes) {
		/*
		 * Fcn is a matrix.
		 */
		for (size_t i = 0; i < numberofedges; i++) {
			// integral over each interval
			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;

			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				// each basis
				tmp = 0.;
				for (size_t l = 0 ; l < numberofqnodes; l++) {
					tmp += Interp[i*numberofqnodes + l] * reference[j + l*numberofnodesperedge] * weights[l];
				}
				pNeumann[pedge_ptr[i*(numberofnodesperedge) + j] - 1] += tmp * length/2.;
			}
		}
	}
	else {
		/*
		 * Fcn is a number
		 */
		for (size_t i = 0; i < numberofedges; i++) {
			// integral over each interval
			vertex_1 = pedge_ptr[i*(numberofnodesperedge) ] - 1;
			vertex_2 = pedge_ptr[i*(numberofnodesperedge) + 1] - 1;

			length = sqrt(
					pow(pnodes_ptr[2*vertex_1] - pnodes_ptr[2*vertex_2], 2)
					+
					pow(pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1] , 2));

			for (size_t j = 0 ; j < numberofnodesperedge; j++) {
				// each basis
				tmp = 0.;
				for (size_t l = 0 ; l < numberofqnodes; l++) {
					tmp += reference[j + l*numberofnodesperedge] * weights[l];
				}
				pNeumann[pedge_ptr[i*(numberofnodesperedge) + j] - 1] += *(Interp)*tmp * length/2.;
			}
		}
	}
}

void Assembler::AssembleMass(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr Fcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;


	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l]*
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = (*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
	else if (mxGetNumberOfElements(Fcn) == numberofelem) {
		for (size_t i = 0; i < numberofelem; i++) {
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			auto r = Interp[i];
			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = r*(*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = *(Interp)*(*pV)*area;

					pI++; pJ++; pV++;
					if (k != j) {
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}

				}
			}
		}
	}
}

void Assembler::Qnodes2D(Real_t*& Coords, MatlabPtr Nodes, MatlabPtr QNodes, MatlabPtr Elems){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(QNodes);


	mwSize vertex_1, vertex_2 , vertex_3;

	for (size_t i = 0; i < numberofelem; i++) {

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		for (size_t l = 0; l < numberofqnodes; l++) {
			Coords[2*(i*numberofqnodes + l)] =
					pnodes_ptr[2 * vertex_1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3]*qnodes_ptr[2*l + 1];
			Coords[2*(i*numberofqnodes + l) + 1] =
					pnodes_ptr[2 * vertex_1 + 1] *(1 - qnodes_ptr[2*l] - qnodes_ptr[2*l + 1]) +
					pnodes_ptr[2 * vertex_2 + 1]*qnodes_ptr[2*l] +
					pnodes_ptr[2 * vertex_3 + 1]*qnodes_ptr[2*l + 1];
		}
	}
}


void Assembler::Qnodes1D(Real_t*& Coords, MatlabPtr Nodes, MatlabPtr QNodes, MatlabPtr Edges){
	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pedge_ptr            = (int32_t*)mxGetPr(Edges);

	auto numberofedges          = mxGetN(Edges);
	auto numberofnodesperedge   = mxGetM(Edges);
	auto numberofqnodes         = mxGetN(QNodes);

	mwSize vertex_1, vertex_2;

	if (mxGetM(Nodes) == 1) {
		// 1D problem
		for (size_t i = 0; i < numberofedges; i++) {
			vertex_1 = pedge_ptr[numberofnodesperedge*i] - 1;
			vertex_2 = pedge_ptr[numberofnodesperedge*i + 1] - 1;

			for (size_t l = 0; l < numberofqnodes; l++) {
				Coords[(i*numberofqnodes + l)] =
						pnodes_ptr[vertex_1] + pnodes_ptr[vertex_2] +
						(pnodes_ptr[vertex_2] - pnodes_ptr[vertex_1])*qnodes_ptr[l];

			}
		}
	}
	else {
		// 2D problem
		for (size_t i = 0; i < numberofedges; i++) {

			vertex_1 = pedge_ptr[numberofnodesperedge*i] - 1;
			vertex_2 = pedge_ptr[numberofnodesperedge*i + 1] - 1;

			for (size_t l = 0; l < numberofqnodes; l++) {
				Coords[2*(i*numberofqnodes + l)] =
						(pnodes_ptr[2*vertex_1] + pnodes_ptr[2*vertex_2] +
						(pnodes_ptr[2*vertex_2] - pnodes_ptr[2*vertex_1])*qnodes_ptr[l])/2.0;
				Coords[2*(i*numberofqnodes + l) + 1] =
						(pnodes_ptr[2*vertex_1 + 1] + pnodes_ptr[2*vertex_2 + 1] +
						(pnodes_ptr[2*vertex_2 + 1] - pnodes_ptr[2*vertex_1 + 1])*qnodes_ptr[l])/2.0;
			}
		}
	}

}



// calculate integral on boundary
void Assembler::AssembleLoad(Real_t*& pLoad, MatlabPtr Nodes,
		MatlabPtr QNodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr Fcn) {

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area, tmp;


	// Fcn cannot be a function handle, too slow

	/* @Revised: Fcn can be a function handle for one pass.
	 * Which means extra space to store all the values. However,
	 * we did not do it because Matlab can handle this easily.
	 */
	auto Fcn_ptr = Matlab_Cast<Real_t>(Fcn);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		if(mxGetNumberOfElements(Fcn) == numberofqnodes * numberofelem) {
			// Fcn has numberofqnodes * numberofelem elements
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Fcn_ptr[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp*area;
			}
		}
		else if (mxGetNumberOfElements(Fcn) == 1) {
			// Fcn is constant
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp += reference[j+ l*numberofnodesperelem]*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += *(Fcn_ptr)*tmp*area;
			}
		}
		else {
			// Other cases does not match
			mexErrMsgTxt("Error:Assembler:AssembleLoad::Failed with unexpected Fcn.\n");
		}
	}
}


void Assembler::AssembleStiff(Real_t* &pI, Real_t* &pJ, Real_t*&pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn)  == numberofelem*numberofqnodes) {
		// Fcn is a matrix
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else if (mxGetNumberOfElements(Fcn)  == numberofelem){
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			auto r = Interp[i];
			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = r*(*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
					}
					*pV = *(Interp)*(*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
}


/*
 * Stiff Kernel as two functions
 */

void Assembler::AssembleStiff(Real_t* &pI, Real_t* &pJ, Real_t*&pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn_X, MatlabPtr Fcn_Y) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp_X             = mxGetPr(Fcn_X);
	auto  Interp_Y             = mxGetPr(Fcn_Y);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn_X)  == numberofelem * numberofqnodes &&
			mxGetNumberOfElements(Fcn_Y) == numberofelem * numberofqnodes) {
		// Fcn is a matrix
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + (
								Interp_X[i*numberofqnodes + l] *
								(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								)
								+
								Interp_Y[i*numberofqnodes + l] *
								(
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
	else {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								(
								*(Interp_X)*
								(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								)
								+
								*(Interp_Y)*
								(
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)
								)*weights[l];
					}
					*pV = (*pV)/4.0/area;
					pI++; pJ++; pV++;
					if (j != k){
						*pI = *(pJ - 1);
						*pJ = *(pI - 1);
						*pV = *(pV - 1);
						pI++; pJ++; pV++;
					}
				}
			}
		}
	}
}


void Assembler::AssembleGradLoad(Real_t* &pLoad, MatlabPtr Nodes, MatlabPtr Elems, 
		MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn_X, MatlabPtr Fcn_Y){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp_X             = mxGetPr(Fcn_X);
	auto  Interp_Y             = mxGetPr(Fcn_Y);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area, tmp;
	Real_t Jacobian[2][2];


	// Fcn cannot be a function handle, too slow

	/* @Revised: Fcn can be a function handle for one pass.
	 * Which means extra space to store all the values. However,
	 * we did not do it because Matlab can handle this easily.
	 */
	auto Fcn_X_ptr = Matlab_Cast<Real_t>(Fcn_X);
	auto Fcn_Y_ptr = Matlab_Cast<Real_t>(Fcn_Y);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
		Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
		Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
		Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

		det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

		area = 0.5*fabs(det);		

		if (mxGetNumberOfElements(Fcn_X) == numberofelem*numberofqnodes &&
			mxGetNumberOfElements(Fcn_Y) == numberofelem*numberofqnodes ){
			// Fcn has numberofqnodes * numberofelem elements
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp += (Fcn_X_ptr[i*numberofqnodes + l]*
				           		(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
				            	 Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
				           	+Fcn_Y_ptr[i*numberofqnodes + l]*
				           		(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
				           	 	 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
				           	)*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp/2.0;
			}
		}
		else if (mxGetNumberOfElements(Fcn_X) == 1 &&
				 mxGetNumberOfElements(Fcn_Y) == 1 ) {
			// Fcn is constant
			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp += (*(Fcn_X_ptr)*
				           		(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
				            	 Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
				           	+*(Fcn_Y_ptr)*
				           		(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
				           	 	 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
				           	)*weights[l];
				}
				pLoad[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp/2.0;
			}
		}
		else {
			// Other cases does not match
			mexErrMsgTxt("Error:Assembler:AssembleGradLoad::Failed with unexpected Fcn.\n");
		}// end iff
	}//end for 
}//end all


void Assembler::AssembleGradXFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX, MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];
					}
					*pV  = *(Interp)*(*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}
}

void Assembler::AssembleGradYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV + Interp[i*numberofqnodes + l] *
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = (*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						*pV = *pV +
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					*pV  = *(Interp)*(*pV)/2.0;
					pI++; pJ++; pV++;
				}
			}
		}//end for
	}
}


void Assembler::AssembleGradXYFunc(Real_t* &pI, Real_t* &pJ, Real_t* &pV,Real_t* &pW,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr Ref, MatlabPtr RefX,MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr Fcn_X, MatlabPtr Fcn_Y){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp_X             = mxGetPr(Fcn_X);
	auto  Interp_Y             = mxGetPr(Fcn_Y);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];


	// evaluate Fcn at all quadrature nodes before calculation
	if (mxGetNumberOfElements(Fcn_X) == numberofelem*numberofqnodes &&
			mxGetNumberOfElements(Fcn_Y) == numberofelem*numberofqnodes ){

		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					*pW = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){


						*pV = *pV + Interp_X[i*numberofqnodes + l] *
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];

						*pW = *pW + Interp_Y[i*numberofqnodes + l] *
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];

					}
					*pV  = (*pV)/2.0;
					*pW  = (*pW)/2.0;
					pI++; pJ++; pV++;pW++;
				}
			}
		}//end for
	}//end if
	else {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;


			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];

			area = 0.5*fabs(det);


			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < numberofnodesperelem; k++){
					*pI = pelem_ptr[i*numberofnodesperelem + j];
					*pJ = pelem_ptr[i*numberofnodesperelem + k];
					*pV = 0.;
					*pW = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){


						*pV = *pV +
						           (Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
						            Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])
									* reference[k+ l*numberofnodesperelem]*
									weights[l];

						*pW = *pW +
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
								 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])
								* reference[k+ l*numberofnodesperelem]*
								weights[l];

					}
					*pV  = *(Interp_X)*(*pV)/2.0;
					*pW  = *(Interp_Y)*(*pW)/2.0;
					pI++; pJ++; pV++;pW++;
				}
			}
		}//end for
	}
}//end all



// build matrix form of load vector.
void Assembler::AssembleLoadMatrix(Real_t*& pI, Real_t*& pJ, Real_t*& pV, MatlabPtr Nodes,
		MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr Fcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
//	auto  qnodes_ptr           = mxGetPr(QNodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(Fcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area, tmp;


	auto Fcn_ptr = Matlab_Cast<Real_t>(Fcn);
	// linear interpolation
	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		if(mxGetNumberOfElements(Fcn) == numberofqnodes * numberofelem) {

			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Fcn_ptr[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem]*weights[l];
				}

				// which node
				*pI = pelem_ptr[i*numberofnodesperelem + j];
				// which element
				*pJ = (i + 1);
				// value
				*pV = tmp * area;

				pI++;pJ++;pV++;

			}
		}
		else {
			mexErrMsgTxt("AssemblerExtension::AssembleLoadMatrix::Dimension does not match.\n");
		}
	}
}

void Assembler::AssembleOverElement(Real_t*& w, MatlabPtr Nodes, MatlabPtr Elems,
		MatlabPtr Ref, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn_S, MatlabPtr Fcn_A, MatlabPtr u,
		MatlabPtr v){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);

	auto  Interp_S             = mxGetPr(Fcn_S);
	auto  Interp_A             = mxGetPr(Fcn_A);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);
	auto ptru                   = mxGetPr(u);
	auto ptrv                   = mxGetPr(v);


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(Fcn_S)  == numberofelem &&
			mxGetNumberOfElements(Fcn_A) == numberofelem) {
		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);


			auto a = Interp_A[i];
			auto s = Interp_S[i];

			int32_t I, J;
			Real_t K1, K2, K;


			w[i] = 0;
			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					I = pelem_ptr[i*numberofnodesperelem + j] - 1;
					J = pelem_ptr[i*numberofnodesperelem + k] - 1;
					K1 = 0.;
					K2 = 0.;
					for (size_t l = 0; l < numberofqnodes; l++){
						K1 = K1 + (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*weights[l];
						K2 = K2 +
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*
								weights[l];
					}
					K1 = s*(K1)/4.0/area;
					K2 = a*(K2) * area;
					K  = K1 + K2;

					w[i] += K * ptru[I] * ptrv[J];
					if (j != k){
						w[i] += K * ptru[J] * ptrv[I];
					}
				}
			}
		}
	}
}


void Assembler::AssembleOverNode(Real_t*& w, MatlabPtr Nodes, MatlabPtr Elems,
		MatlabPtr Ref, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr Fcn_S, MatlabPtr Fcn_A,
		MatlabPtr u, MatlabPtr v) {

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);

	auto  Interp_S             = mxGetPr(Fcn_S);
	auto  Interp_A             = mxGetPr(Fcn_A);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofnodes          = mxGetN(Nodes);
	auto numberofqnodes         = mxGetN(RefX);
	auto ptru                   = mxGetPr(u);
	auto ptrv                   = mxGetPr(v);


	mwSize vertex_1, vertex_2, vertex_3;
	Real_t det, area;
	Real_t Jacobian[2][2];

	// for performance purpose, not to use fill
	memset(w, 0,sizeof(Real_t) * numberofnodes);

	if (mxGetNumberOfElements(Fcn_S)  == numberofnodes &&
			mxGetNumberOfElements(Fcn_A) == numberofnodes) {

		for (size_t i =0; i < numberofelem; i++){
			// Fcn is a constant
			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			int32_t I, J, K;
			Real_t K1, K2, K3;

			for (size_t j = 0; j < numberofnodesperelem; j++){
				for (size_t k = 0; k < j + 1; k++){
					I = pelem_ptr[i*numberofnodesperelem + j] - 1;
					J = pelem_ptr[i*numberofnodesperelem + k] - 1;

					/*
					 * calculate (I,J) element
					 */

					/*
					 * (I,J) element uses all qnodes.
					 */
					for (size_t l = 0; l < numberofqnodes; l++) {
						/*
						 * l-th component
						 */
						K1 =(
							(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] +
							 Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] +
							 Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
							+
							(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] +
							 Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
							(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] +
							 Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
							)*weights[l];
						K2 = reference[j+ l*numberofnodesperelem]*
							 reference[k+ l*numberofnodesperelem]*
							 weights[l];
						for (size_t m = 0; m < numberofnodesperelem; m++) {
							// m th node in this element
							K = pelem_ptr[numberofnodesperelem * i + m] - 1;
							K3  = reference[m + l * numberofnodesperelem] *
									(Interp_S[K] * (K1/4.0/area) + Interp_A[K] * K2 * area);
							w[K] += K3 * ptru[I] * ptrv[J];
							if (j != k) {
								w[K] += K3 * ptru[J] * ptrv[I];
							}
						}// end of m
					}// end of l
				}// end of k
			}//end of j
		}//end of i
	}
}

void Assembler::AssembleMassEnergy(Real_t* &MassEnergy,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr QFcn, MatlabPtr PFcn){

	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(QFcn);
	auto  Pnodes               = mxGetPr(PFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);


	mwSize vertex_1, vertex_2 , vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;


	if (mxGetNumberOfElements(QFcn) == numberofelem*numberofqnodes &&
		mxGetNumberOfElements(PFcn) == numberofpnodes){
		*MassEnergy = 0.;
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			// Due to symmetric property, only need half of the work load.
			tmp = 0.;
			for (size_t j = 0; j < numberofnodesperelem; j++){
				rInd = pelem_ptr[i*numberofnodesperelem + j] - 1;

				for (size_t k = 0; k < numberofnodesperelem; k++){
					cInd = pelem_ptr[i*numberofnodesperelem + k] - 1;

					for (size_t l = 0; l < numberofqnodes; l++){
						qInd = i*numberofqnodes + l;

						tmp += 	Interp[qInd]*weights[l]*
								Pnodes[rInd]*reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*Pnodes[cInd];								
					}
				}
			}
			*MassEnergy = *MassEnergy + tmp*area;
		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");
	}
}

void Assembler::AssembleStiffEnergy(Real_t* &StiffEnergy,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr QFcn, MatlabPtr PFcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(QFcn);
	auto  Pnodes               = mxGetPr(PFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);
	auto numberofpnodes         = mxGetN(Nodes);

	mwSize vertex_1, vertex_2, vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(QFcn) == numberofelem*numberofqnodes &&
		mxGetNumberOfElements(PFcn) == numberofpnodes) {
		// Fcn is a matrix
		*StiffEnergy = 0.;
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			tmp = 0.;
			for (size_t j = 0; j < numberofnodesperelem; j++){
				rInd = pelem_ptr[i*numberofnodesperelem + j] - 1;

				for (size_t k = 0; k < numberofnodesperelem; k++){
					cInd = pelem_ptr[i*numberofnodesperelem + k] - 1;

					for (size_t l = 0; l < numberofqnodes; l++){
						qInd = i*numberofqnodes + l;

						tmp += 	Interp[qInd]*weights[l]*
								Pnodes[rInd]* (
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*Pnodes[cInd];
					}
				}
			}
			*StiffEnergy = *StiffEnergy + tmp/4.0/area;

		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");	
	}
}

void Assembler::AssembleMassEnergyGrad(Real_t* &MassEnergyGrad,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr QFcn, MatlabPtr PFcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(QFcn);
	auto  Pnodes               = mxGetPr(PFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);


	mwSize vertex_1, vertex_2 , vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;


	if (mxGetNumberOfElements(QFcn) == numberofelem*numberofqnodes &&
		mxGetNumberOfElements(PFcn) == numberofpnodes){
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			// Due to symmetric property, only need half of the work load.
			for (size_t j = 0; j < numberofnodesperelem; j++){
				rInd = pelem_ptr[i*numberofnodesperelem + j] - 1;

				tmp = 0.;
				for (size_t k = 0; k < numberofnodesperelem; k++){
					cInd = pelem_ptr[i*numberofnodesperelem + k] - 1;

					for (size_t l = 0; l < numberofqnodes; l++){
						qInd = i*numberofqnodes + l;

						tmp += 	Interp[qInd]*weights[l]*
								reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*Pnodes[cInd];								
					}
				}
				MassEnergyGrad[rInd] +=  tmp*area;
			}			
		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");
	}
}

void Assembler::AssembleStiffEnergyGrad(Real_t* &StiffEnergyGrad,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX,
		MatlabPtr RefY, MatlabPtr Weights, MatlabPtr QFcn, MatlabPtr PFcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  Interp               = mxGetPr(QFcn);
	auto  Pnodes               = mxGetPr(PFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);
	auto numberofpnodes         = mxGetN(Nodes);

	mwSize vertex_1, vertex_2, vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(QFcn) == numberofelem*numberofqnodes &&
		mxGetNumberOfElements(PFcn) == numberofpnodes) {
		// Fcn is a matrix
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t j = 0; j < numberofnodesperelem; j++){
				rInd = pelem_ptr[i*numberofnodesperelem + j] - 1;

				tmp = 0.;
				for (size_t k = 0; k < numberofnodesperelem; k++){
					cInd = pelem_ptr[i*numberofnodesperelem + k] - 1;

					for (size_t l = 0; l < numberofqnodes; l++){
						qInd = i*numberofqnodes + l;

						tmp += 	Interp[qInd]*weights[l]*(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*Pnodes[cInd];
					}
				}
				StiffEnergyGrad[rInd] += tmp/4.0/area;
			}
			

		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");	
	}
}

void Assembler::AssembleMassGrad(Real_t* &MassGrad,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr rPFcn, MatlabPtr cPFcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  rPnodes              = mxGetPr(rPFcn);
 	auto  cPnodes              = mxGetPr(cPFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);


	mwSize vertex_1, vertex_2 , vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;


	if (mxGetNumberOfElements(rPFcn) == numberofpnodes &&
		mxGetNumberOfElements(cPFcn) == numberofpnodes){
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			for (size_t l = 0; l < numberofqnodes; l++){
				tmp = 0.;
				qInd = i*numberofqnodes + l;

				for (size_t j = 0; j < numberofnodesperelem; j++){
					for (size_t k = 0; k < numberofnodesperelem; k++){
						rInd = pelem_ptr[i*numberofnodesperelem + j]-1;
						cInd = pelem_ptr[i*numberofnodesperelem + k]-1;
						
						tmp += 	rPnodes[rInd]*reference[j+ l*numberofnodesperelem]*
								reference[k+ l*numberofnodesperelem]*cPnodes[cInd];
								
					}
				}
				MassGrad[qInd] = tmp * weights[l] * area;
			}
		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");
	}
}

void Assembler::AssembleStiffGrad(Real_t* &StiffGrad,
		MatlabPtr Nodes, MatlabPtr Elems, MatlabPtr RefX, MatlabPtr RefY,
		MatlabPtr Weights, MatlabPtr rPFcn, MatlabPtr cPFcn) {


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  referenceX           = mxGetPr(RefX);
	auto  referenceY           = mxGetPr(RefY);
	auto  weights              = mxGetPr(Weights);
	auto  rPnodes              = mxGetPr(rPFcn);
 	auto  cPnodes              = mxGetPr(cPFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(RefX);
	auto numberofpnodes         = mxGetN(Nodes);

	mwSize vertex_1, vertex_2, vertex_3, rInd, cInd, qInd;
	Real_t det, area, tmp;
	Real_t Jacobian[2][2];

	if (mxGetNumberOfElements(rPFcn) == numberofpnodes &&
		mxGetNumberOfElements(cPFcn) == numberofpnodes) {
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			Jacobian[0][0] = pnodes_ptr[2*vertex_3 + 1] - pnodes_ptr[2*vertex_1 + 1];
			Jacobian[1][1] = pnodes_ptr[2*vertex_2    ] - pnodes_ptr[2*vertex_1    ];
			Jacobian[0][1] = pnodes_ptr[2*vertex_1 + 1] - pnodes_ptr[2*vertex_2 + 1];
			Jacobian[1][0] = pnodes_ptr[2*vertex_1    ] - pnodes_ptr[2*vertex_3    ];

			// Orientation corrected.
			det = Jacobian[0][0] * Jacobian[1][1] - Jacobian[0][1] * Jacobian[1][0];
			area = 0.5*fabs(det);

			// Due to symmetric property, half of work load can be reduced
			for (size_t l = 0; l < numberofqnodes; l++){
				tmp = 0;
				qInd = i*numberofqnodes + l;

				for (size_t j = 0; j < numberofnodesperelem; j++){
					for (size_t k = 0; k < j + 1; k++){
						rInd = pelem_ptr[i*numberofnodesperelem + j]-1;
						cInd = pelem_ptr[i*numberofnodesperelem + k]-1; 						
							 
						tmp = tmp + rPnodes[rInd] *(
								(Jacobian[0][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[0][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[0][1]*referenceY[k+ l*numberofnodesperelem])
								+
								(Jacobian[1][0]*referenceX[j+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[j+ l*numberofnodesperelem])*
								(Jacobian[1][0]*referenceX[k+ l*numberofnodesperelem] + Jacobian[1][1]*referenceY[k+ l*numberofnodesperelem])
								)*cPnodes[cInd];
					
					}
				}
				StiffGrad[qInd] = tmp * weights[l] /4.0 / area;
			}
		}
	}else{
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleMassGrad::Failed with unexpected Fcn.\n");
	}
}

void Assembler::AssembleQ2IntTrans(Real_t* &qWeights,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);


	mwSize vertex_1, vertex_2 , vertex_3;
	Real_t det, area;


	for (size_t i =0; i < numberofelem; i++){

		vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
		vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
		vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

		det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
				(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
		area = 0.5*fabs(det);

		for (size_t l = 0; l < numberofqnodes; l++){
			qWeights[i*numberofqnodes + l] += weights[l]*area;								
		}
	}

}

void Assembler::AssembleLoadTrans(Real_t* &qLoad,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr Weights, MatlabPtr PFcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  weights              = mxGetPr(Weights);
	auto  Pnodes              = mxGetPr(PFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);


	mwSize vertex_1, vertex_2 , vertex_3, pInd, qInd;
	Real_t det, area, tmp;


	if (mxGetNumberOfElements(PFcn) == numberofpnodes){
		for (size_t i =0; i < numberofelem; i++){

			vertex_1 = pelem_ptr[numberofnodesperelem*i] - 1;
			vertex_2 = pelem_ptr[numberofnodesperelem*i + 1] - 1;
			vertex_3 = pelem_ptr[numberofnodesperelem*i + 2] - 1;

			det = (pnodes_ptr[vertex_2*2] - pnodes_ptr[vertex_1*2])*(pnodes_ptr[vertex_3*2 + 1] - pnodes_ptr[vertex_1*2 + 1]) -
					(pnodes_ptr[vertex_2*2 + 1] - pnodes_ptr[vertex_1*2 + 1])*(pnodes_ptr[vertex_3*2] - pnodes_ptr[vertex_1*2]);
			area = 0.5*fabs(det);

			// Due to symmetric property, only need half of the work load.
			for (size_t l = 0; l < numberofqnodes; l++) {
				tmp = 0.;
				qInd = i*numberofqnodes + l;

				for (size_t j = 0; j < numberofnodesperelem; j++) {
					pInd = 	pelem_ptr[i*numberofnodesperelem + j] - 1;
					tmp +=  Pnodes[pInd]*reference[j+ l*numberofnodesperelem];

				}
				qLoad[qInd] = tmp*area*weights[l];
			}

		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleLoadTrans::Failed with unexpected Fcn.\n");
	}
}

void Assembler::AssembleP2QTrans(Real_t* &PFcn,
		MatlabPtr Nodes, MatlabPtr Elems,MatlabPtr Ref,
		MatlabPtr QFcn){


	auto  pnodes_ptr           = mxGetPr(Nodes);
	auto  pelem_ptr            = (int32_t*)mxGetPr(Elems);
	auto  reference            = mxGetPr(Ref);
	auto  Qnodes              = mxGetPr(QFcn);

	auto numberofelem           = mxGetN(Elems);
	auto numberofnodesperelem   = mxGetM(Elems);
	auto numberofqnodes         = mxGetN(Ref);
	auto numberofpnodes         = mxGetN(Nodes);

	mwSize pInd, qInd;
	Real_t tmp;


	if (mxGetNumberOfElements(QFcn) == numberofqnodes * numberofelem){
		for (size_t i =0; i < numberofelem; i++){

			for (size_t j = 0; j < numberofnodesperelem; j++) {
				tmp = 0.;
				for (size_t l = 0; l < numberofqnodes; l++) {
					tmp +=  Qnodes[i*numberofqnodes + l]*reference[j+ l*numberofnodesperelem];
				}
				PFcn[pelem_ptr[i*numberofnodesperelem + j] - 1] += tmp;
			}

		}
	}else {
		// Other cases does not match
		mexErrMsgTxt("Error:Assembler:AssembleP2QTrans::Failed with unexpected Fcn.\n");
	}
}


template class mexplus::Session<Assembler>;


namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<Assembler>::create(new Assembler()));
}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<Assembler>::destroy(input.get(0));
}

MEX_DEFINE(reference2D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofpoints = mxGetN(prhs[1]);
	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);

	assembler->Reference(plhs[0], plhs[1], plhs[2], CAST(prhs[1]), CAST(prhs[2]));


}

MEX_DEFINE(reference1D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 3);
	OutputArguments output(nlhs, plhs, 2);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofpoints = static_cast<size_t>(*Matlab_Cast<Real_t>(CAST(prhs[1])) + 1);

	size_t numberofqpoints = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofpoints, numberofqpoints, mxDOUBLE_CLASS, mxREAL);
	assembler->Reference(plhs[0], plhs[1], CAST(prhs[1]), CAST(prhs[2]));

}

MEX_DEFINE(assems)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler->AssembleStiff(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));
}

// advanced interface with matrix kernel
MEX_DEFINE(assemsm)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler->AssembleStiff(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));
}

MEX_DEFINE(assema)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler->AssembleMass(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]));


}

MEX_DEFINE(assemsa)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 9);
	OutputArguments output(nlhs, plhs, 4);
	Assembler* assembler = Session<Assembler>::get(input.get(0));
	/*
	 * do both at same time, if it is needed.(always do).
	 */
}

MEX_DEFINE(asseml) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pLoad = mxGetPr(plhs[0]);
	assembler->AssembleLoad(pLoad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));

}

MEX_DEFINE(qnodes2D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofelem           = mxGetN(prhs[3]);
	size_t numberofqnodes         = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(2, numberofelem * numberofqnodes,  mxDOUBLE_CLASS, mxREAL);
	Real_t* Coords = mxGetPr(plhs[0]);
	assembler->Qnodes2D(Coords,CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]));
}

MEX_DEFINE(qnodes1D)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 4);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));


	size_t numberofedges          = mxGetN(prhs[3]);
	size_t numberofqnodes         = mxGetN(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[1]), numberofedges * numberofqnodes,  mxDOUBLE_CLASS, mxREAL);
	Real_t* Coords = mxGetPr(plhs[0]);
	assembler->Qnodes1D(Coords,CAST(prhs[1]), CAST(prhs[2]),
			CAST(prhs[3]));
}

MEX_DEFINE(assemrbc) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pNeumann = mxGetPr(plhs[0]);
	assembler->AssembleBC(pNeumann, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));

}

MEX_DEFINE(assemlbc) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes         = mxGetN(prhs[1]);
	size_t numberofedges         = mxGetN(prhs[2]);
	size_t numberofnodesperedges = mxGetM(prhs[2]);

	plhs[0] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericMatrix(numberofedges * numberofnodesperedges * numberofnodesperedges, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);
	Real_t* pJ = mxGetPr(plhs[1]);
	Real_t* pV = mxGetPr(plhs[2]);
	assembler->AssembleBC(pI, pJ, pV,  CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]));

}

MEX_DEFINE(assemgradl) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 9);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   = mxGetN(prhs[1]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes,1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pLoad = mxGetPr(plhs[0]);
	
	assembler->AssembleGradLoad(pLoad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]), CAST(prhs[8]));
}

MEX_DEFINE(assemex_gradfunc_x) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradXFunc(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));
}

MEX_DEFINE(assemex_gradfunc_y)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleGradYFunc(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));
}

MEX_DEFINE(assemex_gradfunc_xy)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 9);
	OutputArguments output(nlhs, plhs, 4);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);


	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	plhs[3] = mxCreateNumericMatrix(numberofnodesperelem * numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pW = mxGetPr(plhs[3]);

	assembler_ex->AssembleGradXYFunc(pI, pJ, pV, pW, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]),CAST(prhs[8]));
}

MEX_DEFINE(assemex_massenergy) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , 1 ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* MassEnergy = mxGetPr(plhs[0]);
	assembler->AssembleMassEnergy(MassEnergy, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));

}

MEX_DEFINE(assemex_stiffenergy) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , 1 ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* StiffEnergy = mxGetPr(plhs[0]);
	assembler->AssembleStiffEnergy(StiffEnergy, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));

}

MEX_DEFINE(assemex_massenergygrad) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   		= mxGetN(prhs[1]);
	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( numberofnodes , 1 ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* MassEnergyGrad = mxGetPr(plhs[0]);
	assembler->AssembleMassEnergyGrad(MassEnergyGrad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));

}

MEX_DEFINE(assemex_stiffenergygrad) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes   		= mxGetN(prhs[1]);
	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( numberofnodes , 1 ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* StiffEnergyGrad = mxGetPr(plhs[0]);
	assembler->AssembleStiffEnergyGrad(StiffEnergyGrad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));

}

MEX_DEFINE(assemex_massgrad) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , numberofqnodes * numberofelem ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* MassGrad = mxGetPr(plhs[0]);
	assembler->AssembleMassGrad(MassGrad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]));

}

MEX_DEFINE(assemex_stiffgrad) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 8);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , numberofqnodes * numberofelem ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* StiffGrad = mxGetPr(plhs[0]);
	assembler->AssembleStiffGrad(StiffGrad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]));

}

MEX_DEFINE(assemex_q2itrans) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 5);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , numberofqnodes * numberofelem ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* qWeights = mxGetPr(plhs[0]);
	assembler->AssembleQ2IntTrans(qWeights, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]));

}

MEX_DEFINE(assemex_loadtrans) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodesperelem = mxGetM(prhs[2]);
	size_t numberofelem         = mxGetN(prhs[2]);
	size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( 1 , numberofqnodes * numberofelem ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* qLoad = mxGetPr(plhs[0]);
	assembler->AssembleLoadTrans(qLoad, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]));

}

MEX_DEFINE(assemex_p2qtrans) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 5);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes        = mxGetN(prhs[1]);
	// size_t numberofnodesperelem = mxGetM(prhs[2]);
	// size_t numberofelem         = mxGetN(prhs[2]);
	// size_t numberofqnodes       = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix( numberofnodes , 1 ,  mxDOUBLE_CLASS, mxREAL);
	Real_t* PFcn = mxGetPr(plhs[0]);
	assembler->AssembleP2QTrans(PFcn, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]));

}


MEX_DEFINE(assemex_lm) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 6);
	OutputArguments output(nlhs, plhs, 3);
	Assembler* assembler_ex = Session<Assembler>::get(input.get(0));

	auto numberofnodesperelem = mxGetM(prhs[2]);
	auto numberofelem    =  mxGetN(prhs[2]);


//	MatlabPtr Nodes,
//			MatlabPtr QNodes, MatlabPtr Elems,MatlabPtr Ref,
//			MatlabPtr Weights, MatlabPtr Fcn

	plhs[0] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem , 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pI = mxGetPr(plhs[0]);

	plhs[1] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem , 1,  mxDOUBLE_CLASS, mxREAL);
	Real_t* pJ = mxGetPr(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(numberofnodesperelem * numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* pV = mxGetPr(plhs[2]);

	assembler_ex->AssembleLoadMatrix(pI, pJ, pV, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]));

}

MEX_DEFINE(assemloe)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 11);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes          = mxGetN(prhs[1]);
	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofelem, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* w = mxGetPr(plhs[0]);

	assembler->AssembleOverElement(w, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]), CAST(prhs[8]),
			CAST(prhs[9]), CAST(prhs[10]));

}


MEX_DEFINE(assemlon)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 11);
	OutputArguments output(nlhs, plhs, 1);
	Assembler* assembler = Session<Assembler>::get(input.get(0));

	size_t numberofnodes          = mxGetN(prhs[1]);
	size_t numberofelem           = mxGetN(prhs[2]);
	size_t numberofnodesperelem   = mxGetM(prhs[2]);
	size_t numberofqnodes         = mxGetM(prhs[3]);

	plhs[0] = mxCreateNumericMatrix(numberofnodes, 1, mxDOUBLE_CLASS, mxREAL);
	Real_t* w = mxGetPr(plhs[0]);

	assembler->AssembleOverNode(w, CAST(prhs[1]),
			CAST(prhs[2]), CAST(prhs[3]),
			CAST(prhs[4]), CAST(prhs[5]),
			CAST(prhs[6]), CAST(prhs[7]), CAST(prhs[8]),
			CAST(prhs[9]), CAST(prhs[10]));

}
}


MEX_DISPATCH


