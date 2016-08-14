#include "LET.h"
#include <mexplus.h>
#include "utils.h"

using namespace std;
using namespace mexplus;



static std::vector<double> X {
        -0.90617984593866396370032134655048,
        -0.53846931010568310771446931539685,
        0,
        0.53846931010568310771446931539685,
        0.90617984593866396370032134655048

};
static std::vector<double> W {
        0.23692688505618908489935847683228,
        0.47862867049936647090291330641776,
        0.56888888888888888888888888888889,
        0.47862867049936647090291330641776,
        0.23692688505618908489935847683228
};
inline double distance(double x0, double y0, double x1, double y1) noexcept {
    return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

inline double integral_corner(double sigma_t, double length) {
    double ret = 0.;
    double mid = (sqrt(2) + 1.0) * length/4.0;
    double rad = (sqrt(2) - 1.0) * length/4.0;
    for (int i = 0; i < X.size(); ++i) {
        double val = (mid + rad * X[i]);
        ret += exp(-sigma_t * val) * (M_PI/2.0 - 2 * acos(length/2.0/val)) * W[i];
    }
    ret *= 4.0 * rad;
    return ret;
}

class kernel_kifmm: public LET {
public:
    // involves some grid data to locate information
    unsigned long side;
    double start_x;
    double start_y;
    double *mu_t;

    double integral_helper( double x0, double y0, double x1, double y1) noexcept {
        return distance(x0, y0, x1, y1) * this->getAttribute((x0 + x1)/2 , (y0 + y1)/2);
    }

    size_t getRow(double y) noexcept {
        return size_t(floor(side * y));
    }

    size_t getCol(double x) noexcept {
        return size_t(floor(side * x));
    }

    double integral(double x0, double y0, double x1, double y1) noexcept {
        auto col0 = getRow(x0);
        auto col1 = getRow(x1);
        auto row0 = getCol(y0);
        auto row1 = getCol(y1);

        // 9 cases
        if ((row0 == row1) && (col0 == col1)) {
            return integral_helper(x0, y0, x1, y1);
        }
        else if ((row0 == row1 + 1) && (col0 == col1)) {
            auto row0 = getCol(y0);
            double ybar = double(row0)/side;
            double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
            return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
        }
        else if ((row0 == row1 - 1) && (col0 == col1)) {
            auto row1 = getCol(y1);
            double ybar = double(row1)/side;
            double xbar = ((y1 - ybar) * x0 + (ybar - y0) * x1)/(y1 - y0);
            return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
        }
        else if ((col0 == col1 + 1) && (row0 == row1)) {
            auto col0 = getCol(x0);
            double xbar = double(col0)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
        }
        else if ((col0 == col1 - 1) && (row0 == row1)) {
            auto col1 = getCol(x1);
            double xbar = double(col1)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            return integral_helper(x0, y0, xbar, ybar) + integral_helper(xbar, ybar, x1, y1);
        }
        else if ((col0 == col1 + 1) && (row0 == row1 + 1)) {
            double xbar = double(col0)/side;
            double ybar2 = double(row0)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

            if (xbar < xbar2) {
                return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
            }
            else {
                return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

            }
        }
        else if ((col0 == col1 + 1) && (row0 == row1 - 1)) {
            double xbar = double(col0)/side;
            double ybar2 = double(row1)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

            if (xbar < xbar2) {
                return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
            }
            else {
                return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

            }

        }
        else if ((col0 == col1 - 1) && (row0 == row1 + 1)) {
            double xbar = double(col1)/side;
            double ybar2 = double(row0)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

            if (xbar > xbar2) {
                return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
            }
            else {
                return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);
            }
        }
        else if ((col0 == col1 - 1) && (row0 == row1 - 1)) {
            double xbar = double(col1)/side;
            double ybar2 = double(row1)/side;
            double ybar = ((x1 - xbar) * y0 + (xbar - x0) * y1)/(x1 - x0);
            double xbar2 = ((y1 - ybar2) * x0 + (ybar2 - y0) * x1)/(y1 - y0);

            if (xbar > xbar2) {
                return integral_helper(x1, y1, xbar, ybar) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar2, ybar2, x0, y0);
            }
            else {
                return integral_helper(x1, y1, xbar2, ybar2) + \
						integral_helper(xbar, ybar, xbar2, ybar2) + integral_helper(xbar, ybar, x0, y0);

            }
        }
        else {
            double xm = (x0 + x1)/2;
            double ym = (y0 + y1)/2;
            return integral(x0, y0, xm, ym) + integral(xm, ym, x1, y1);
        }
    }

    double integral_approximate(double x0, double y0, double x1, double y1) noexcept {
        // use precise integral.
        if (distance(x0, y0, x1, y1) < 1.0/side) {
            return integral_helper(x0, y0, x1, y1);
        }
        else {
            auto xm = (x0 + x1)/2;
            auto ym = (y0 + y1)/2;
            return integral_approximate( x0, y0, xm, ym) +
                   integral_approximate( xm, ym, x1, y1);
        }
    }
    double getAttribute(double x, double y);
    virtual double kernelFunction(Point& r0, Point& r1);
};
double kernel_kifmm::kernelFunction(Point& r0, Point& r1){
    double rSquare	= distance(r0.x, r0.y, r1.x, r1.y) ;
    if (rSquare  == 0.){
//    	return 0.;
        return (1.0 - exp(- this->getAttribute(r0.x, r0.y)/side/2.0))/ this->getAttribute(r0.x, r0.y) +
               integral_corner(this->getAttribute(r0.x, r0.y), 1.0/side)/(2*M_PI);
    }
    else{
        return exp(-this->integral(r0.x, r0.y, r1.x, r1.y))/rSquare/(2*M_PI)/side/side;
    }

}

double kernel_kifmm::getAttribute(double x, double y) {
    if (x < start_x || x >= 1 + start_x || y < start_y || y >= 1 + start_y) {
        return 0.;
    }


    auto col = (unsigned long)((x - start_x) * side);
    auto row = (unsigned long)((y - start_y) * side);

    return mu_t[row * side + col];
}


template class mexplus::Session<kernel_kifmm>;

namespace {

MEX_DEFINE(new)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 0);
	OutputArguments output(nlhs, plhs, 1);
	output.set(0, Session<kernel_kifmm>::create(new kernel_kifmm()));

}

MEX_DEFINE(delete) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	Session<kernel_kifmm>::destroy(input.get(0));
}

MEX_DEFINE(calc_cache) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);

	auto kern = Session<kernel_kifmm>::get(input.get(0));

	auto ncheb_ptr = mxGetPr(prhs[1]);
	auto charges_ptr = mxGetPr(prhs[2]);

	auto location_ptr = mxGetPr(prhs[3]);
	auto numberoflocation = mxGetN(prhs[3]);

	auto N_ptr = mxGetPr(prhs[4]);
	auto rank_ptr = mxGetPr(prhs[5]);

	auto mu_t_ptr = mxGetPr(prhs[6]);

	vector<Point> location;

	location.resize(numberoflocation);

	for(size_t i = 0; i < numberoflocation; i++) {
		location[i].x = location_ptr[2 * i];
		location[i].y = location_ptr[2 * i + 1];
	}

	plhs[0] = mxCreateNumericMatrix((unsigned long)*N_ptr, 1 ,mxDOUBLE_CLASS, mxREAL);


	kern->side = (unsigned long)(sqrt(*N_ptr));
	kern->start_x = 0.;
	kern->start_y = 0.;
	kern->mu_t = mu_t_ptr;

	//kern->tic();
	kern->initialize((unsigned short)*ncheb_ptr , location, charges_ptr, (unsigned long)*N_ptr, (unsigned)*rank_ptr, Option::Cache);
	//kern->toc();


    //std::cout << "Up  : elapsed time at " << kern-> elapsed << " seconds"  << endl;
    kern->resetPotential(kern->root);
    VectorXd potentialMatrix;
    potentialMatrix = VectorXd::Zero(kern->N);

    //kern->tic();
    kern->downPass(kern->root, potentialMatrix);

    //kern->toc();
    //std::cout << "Down: elapsed time at " << kern-> elapsed << " seconds"  << endl;
    Map<VectorXd>(mxGetPr(plhs[0]), potentialMatrix.rows()) = potentialMatrix;

}

MEX_DEFINE(calc_fast) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	InputArguments input(nrhs, prhs, 7);
	OutputArguments output(nlhs, plhs, 1);

	auto kern = Session<kernel_kifmm>::get(input.get(0));

	auto ncheb_ptr = mxGetPr(prhs[1]);
	auto charges_ptr = mxGetPr(prhs[2]);

	auto location_ptr = mxGetPr(prhs[3]);
	auto numberoflocation = mxGetN(prhs[3]);

	auto N_ptr = mxGetPr(prhs[4]);
	auto rank_ptr = mxGetPr(prhs[5]);

	auto mu_t_ptr = mxGetPr(prhs[6]);

	vector<Point> location;

	location.resize(numberoflocation);

	for(size_t i = 0; i < numberoflocation; i++) {
		location[i].x = location_ptr[2 * i];
		location[i].y = location_ptr[2 * i + 1];
	}

	plhs[0] = mxCreateNumericMatrix((unsigned long)*N_ptr, 1 ,mxDOUBLE_CLASS, mxREAL);


	kern->side = (unsigned long)(sqrt(*N_ptr));
	kern->start_x = 0.;
	kern->start_y = 0.;
	kern->mu_t = mu_t_ptr;

//	kern->tic();
	kern->initialize((unsigned short)*ncheb_ptr , location, charges_ptr, (unsigned long)*N_ptr, (unsigned)*rank_ptr, Option::Fast);
//	kern->toc();


//    std::cout << "Up  : elapsed time at " << kern-> elapsed << " seconds"  << endl;
    kern->resetPotential(kern->root);
    VectorXd potentialMatrix;
    potentialMatrix = VectorXd::Zero(kern->N);

//    kern->tic();
    kern->downPass(kern->root, potentialMatrix);

//    kern->toc();
//    std::cout << "Down: elapsed time at " << kern-> elapsed << " seconds"  << endl;
    Map<VectorXd>(mxGetPr(plhs[0]), potentialMatrix.rows()) = potentialMatrix;

}

MEX_DEFINE(debug)  (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
	InputArguments input(nrhs, prhs, 1);
	OutputArguments output(nlhs, plhs, 0);
	auto kern = Session<kernel_kifmm>::get(input.get(0));

	std::cout << "==== Debug information ====" << endl << endl;
	std::cout << "= max level: "<< kern->maxLevel << endl;
	std::cout << "= root particles number: " << kern->N << endl;

	queue<pair<LET_Node*, vector<int>>> q;
	q.push({kern->root, {0}});
	while(!q.empty()) {
		// BFS
		auto cur = q.front();
		q.pop();
		int length = cur.second.size();
		cout << "= box path : ";
		for (int k = 0; k < length; ++k) {
			cout << cur.second[k] << " ";
		}

		for (int k = 0; k < kern->maxLevel - length + 1; ++k) {
			cout << "  ";
		}

		cout << "at position (" <<  cur.first->center.x <<", "<< cur.first->center.y << ") particles: " << cur.first->N << endl;

		for (int i = 0; i < 4; ++i) {
			if (cur.first->child[i] != nullptr) {
				vector<int> tmp = cur.second;
				tmp.push_back(cur.first->child[i]->nodeIndex);
				q.push({cur.first->child[i], tmp});
			}
		}
	}

}

}
MEX_DISPATCH
