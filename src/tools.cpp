#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth){
  //cout<<"three"<<endl;
	VectorXd rmse(4);
	rmse << 0,0,0,0;

  double px_sqrd_err;
	double py_sqrd_err;
	double vx_sqrd_err;
	double vy_sqrd_err;

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0)
			{
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
			}
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); i++){
        // ... your code here
        //cout<<px_sqrd_err + pow((estimations[i][0] - ground_truth[i][0]),2)<<endl;
        px_sqrd_err = px_sqrd_err + pow((estimations[i][0] - ground_truth[i][0]),2);
				py_sqrd_err = py_sqrd_err + pow((estimations[i][1] - ground_truth[i][1]),2);
				vx_sqrd_err = vx_sqrd_err + pow((estimations[i][2] - ground_truth[i][2]),2);
				vy_sqrd_err = vy_sqrd_err + pow((estimations[i][3] - ground_truth[i][3]),2);
	}

	//calculate the mean
	// ... your code here
    rmse <<  px_sqrd_err / estimations.size(),
						py_sqrd_err / estimations.size(),
						vx_sqrd_err / estimations.size(),
						vy_sqrd_err / estimations.size();


	//calculate the squared root
	// ... your code here
	for(int n=0; n < rmse.size(); n++){

	    rmse[n] = sqrt(rmse[n]);

	}
	/*for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];
		//cout<<"Estimation"<<estimations[i]<<endl;
		//cout<<"G Truth"<<ground_truth[i]<<endl;
		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();
	//cout <<"RMSE Size"<<estimations.size()<<endl;
	//calculate the squared root
	rmse = rmse.array().sqrt();*/

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//TODO: YOUR CODE HERE

	//check division by zero
	if(((pow((px*px+ py*py),1.5)) == 0) || ((sqrt(px*px+ py*py)) == 0))
	{
	    cout<<"Devision by zero while calculating Jacobian elements"<<endl;
	    return Hj= Eigen::MatrixXd::Identity(3, 4);
	}

	else
	{
	    //compute the Jacobian matrix
        Hj << px/(sqrt(px*px+ py*py)), py/(sqrt(px*px+ py*py)), 0, 0,
        -py/(px*px+ py*py), px/(px*px+ py*py), 0, 0,
        py*(vx*py - vy*px) / (pow((px*px+ py*py),1.5)), px*(vy*px - vx*py) / (pow((px*px+ py*py),1.5)), px/(sqrt(px*px+ py*py)), py/(sqrt(px*px+ py*py));


	    return Hj;
	}

}
