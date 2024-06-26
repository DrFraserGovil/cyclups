#pragma once
#include "../Assumptions/Assumptions.h"
#include "prediction.h"
#include "Eigen"
#include "../Constraints/Constraints.h"
#include "OptimiserProperties.h"
namespace cyclups
{
	


	struct Container
	{
		// Vector X;
		std::vector<Vector> a_blups;
		std::vector<double> p_blps;
		std::vector<double> p_blups;
		std::vector<Vector> ks;
		Vector Delta;
		double Beta;
		Vector Bp_blups;
		Eigen::LDLT<Matrix> BBt; 
		Matrix K;
		Eigen::PartialPivLU<Matrix> Binv;
		Eigen::PartialPivLU<Matrix> B_trans_inv;
		Vector p;
		double Prefactor;
	};

	

	class Predictor
	{
		public:
			OptimiserProperties Optimiser;


			//needs to be templated because it can accept either a ConstraintSet or a Constraint (or any subclass thereof).
			template<class T>
			Predictor(kernel::Kernel k, basis::Basis b, T &constraint) : Kernel(k), Basis(b), Constraint(constraint){};

		
			Prediction Predict(cvec predictX, const PairedData & data);

			Prediction RegularisedPrediction(cvec predictX, const PairedData & data, RegularisingFunction Regulariser);

			//the Retire function kills off the Constraint object, allowing its data to be returned to its original owner, and hence used (i.e. in a new predictor). This calls the Destructor of the Constraint object without forcing the Predictor to go out of scope. 
			void Retire();
		private:
			Container Store; // a place to put lots of useful storage variables

			kernel::Kernel Kernel;
			basis::Basis Basis;
			constraint::ConstraintSet Constraint;
			RegularisingFunction R;
			bool UsingRegulariser = false;
			void Initialise(cvec PredictX, const PairedData & data);
			void Optimise(cvec predictX, const PairedData & data);
			double ComputeScore(cvec predictX);
			void BulkUp(cvec predictX);
	};
}