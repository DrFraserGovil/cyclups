#pragma once
#include "../ConstraintSet.h"
namespace cyclups::constraint
{
	class GreaterThan : public Constraint
	{
		public:
			GreaterThan(double value)
			{
				usingConst = true;
				constVal = value;
			}
			GreaterThan(double (* fnc)(double) )
			{
				usingFunction = true;
				fncVal = fnc;
			}
			GreaterThan(double value, bool (* domain)(double))
			{
				usingConst = true;
				constVal = value;
				usingDomain = true;
				inDomain = domain;
			}
			GreaterThan(double (* fnc)(double), bool (* domain)(double))
			{
				usingFunction = true;
				fncVal = fnc;
				usingDomain = true;
				inDomain = domain;
			}
		private:
			
			double (* fncVal)(double);
			bool usingFunction = false;
			bool usingConst = false;
			double constVal;
			

			InitialiseContainer Initialiser(cvec t)
			{
				cvec domainT = ApplyDomain(t);

				int n = domainT.size();
				Matrix B;
				if (domainT.size() == t.size())
				{
					B = Matrix::Identity(n,n);
				}
				else
				{
					B = Matrix::Zero(n,t.size());
					int q = 0;
					for (int i = 0; i < t.size(); ++i)
					{
						if (q< domainT.size() && t[i] == domainT[q])
						{
							B(q,i) = 1;
							++q;
							break;
						}
					}
				}
				
				transformOperator f = [](Vector & output, const Vector & input,std::vector<double> & params){for (int i =0; i < output.size(); ++i){
					output[i] = params[i] + exp(input[i]);
				}};
				transformOperator grad= [](Vector & output, const Vector & input,std::vector<double> & params){for (int i =0; i < output.size(); ++i){
					output[i] = exp(input[i]);
				}};
				transformOperator inv = [](Vector & output, const Vector & input,std::vector<double> & params){
					for (int i =0; i < output.size(); ++i){
					double buffered = std::max(input[i] - params[i],1e-2);
					output[i] = log(buffered);
				}};
				ConstraintVector c = ConstraintVector::Optimise(n,n,SeparableTransform(f,grad,inv));
				if (usingConst)
				{
					std::vector<double> v (n,constVal);
					c.SetParams(v);
				}
				if (usingFunction)
				{
					std::vector<double> v(n);
					for (int i = 0; i < n; ++i)
					{
						v[i] = fncVal(domainT[i]);
					}
					c.SetParams(v);
				}
				return InitialiseContainer(c,B); 
			}
	};

	inline GreaterThan Positive()
	{
		return GreaterThan(0.); //has to be 0. as just 0 is the Null pointer which causes ambiguous conversion with the function pointer option
	}
	inline GreaterThan Positive(bool (*domain)(double))
	{
		return GreaterThan(0.,domain); //has to be 0. as just 0 is the Null pointer which causes ambiguous conversion with the function pointer option
	}
}