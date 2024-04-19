#pragma once
#include <vector>
#include "JSL.h"

namespace cyclups
{
	//some simple datatypes for ferrying around XY coordinates in a recognisable manner

	struct Pair
	{
		double X;
		double Y;
		double Error;
		Pair(double x, double y)
		{
			X = x;
			Y = y;
			Error = 0;
		}
		Pair(double x, double y,double e)
		{
			X = x;
			Y = y;
			Error = e;
		}
	};

	class PairedData
	{
		public:
			std::vector<double> X;
			std::vector<double> Y;
			std::vector<double> E;
		
			PairedData(int n);	
			PairedData(std::vector<double> x, std::vector<double> y);
			PairedData(std::vector<double> x, std::vector<double> y,std::vector<double> errs);

			Pair operator [](int i) const;
			std::vector<Pair> GetPairs();
	};
}