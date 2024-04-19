#include "dataArrays.h"

namespace cyclups
{
	

	PairedData::PairedData(int n)
	{
		X.resize(n,0);
		Y.resize(n,0);
		E.resize(n,0);
	}

	PairedData::PairedData(std::vector<double> x, std::vector<double> y)
	{
		if (x.size() != y.size())
		{
			JSL::Error("Paired Vectors must be the same size");
		}
		X = x;
		Y = y;
		E.resize(x.size(),0);
	}
	PairedData::PairedData(std::vector<double> x, std::vector<double> y, std::vector<double> errs)
	{
		if (x.size() != y.size() || y.size() != errs.size())
		{
			JSL::Error("Paired Vectors must be the same size");
		}
		X = x;
		Y = y;
		E = errs;
	}

	Pair PairedData::operator[](int i) const
	{
		return Pair(X[i],Y[i],E[i]);
	}

	std::vector<Pair> PairedData::GetPairs()
	{
		std::vector<Pair> out;
		for (int i = 0; i < X.size(); ++i)
		{
			out.push_back(Pair(X[i],Y[i],E[i]));
		}
	}
}
