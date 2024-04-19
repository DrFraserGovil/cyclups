#pragma once
#include "dataArrays.h"
#include "../customTypes.h"
#include "JSL.h"
#include <random>
namespace cyclups
{


	namespace generator
	{
		extern std::default_random_engine randomiser;

		PairedData Sample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity);
		
		PairedData UniformXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity);

		PairedData NoisyXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity);

		PairedData RandomXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity);

	}
}