#include "generator.h"

namespace cyclups
{


	namespace generator
	{
		std::default_random_engine randomiser;
		std::normal_distribution<double> gaussian (0.0,1.);
		std::uniform_real_distribution<double> uniform(0,1);
		PairedData generateSamples(std::vector<double> x, functionPointer f, double noise,double heteroskedacity)
		{
			int n = x.size();
			std::vector<double> y(n);
			std::vector<double> e(n);
			double minNoise = std::max(1e-5,noise * (1.0 - heteroskedacity));
			double maxNoise = noise* (1.0 + heteroskedacity);
			for (int i = 0; i < n; ++i)
			{
				double r = uniform(randomiser);
				double rx = (x[i] - x[0])/(x[x.size()-1] - x[0]);
				double snoise = r*rx*(maxNoise - minNoise) + minNoise;
				std::cout << snoise << std::endl;
				y[i] = f(x[i]) + gaussian(randomiser) *snoise;
				e[i] = snoise;
			}
			return PairedData(x,y,e);
		}


		PairedData Sample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity)
		{
			return UniformXSample(N,f,xMin, xMax, noise,heteroskedacity);
		}

		PairedData UniformXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity)
		{
			std::vector<double> x = JSL::Vector::linspace(xMin,xMax,N);

			return generateSamples(x,f,noise,heteroskedacity);
		}

		PairedData RandomXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity)
		{
			std::vector<double> x = JSL::Vector::RandVec(N,xMin,xMax);
			std::sort(x.begin(),x.end());
			return generateSamples(x,f,noise,heteroskedacity);
		}
		PairedData NoisyXSample(int N, functionPointer f, double xMin, double xMax, double noise,double heteroskedacity)
		{
			std::vector<double> x = JSL::Vector::linspace(xMin,xMax,N);
			double dx = (x[1] - x[0]);
			for (int i = 0; i < N; ++i)
			{
				x[i] += gaussian(randomiser) * dx/2;
			}
			std::sort(x.begin(),x.end());
			return generateSamples(x,f,noise,heteroskedacity);
		}
	}
}