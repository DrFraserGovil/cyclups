#include "OptimiserProperties.h"

void cyclups::OptimiserProperties::Clear()
{
	PrevScore = 0;
	GradientMemory= 0;
	ScoreMemory = 0;
	TrueAlpha = alpha;
	MaxAlpha = alpha*AlphaFactor;
	MinAlpha = alpha/AlphaFactor;
	TriggeringStep = 0;
	NegativeCounter = 0;
	ReachedMaxSteps = false;
	GradientConverged = false;
	ScoreConverged = false;
	Converged = false;
}
void cyclups::OptimiserProperties::CheckConvergence(int l, double gradnorm)
{
				
	if (l >= MaxSteps)
	{
		ReachedMaxSteps = true;
	}
	
	double earlyCorrector = 1.0/(1.0 - pow(ConvergenceMemory,l+1));
	
	GradientMemory =  (ConvergenceMemory * GradientMemory + (1.0 - ConvergenceMemory) * gradnorm);
	if (GradientMemory * earlyCorrector < ConvergedGradient)
	{
		triggeringGradient = gradnorm;
		GradientConverged = true;
	}
	
	Converged = (ReachedMaxSteps || GradientConverged || ScoreConverged) && l > MinSteps;
	if (Converged)
	{
		TriggeringStep = l;
		alpha = TrueAlpha; //ensures it always goes back into the correct state
	}		
	// if (l%10 == 0)
	// {
	// std::cout << l << "  " << GradientMemory*earlyCorrector << "/" << ConvergedGradient << "   " << ScoreMemory*earlyCorrector << "/" << ConvergedScore << std::endl;}
}

void cyclups::OptimiserProperties::CheckConvergence(int l, double gradnorm, double score)
{
	//do score bit
	
	if (l > 0)
	{
		double alphaCorrector = (alpha/TrueAlpha);
		double earlyCorrector = 1.0/(1.0 - pow(ConvergenceMemory,l+1));
		double scoreDelta = abs((score - PrevScore)/PrevScore);
		ScoreMemory = (ConvergenceMemory * ScoreMemory + (1.0 - ConvergenceMemory) * alphaCorrector *scoreDelta);
		if (ScoreMemory* earlyCorrector < ConvergedScore)
		{
			triggeringScore = scoreDelta;
			ScoreConverged = true;
		}
	}
	UpdateAlpha(score);
	PrevScore = score;
	CheckConvergence(l,gradnorm);
}

void cyclups::OptimiserProperties::PrintReason()
{
	std::cout << "The Optimiser halted at step " << TriggeringStep << " because:\n";
	if (ReachedMaxSteps)
	{
		std::cout << "\t-Reached max iteration count.\n";
	}
	if (GradientConverged)
	{
		std::cout << "\t-Mean-Gradient converged below " << ConvergedGradient << "(" << triggeringGradient << ")\n";
	}
	if (ScoreConverged)
	{
		std::cout << "\t-Mean-Score has not changed by more than " << 100*ConvergedScore << "% (" <<  triggeringScore << ")\n";
	}
}

void cyclups::OptimiserProperties::UpdateAlpha(double score)
{			
	if (score > PrevScore)
	{
		NegativeCounter +=2;
		if (NegativeCounter >= 10)
		{
			NegativeCounter = 0;
			alpha *= 0.7;
			if (alpha < MinAlpha)
			{
				alpha = MinAlpha;
				MinAlpha *= 0.99;
			}
			// alpha = std::max(MinAlpha,alpha *0.7);
		}
	}
	else
	{
		NegativeCounter -=1;
		if (NegativeCounter == -2)
		{
			NegativeCounter = 0;
			alpha *= 1.07;
			if (alpha > MaxAlpha)
			{
				alpha = MaxAlpha;
				MaxAlpha *=1.01;
			}
			// alpha = std::min(MaxAlpha,alpha * 1.03);
		}
	}

};
