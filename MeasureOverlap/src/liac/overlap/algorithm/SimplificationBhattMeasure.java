/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package liac.overlap.algorithm;

import liac.overlap.core.GaussianDistribution;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author liac01
 */
public class SimplificationBhattMeasure extends AlgorithmMeasure {

	public SimplificationBhattMeasure(GaussianDistribution first, GaussianDistribution second)
	{
		this.first = first;
		this.second = second;		
	}
	
	@Override
	public double overlapRate()
	{
		return this.bhattCoefficient(this.first.getMean(), this.second.getMean(), 
				                     this.first.getCov(), this.second.getCov(), 
									 this.first.getDet(), this.second.getDet());
	}
	
	private double bhattDistance(SimpleMatrix meanX, SimpleMatrix meanY, SimpleMatrix covX, SimpleMatrix covY, double detX, double detY)
	{
		SimpleMatrix distance = meanX.minus(meanY);
		
		SimpleMatrix covariance = covX.plus(covY).divide(2);

		SimpleMatrix bhattDistance = ((distance.transpose()).mult(covariance.invert()).mult(distance)).divide(8);

		return bhattDistance.get(0, 0) + Math.log1p(covariance.determinant() / ((Math.pow(detX * detY, 0.5)))) / 2;
	}

	private double bhattCoefficient(SimpleMatrix meanX, SimpleMatrix meanY, SimpleMatrix covX, SimpleMatrix covY, double detX, double detY)
	{
		return 1 / Math.exp(this.bhattDistance(meanX, meanY, covX, covY, detX, detY));
	}
}
