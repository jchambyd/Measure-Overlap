/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package liac.overlap.sample;

import liac.overlap.algorithm.AlgorithmMeasure;
import liac.overlap.algorithm.BhattacharyyaMeasure;
import liac.overlap.algorithm.OLRMeasure;
import liac.overlap.core.GaussianDistribution;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author liac01
 */
public class Main {	

	public static void main(String [] args)
	{
		SimpleMatrix mu1 = new SimpleMatrix(2, 1);
		SimpleMatrix mu2 = new SimpleMatrix(2, 1);

		mu2.set(0, 0, 2);

		SimpleMatrix sigma1 = new SimpleMatrix(2, 2);
		SimpleMatrix sigma2 = new SimpleMatrix(2, 2);

		sigma1.set(0, 0, 1);
		sigma1.set(1, 1, 1);

		sigma2.set(0, 0, 2.17);
		sigma2.set(0, 1, 1.82);
		sigma2.set(1, 0, 1.82);
		sigma2.set(1, 1, 2.17);
		
		GaussianDistribution first = new GaussianDistribution(sigma1, mu1, 0.5);
		GaussianDistribution second = new GaussianDistribution(sigma2, mu2, 0.5);
		
		AlgorithmMeasure loAlgorithm = new OLRMeasure(first, second);
		AlgorithmMeasure loAlgorithm1 = new BhattacharyyaMeasure(first, second);
		
		System.out.println("The OLR of this mixture Gaussian: " + loAlgorithm.overlapRate());
		System.out.println("The Bhatt of this mixture Gaussian: " + loAlgorithm1.overlapRate());
	}
}
