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
public abstract class AlgorithmMeasure {

	protected GaussianDistribution first;
	protected GaussianDistribution second;

	public abstract double overlapRate();

	protected double getPDF(SimpleMatrix x)
	{
		return this.first.getPrior() * this.first.mvnpdf(x) + this.second.getPrior() * this.second.mvnpdf(x);
	}

	/**
	 * @return the first
	 */
	public GaussianDistribution getFirst()
	{
		return first;
	}

	/**
	 * @param first the first to set
	 */
	public void setFirst(GaussianDistribution first)
	{
		this.first = first;
	}

	/**
	 * @return the second
	 */
	public GaussianDistribution getSecond()
	{
		return second;
	}

	/**
	 * @param second the second to set
	 */
	public void setSecond(GaussianDistribution second)
	{
		this.second = second;
	}
}
