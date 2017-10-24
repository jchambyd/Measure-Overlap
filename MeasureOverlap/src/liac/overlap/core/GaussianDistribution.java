/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package liac.overlap.core;

import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author liac01
 */
public class GaussianDistribution {

	private SimpleMatrix cov;
	private SimpleMatrix invCov;
	private SimpleMatrix mean;
	private double prior;
	private double det;

	public GaussianDistribution(SimpleMatrix cov, SimpleMatrix mean, double prior)
	{
		this.cov = cov;
		this.invCov = this.cov.invert();
		this.prior = prior;
		this.mean = mean;
		this.det = this.cov.determinant();
	}

	/**
	 * Calculate multivariate probability density function
	 *
	 * @param x input vector
	 * @return a densidade de probabilidade
	 */
	public double mvnpdf(SimpleMatrix x)
	{
		double dim = x.getNumElements();
		SimpleMatrix distance = x.minus(this.getMean());

		double pdf = Math.exp(-0.5 * distance.transpose().dot(this.invCov.mult(distance)))
				/ (Math.pow(2 * Math.PI, dim / 2.0) * Math.sqrt(this.det));

		pdf = Double.isNaN(pdf) ? 0 : pdf;
		pdf = Double.isInfinite(pdf) ? Double.MAX_VALUE : pdf;

		return pdf;
	}

	/**
	 * @return the cov
	 */
	public SimpleMatrix getCov()
	{
		return cov;
	}

	/**
	 * @param cov the cov to set
	 */
	public void setCov(SimpleMatrix cov)
	{
		this.cov = cov;
	}

	/**
	 * @return the invCov
	 */
	public SimpleMatrix getInvCov()
	{
		return invCov;
	}

	/**
	 * @param invCov the invCov to set
	 */
	public void setInvCov(SimpleMatrix invCov)
	{
		this.invCov = invCov;
	}

	/**
	 * @return the mean
	 */
	public SimpleMatrix getMean()
	{
		return mean;
	}

	/**
	 * @param mean the mean to set
	 */
	public void setMean(SimpleMatrix mean)
	{
		this.mean = mean;
	}

	/**
	 * @return the prior
	 */
	public double getPrior()
	{
		return prior;
	}

	/**
	 * @param prior the prior to set
	 */
	public void setPrior(double prior)
	{
		this.prior = prior;
	}

	/**
	 * @return the det
	 */
	public double getDet()
	{
		return det;
	}

	/**
	 * @param det the det to set
	 */
	public void setDet(double det)
	{
		this.det = det;
	}
}
