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
public class OLRMeasure extends AlgorithmMeasure {

	private double steps = 100;

	public OLRMeasure(GaussianDistribution first, GaussianDistribution second)
	{
		this.first = first;
		this.second = second;
	}

	// Just for two dimension, since (x1, ?), is calculated (x1, x2) and return x2 
	private double RC(SimpleMatrix x)
	{
		double E = this.first.getInvCov().get(0, 0) * (x.get(0) - this.first.getMean().get(0));
		double F = this.first.getInvCov().get(1, 0) * (x.get(0) - this.first.getMean().get(0));
		double G = this.second.getInvCov().get(0, 0) * (x.get(0) - this.second.getMean().get(0));
		double H = this.second.getInvCov().get(1, 0) * (x.get(0) - this.second.getMean().get(0));

		double I = E * this.second.getInvCov().get(1, 1) - F * this.second.getInvCov().get(0, 1);
		double J = H * this.first.getInvCov().get(0, 1) - G * this.first.getInvCov().get(1, 1);
		double K = this.first.getInvCov().get(0, 1) * this.second.getInvCov().get(1, 1)
				- this.second.getInvCov().get(0, 1) * this.first.getInvCov().get(1, 1);

		double M = E * H - F * G;
		double a = K;
		double b = I + J - K * (this.second.getMean().get(1) + this.first.getMean().get(1));
		double c = M - I * this.second.getMean().get(1) - J * this.first.getMean().get(1) + K * (this.second.getMean().get(1) * this.first.getMean().get(1));

		if ((Math.pow(b, 2) - 4 * a * c) < 0) {
			return Double.NaN;
		}

		return Math.max((-b + Math.sqrt(Math.pow(b, 2) - 4 * a * c)) / (2 * a), (-b - Math.sqrt(Math.pow(b, 2) - 4 * a * c)) / (2 * a));
	}

	@Override
	public double overlapRate()
	{
		int dimension = this.first.getMean().getNumElements();
		double delta = 0;
		SimpleMatrix step = new SimpleMatrix(dimension, 1);

		for (int i = 0; i < dimension; i++) {
			delta += Math.pow(this.first.getMean().get(i) - this.second.getMean().get(i), 2);
		}

		delta = Math.sqrt(delta) / this.steps;

		for (int i = 0; i < dimension; i++) {
			step.set(i, delta * (this.first.getMean().get(i) - this.second.getMean().get(i)));
		}

		SimpleMatrix p_pre = new SimpleMatrix(this.first.getMean());
		SimpleMatrix p = p_pre.minus(step);
		SimpleMatrix p_next = p.minus(step);

		double pdf_p, pdf_p_pre, pdf_p_next;
		// Small value of two peaks
		double p_sub_max = Math.min(this.getPDF(this.first.getMean()), this.getPDF(this.second.getMean()));
		double p_min = p_sub_max;
		int index = 0;

		while (index < this.steps) {
			// Calculate pdf values
			pdf_p_pre = this.getPDF(p_pre);
			pdf_p = this.getPDF(p);
			pdf_p_next = this.getPDF(p_next);

			if (pdf_p < pdf_p_pre && pdf_p < pdf_p_next) {
				p_min = pdf_p;
			}

			p_pre = new SimpleMatrix(p);
			p = new SimpleMatrix(p_next);
			// Next point on ridge curve
			p_next = p.minus(step);

			index += 1;
		}

		return (p_min < p_sub_max) ? p_min / p_sub_max : 1;
	}
}
