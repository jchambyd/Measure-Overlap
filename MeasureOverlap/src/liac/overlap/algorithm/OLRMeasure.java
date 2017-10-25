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
	
	private double RC(SimpleMatrix x)
	{
		double E = this.first.getInvCov().get(0,0) * (x.get(0) - this.first.getMean().get(0));
		double F = this.first.getInvCov().get(1,0) * (x.get(0) - this.first.getMean().get(0));
		double G = this.second.getInvCov().get(0,0) * (x.get(0) - this.second.getMean().get(0));
		double H = this.second.getInvCov().get(1,0) * (x.get(0) - this.second.getMean().get(0));
		
        double I = E * this.second.getInvCov().get(1, 1) - F * this.second.getInvCov().get(0, 1);
		double J = H * this.first.getInvCov().get(0, 1) - G * this.first.getInvCov().get(1, 1);
		double K = this.first.getInvCov().get(0, 1) * this.second.getInvCov().get(1, 1) - 
				   this.second.getInvCov().get(0, 1) * this.first.getInvCov().get(1, 1);
        
        double M = E * H - F * G;
		double a = K;
        double b = I + J - K * (this.second.getMean().get(1) + this.first.getMean().get(1));
        double c = M - I * this.second.getMean().get(1) - J * this.first.getMean().get(1) + K * (this.second.getMean().get(1) * this.first.getMean().get(1));

        if ((Math.pow(b, 2) - 4 * a * c) < 0)
            return Double.NaN;

        return Math.max((-b + Math.sqrt(Math.pow(b, 2) - 4 * a * c)) / (2 * a), (-b - Math.sqrt(Math.pow(b, 2) - 4*a*c)) / (2*a));
	}
	
	
	@Override
	public double overlapRate()
	{
		double e = Math.sqrt(Math.pow(this.first.getMean().get(0) - this.second.getMean().get(0), 2) + Math.pow(this.first.getMean().get(1) - this.second.getMean().get(1), 2)) / this.steps;
		
		double x_step = e * (this.first.getMean().get(0) - this.second.getMean().get(0)); // Each step for x
        double y_step = e * (this.first.getMean().get(1) - this.second.getMean().get(1)); // Each step for y
		
		double p_x = this.first.getMean().get(0), p_y = this.first.getMean().get(1);
		
		SimpleMatrix p = new SimpleMatrix(2, 1);
		p.set(0, 0, p_x);
		
		do
		{
			p_x = p_x - x_step;
			p.set(0, p_x);
			p_y = p_y - y_step; //this.RC(p);
		} while (Double.isNaN(p_y));
		
		p.set(1, 0, p_y);
		
		SimpleMatrix p_pre = new SimpleMatrix(this.first.getMean());
		SimpleMatrix p_next = new SimpleMatrix(2, 1);
		
		double pdf_p, pdf_p_pre, pdf_p_next;
		double p_sub_max = Math.min(this.getPDF(this.first.getMean()), this.getPDF(this.second.getMean()));
		double p_min = p_sub_max;
		int index = 0;
		
		while (index < this.steps)
		{	
			// Next point on ridge curve
			p_next.set(0, 0, p.get(0) - x_step);
			p_y = p.get(1, 0) - y_step; //this.RC(p_next);
						
            if (!Double.isNaN(p_y))
			{
				p_next.set(1, 0, p_y);
				// Calculate pdf values
				pdf_p = this.getPDF(p); 
				pdf_p_pre = this.getPDF(p_pre);
				pdf_p_next = this.getPDF(p_next);
				
                if (pdf_p < pdf_p_pre && pdf_p < pdf_p_next)
                    p_min = pdf_p;
			}
            p_pre = new SimpleMatrix(p);
            p = new SimpleMatrix(p_next);
			index += 1;
		}
		
		return (p_min < p_sub_max) ? p_min / p_sub_max : 1;
	}
}
