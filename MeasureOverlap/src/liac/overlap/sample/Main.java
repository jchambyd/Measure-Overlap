/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package liac.overlap.sample;

import java.util.ArrayList;
import liac.overlap.algorithm.AlgorithmMeasure;
import liac.overlap.algorithm.BhattacharyyaMeasure;
import liac.overlap.algorithm.OLRMeasure;
import liac.overlap.core.GaussianDistribution;
import liac.overlap.loader.DataLoader;
import liac.overlap.loader.DataLoaderException;
import liac.overlap.loader.Dataset;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author liac01
 */
public class Main {

	public static void main(String[] args) throws DataLoaderException
	{
		Dataset dataset = DataLoader.loadARFF("data/iris.arff");

		ArrayList<SimpleMatrix> means = Main.getMeans(dataset);
		ArrayList<SimpleMatrix> covs = Main.getCovarianceMatrix(dataset, means);

		Main.getOverlapRate(covs, means, 1, 2, dataset);
	}

	public static void mxSimpleTest()
	{
		SimpleMatrix mu1 = new SimpleMatrix(2, 1);
		SimpleMatrix mu2 = new SimpleMatrix(2, 1);

		mu2.set(0, 0, 3);

		SimpleMatrix sigma1 = new SimpleMatrix(2, 2);
		SimpleMatrix sigma2 = new SimpleMatrix(2, 2);

		sigma1.set(0, 0, 1);
		sigma1.set(1, 1, 1);

		sigma2.set(0, 0, 2.17);
		sigma2.set(0, 1, 1.82);
		sigma2.set(1, 0, 1.82);
		sigma2.set(1, 1, 2.17);

		GaussianDistribution first = new GaussianDistribution(sigma1, mu1, 0.3);
		GaussianDistribution second = new GaussianDistribution(sigma2, mu2, 0.7);

		AlgorithmMeasure loAlgorithm = new OLRMeasure(first, second);
		AlgorithmMeasure loAlgorithm1 = new BhattacharyyaMeasure(first, second);

		System.out.println("The OLR of this mixture Gaussian: " + loAlgorithm.overlapRate());
		System.out.println("The Bhatt of this mixture Gaussian: " + loAlgorithm1.overlapRate());
	}

	public static double getOverlapRate(ArrayList<SimpleMatrix> covs, ArrayList<SimpleMatrix> means, int classOne, int classTwo, int featOne, int featTwo)
	{
		SimpleMatrix mu1 = new SimpleMatrix(2, 1);
		SimpleMatrix mu2 = new SimpleMatrix(2, 1);

		mu1.set(0, 0, means.get(classOne).get(featOne, 0));
		mu1.set(1, 0, means.get(classOne).get(featTwo, 0));

		mu2.set(0, 0, means.get(classTwo).get(featOne, 0));
		mu2.set(1, 0, means.get(classTwo).get(featTwo, 0));

		SimpleMatrix sigma1 = new SimpleMatrix(2, 2);
		SimpleMatrix sigma2 = new SimpleMatrix(2, 2);

		sigma1.set(0, 0, covs.get(classOne).get(featOne, featOne));
		sigma1.set(0, 1, covs.get(classOne).get(featOne, featTwo));
		sigma1.set(1, 0, covs.get(classOne).get(featTwo, featOne));
		sigma1.set(1, 1, covs.get(classOne).get(featTwo, featTwo));

		sigma2.set(0, 0, covs.get(classTwo).get(featOne, featOne));
		sigma2.set(0, 1, covs.get(classTwo).get(featOne, featTwo));
		sigma2.set(1, 0, covs.get(classTwo).get(featTwo, featOne));
		sigma2.set(1, 1, covs.get(classTwo).get(featTwo, featTwo));

		GaussianDistribution first = new GaussianDistribution(sigma1, mu1, 0.5);
		GaussianDistribution second = new GaussianDistribution(sigma2, mu2, 0.5);

		AlgorithmMeasure loAlgorithm = new OLRMeasure(first, second);
		AlgorithmMeasure loAlgorithm1 = new BhattacharyyaMeasure(first, second);

		System.out.println("The OLR of this mixture Gaussian: " + loAlgorithm.overlapRate());
		System.out.println("The Bhatt of this mixture Gaussian: " + loAlgorithm1.overlapRate());

		return 0;
	}

	public static double getOverlapRate(ArrayList<SimpleMatrix> covs, ArrayList<SimpleMatrix> means, int classOne, int classTwo, Dataset dataset)
	{
		int dimension = dataset.getInputSize();

		SimpleMatrix mu1 = means.get(classOne).extractMatrix(0, dimension, 0, 1);
		SimpleMatrix mu2 = means.get(classTwo).extractMatrix(0, dimension, 0, 1);

		SimpleMatrix sigma1 = covs.get(classOne).extractMatrix(0, dimension, 0, dimension);
		SimpleMatrix sigma2 = covs.get(classTwo).extractMatrix(0, dimension, 0, dimension);

		GaussianDistribution first = new GaussianDistribution(sigma1, mu1, 0.5);
		GaussianDistribution second = new GaussianDistribution(sigma2, mu2, 0.5);

		AlgorithmMeasure loAlgorithm = new OLRMeasure(first, second);
		AlgorithmMeasure loAlgorithm1 = new BhattacharyyaMeasure(first, second);

		System.out.println("The OLR of this mixture Gaussian: " + loAlgorithm.overlapRate());
		System.out.println("The Bhatt of this mixture Gaussian: " + loAlgorithm1.overlapRate());

		return 0;
	}

	public static ArrayList<SimpleMatrix> getMeans(Dataset dataset)
	{
		ArrayList<SimpleMatrix> means = new ArrayList<>();
		SimpleMatrix ds = new SimpleMatrix(dataset.getDataset());
		SimpleMatrix mean = new SimpleMatrix(ds.numRows(), 1);
		int numClasses = dataset.getNumClasses(), numAttr = dataset.getInputSize();
		int nums[] = new int[numClasses];

		for (int i = 0; i < numClasses; i++) {
			means.add(mean);
		}

		//Calculating mean vectors
		for (int i = 0; i < ds.numCols(); i++) {
			for (int j = 0; j < numClasses; j++) {
				if (ds.get(numAttr + j, i) == 1) {
					means.set(j, means.get(j).plus(ds.extractVector(false, i)));
					++nums[j];
					break;
				}
			}
		}
		for (int i = 0; i < numClasses; i++) {
			means.set(i, means.get(i).scale(1.0 / nums[i]));
		}

		return means;
	}

	public static ArrayList<SimpleMatrix> getCovarianceMatrix(Dataset dataset, ArrayList<SimpleMatrix> means)
	{
		ArrayList<SimpleMatrix> covs = new ArrayList<>();
		SimpleMatrix ds = new SimpleMatrix(dataset.getDataset());
		int numClasses = dataset.getNumClasses(), numAttr = dataset.getInputSize();
		int nums[] = new int[numClasses];

		for (int i = 0; i < numClasses; i++) {
			covs.add(new SimpleMatrix(ds.numRows(), ds.numRows()));
		}

		//Calculating mean vectors
		for (int i = 0; i < ds.numCols(); i++) {
			for (int j = 0; j < numClasses; j++) {
				if (ds.get(numAttr + j, i) == 1) {
					SimpleMatrix tmp = new SimpleMatrix(ds.extractVector(false, i));
					tmp = tmp.minus(means.get(j));
					covs.set(j, covs.get(j).plus(tmp.mult(tmp.transpose())));
					++nums[j];
					break;
				}
			}
		}
		for (int i = 0; i < numClasses; i++) {
			covs.set(i, covs.get(i).scale(1.0 / (nums[i] - 1)));
		}

		return covs;
	}

	public static SimpleMatrix getCovarianceMatrix(SimpleMatrix ds)
	{
		SimpleMatrix result = new SimpleMatrix(ds.numRows(), ds.numRows());
		SimpleMatrix mean = new SimpleMatrix(ds.numRows(), 1);

		//Calculating the mean vector
		for (int i = 0; i < ds.numCols(); i++) {
			mean = mean.plus(ds.extractVector(false, i));
		}
		mean = mean.scale(1.0 / ds.numCols());

		//Calculating covariance matrix
		for (int i = 0; i < ds.numCols(); i++) {
			SimpleMatrix tmp = new SimpleMatrix(ds.extractVector(false, i));
			tmp = tmp.minus(mean);
			result = result.plus(tmp.mult(tmp.transpose()));
		}
		result = result.scale(1.0 / (ds.numCols() + 1));

		return result;
	}
}
