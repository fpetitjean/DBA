/*******************************************************************************
 * Copyright (C) 2018 Francois PETITJEAN 
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/ 

import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;


/**
 * This toy class show the use of DBA.
 * @author Francois Petitjean
 */
public class DBA{
	static final long serialVersionUID = 1L;

	private final static int NIL = -1;
	private final static int DIAGONAL = 0;
	private final static int LEFT = 1;
	private final static int UP = 2;

	/**
	 * Performs the DBA averaging by first finding the median over a sample, then doing n iterations of the update
	 * @param sequences set of sequences to average
	 */
	public static double[] performDBA(double[][] sequences) {
		
		int maxLength=0;
		for (int i = 0; i < sequences.length; i++) {
			maxLength = Math.max(maxLength,sequences[i].length);
		}
		double[][]costMatrix = new double[maxLength][maxLength];
		int[][]pathMatrix = new int[maxLength][maxLength];
		int medoidIndex = approximateMedoidIndex(sequences,costMatrix);
		double[]center = Arrays.copyOf(sequences[medoidIndex], sequences[medoidIndex].length);
		
		for (int i = 0; i < 15; i++) {
			center = DBAUpdate(center, sequences,costMatrix,pathMatrix);
		}
		return center;
	}
	
	private static int approximateMedoidIndex(double[][] sequences,double[][]mat) {
		/*
		 * we are finding the medoid, as this can take a bit of time, 
		 * if there is more than 50 time series, we sample 50 as possible 
		 * medoid candidates
		 */
		ArrayList<Integer>allIndices = new ArrayList<>();
		for (int i = 0; i < sequences.length; i++) {
			allIndices.add(i);
		}
		Collections.shuffle(allIndices);
		ArrayList<Integer>medianIndices = new ArrayList<>();
		for (int i = 0; i < sequences.length && i<50; i++) {
			medianIndices.add(allIndices.get(i));
		}
		
		int indexMedoid = -1;
		double lowestSoS = Double.MAX_VALUE;
		
		for (int medianCandidateIndex:medianIndices) {
			double[] possibleMedoid = sequences[medianCandidateIndex];
			double tmpSoS = sumOfSquares(possibleMedoid, sequences,mat);
			if (tmpSoS < lowestSoS) {
				indexMedoid = medianCandidateIndex;
				lowestSoS = tmpSoS;
			}
		}
		return indexMedoid;
	}
	
	private static double sumOfSquares(double[]sequence,double[][]sequences,double[][]mat) {
		double sos = 0.0;
		for (int i = 0; i < sequences.length; i++) {
			double dist = DTW(sequence,sequences[i],mat);
			sos += dist*dist;
		}
		return sos;
	}
	
	public static double DTW(double[]S,double []T,double[][]costMatrix) {
		int i, j;
		costMatrix[0][0] = squaredDistance(S[0],T[0]);
		for (i = 1; i < S.length; i++) {
			costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(S[i], T[0]);
		}
		for (j = 1; j < T.length; j++) {
			costMatrix[0][j] = costMatrix[0][j - 1]	+ squaredDistance(S[0],T[j]);
		}
		for (i = 1; i < S.length; i++) {
			for (j = 1; j < T.length; j++) {
				costMatrix[i][j] = Min3(costMatrix[i - 1][j - 1],
						costMatrix[i][j - 1], costMatrix[i - 1][j])
						+ squaredDistance(S[i],T[j]);
			}
		}
		
		return sqrt(costMatrix[S.length - 1][T.length - 1]);
	}
	
	private static double[] DBAUpdate(double[] C, double[][] sequences,double[][]costMatrix,int[][]pathMatrix) {
		double[]updatedMean = new double[C.length];
		int[]nElementsForMean= new int[C.length];
		
		int i, j, indiceRes;
		double res = 0.0;
		int centerLength = C.length;
		int seqLength;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = squaredDistance(C[0], T[0]);
			pathMatrix[0][0] = DBA.NIL;

			for (i = 1; i < centerLength; i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(C[i], T[0]);
				pathMatrix[i][0] = DBA.UP;
			}
			for (j = 1; j < seqLength; j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + squaredDistance(T[j], C[0]);
				pathMatrix[0][j] = DBA.LEFT;
			}

			for (i = 1; i < centerLength; i++) {
				for (j = 1; j < seqLength; j++) {
					indiceRes = DBA.ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							break;
						case UP:
							res = costMatrix[i - 1][j];
							break;
					}
					costMatrix[i][j] = res + squaredDistance(C[i], T[j]);
				}
			}

			i = centerLength - 1;
			j = seqLength - 1;

			while(pathMatrix[i][j]!=DBA.NIL) {
				updatedMean[i]+=T[j];
				nElementsForMean[i]++;
				switch (pathMatrix[i][j]) {
					case DIAGONAL:
						i = i - 1;
						j = j - 1;
						break;
					case LEFT:
						j = j - 1;
						break;
					case UP:
						i = i - 1;
						break;
				}
			}
			assert(i!=0 || j!=0);
			updatedMean[i]+=T[j];
			nElementsForMean[i]++;
		}

		for (int t = 0; t < centerLength; t++) {
			updatedMean[t] /= nElementsForMean[t];
		}
		
		return updatedMean;

	}



	private static double Min3(final double a, final double b, final double c) {
		if (a < b) {
			if (a < c) {
				return a;
			} else {
				return c;
			}
		} else {
			if (b < c) {
				return b;
			} else {
				return c;
			}
		}
	}

	private static int ArgMin3(final double a, final double b, final double c) {
		if (a < b) {
			if (a < c) {
				return 0;
			} else {
				return 2;
			}
		} else {
			if (b < c) {
				return 1;
			} else {
				return 2;
			}
		}
	}

	private static double squaredDistance(double a, double b) {
		double diff = a-b;
		return diff*diff;
	}

	public static void main(String [] args){
		Random r = new Random(3071980);
		double [][]sequences = new double[50][];
		for(int i=0;i<sequences.length;i++){
			sequences[i] = new double[20];
			for(int j=0;j<sequences[i].length;j++){
				sequences[i][j] = Math.cos(r.nextDouble()*j/20.0*Math.PI) ;
			}
		}
		double [] averageSequence = performDBA(sequences);
		
		System.out.print("[");
		for(int j=0;j<averageSequence.length;j++){
			System.out.print(averageSequence[j]+" ");
		}
		System.out.println("]");
	}
}
