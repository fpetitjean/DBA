/*******************************************************************************
 * Copyright (C) 2016 Chang Wei Tan 
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
 * This toy class show the use of DBA with warping window.	
 * @author Chang Wei Tan
 */
public class DBAWarpingWindow {
	static final long serialVersionUID = 1L;

	private final static int NIL = -1;
	private final static int DIAGONAL = 0;
	private final static int LEFT = 1;
	private final static int UP = 2;

	/**
	 * This attribute is used in order to initialize only once the matrixes
	 */
	private final static int MAX_SEQ_LENGTH = 2000;

	/**
	 * store the cost of the alignment
	 */
	private static double[][] costMatrix = new double[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	
	/**
	 * store the warping path
	 */
	private static int[][] pathMatrix = new int[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];

	/**
	 * Performs the DBA averaging by first finding the median over a sample, then doing n iterations of the update
	 * @param sequences set of sequences to average
	 * @param w warping window
	 */
	public static double[] performDBA(double[][] sequences,int w) {
		int medoidIndex = approximateMedoidIndex(sequences,w);
		double[]center = Arrays.copyOf(sequences[medoidIndex], sequences[medoidIndex].length);
		
		for (int i = 0; i < 15; i++) {
			center = DBAUpdate(center, sequences,w);
		}
		return center;
	}
	
	public static int approximateMedoidIndex(double[][] sequences,int w) {
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
			double tmpSoS = sumOfSquares(possibleMedoid, sequences,w);
			if (tmpSoS < lowestSoS) {
				indexMedoid = medianCandidateIndex;
				lowestSoS = tmpSoS;
			}
		}
		return indexMedoid;
	}
	
	private static double sumOfSquares(double[]sequence,double[][]sequences, int w) {
		double sos = 0.0;
		for (int i = 0; i < sequences.length; i++) {
			double dist = DTW(sequence,sequences[i],w);
			sos += dist*dist;
		}
		return sos;
	}
	
	public static double DTW(double[]S,double []T, int w) {
		int i, j;
		costMatrix[0][0] = squaredDistance(S[0],T[0]);
		for (i = 1; i < Math.min(S.length,w+1); i++) {
			costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(S[i], T[0]);
		}
		for (j = 1; j < Math.min(T.length,w+1); j++) {
			costMatrix[0][j] = costMatrix[0][j - 1]	+ squaredDistance(S[0],T[j]);
		}
		if (j < T.length)
			costMatrix[0][j] = Double.POSITIVE_INFINITY;
		
		for (i = 1; i < S.length; i++) {
			int jStart = Math.max(1, i-w);
			int jStop = Math.min(T.length, i+w+1);
			int indexInftyLeft = i-w-1;
			if(indexInftyLeft>=0)costMatrix[i][indexInftyLeft] = Double.POSITIVE_INFINITY;
			for (j = jStart; j < jStop; j++) {
				costMatrix[i][j] = Min3(costMatrix[i - 1][j - 1],
						costMatrix[i][j - 1], costMatrix[i - 1][j])
						+ squaredDistance(S[i],T[j]);
			}
			if (jStop < T.length)
				costMatrix[i][jStop] = Double.POSITIVE_INFINITY;
		}
		
		return sqrt(costMatrix[S.length - 1][T.length - 1]);
	}
	

	/**
	 * Dtw Barycenter Averaging (DBA)
	 * @param C average sequence to update
	 * @param sequences set of sequences to average
	 * @param w warping window size for DTW
	 */
	public static synchronized double[] DBAUpdate(double[] C, double[][] sequences, int w) {
		double[]updatedMean = new double[C.length];
		int[]nElementsForMean= new int[C.length];
		
		int i, j, indiceRes;
		double res = 0.0;
		int centerLength = C.length;
		int seqLength, jStart, jStop;;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = squaredDistance(C[0], T[0]);
			pathMatrix[0][0] = NIL;

			for (i = 1; i < Math.min(centerLength, 1+w); i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(C[i], T[0]);
				pathMatrix[i][0] = UP;
			}
			for (j = 1; j < Math.min(seqLength, 1+w); j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + squaredDistance(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
			}
			if (j < seqLength)
				costMatrix[0][j] = Double.POSITIVE_INFINITY;

			for (i = 1; i < centerLength; i++) {
				jStart = Math.max(1, i-w);
				jStop = Math.min(seqLength, i+w+1);
				int indexInftyLeft = i-w-1;
				if(indexInftyLeft>=0)costMatrix[i][indexInftyLeft] = Double.POSITIVE_INFINITY;
				for (j = jStart; j < jStop; j++) {
					indiceRes = ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
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
				if (jStop < seqLength)
					costMatrix[i][jStop] = Double.POSITIVE_INFINITY;
			}

			i = centerLength - 1;
			j = seqLength - 1;
			
			while(pathMatrix[i][j]!=NIL) {
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
			if(i!=0 ||j!=0)throw new RuntimeException("This should never occur, please contact Francois <francois.petitjean@monash.edu>");
			//i and j should be 0
			updatedMean[i]+=T[j];
			nElementsForMean[i]++;
		}
		
		for (int t = 0; t < centerLength; t++) {
			updatedMean[t] /= nElementsForMean[t];
		}
		
		return updatedMean;
	}


	public static double Min3(final double a, final double b, final double c) {
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

	public static int ArgMin3(final double a, final double b, final double c) {
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

	public static double squaredDistance(double a, double b) {
		return (a - b) * (a - b);
	}


	public static double barycenter(final Object... tab) {
		if (tab.length < 1) {
			throw new RuntimeException("empty double tab");
		}
		double sum = 0.0;
		sum = 0.0;
		for (Object o : tab) {
			sum += ((Double) o);
		}
		return sum / tab.length;
	}

	public static void main(String [] args){
		int w = 5;
		Random r = new Random(3071980);
		double [][]sequences = new double[50][];
		for(int i=0;i<sequences.length;i++){
			sequences[i] = new double[20];
			for(int j=0;j<sequences[i].length;j++){
				sequences[i][j] = Math.cos(r.nextDouble()*j/20.0*Math.PI) ;
			}
		}
		double [] averageSequence = performDBA(sequences,w);
		
		System.out.print("[");
		for(int j=0;j<averageSequence.length;j++){
			System.out.print(averageSequence[j]+" ");
		}
		System.out.println("]");
	}
}
