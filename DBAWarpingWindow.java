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

import java.util.ArrayList;

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
	private final static int MAX_SEQ_LENGTH = 20;

	/**
	 * store the cost of the alignment
	 */
	private static double[][] costMatrix = new double[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];
	
	/**
	 * store the warping path
	 */
	private static int[][] pathMatrix = new int[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];

	/**
	 * Store the length of the optimal path in each cell
	 */
	private static int[][] optimalPathLength = new int[MAX_SEQ_LENGTH][MAX_SEQ_LENGTH];

	/**
	 * Dtw Barycenter Averaging (DBA)
	 * @param C average sequence to update
	 * @param sequences set of sequences to average
	 * @param w warping window size for DTW
	 */
	public static double[] DBA_update(double[] C, double[][] sequences, int w) {
		final ArrayList<Double>[] tupleAssociation = new ArrayList[C.length];
		for (int i = 0; i < tupleAssociation.length; i++) {
			tupleAssociation[i] = new ArrayList<Double>(sequences.length);
		}
		int nbTuplesAverageSeq, i, j, indiceRes;
		double res = 0.0;
		int centerLength = C.length;
		int seqLength, jStart, jStop;;

		for (double[] T : sequences) {
			seqLength = T.length;

			costMatrix[0][0] = distanceTo(C[0], T[0]);
			pathMatrix[0][0] = NIL;
			optimalPathLength[0][0] = 0;

			for (i = 1; i < Math.min(centerLength, 1+w); i++) {
				costMatrix[i][0] = costMatrix[i - 1][0] + distanceTo(C[i], T[0]);
				pathMatrix[i][0] = UP;
				optimalPathLength[i][0] = i;
			}
			for (j = 1; j < Math.min(seqLength, 1+w); j++) {
				costMatrix[0][j] = costMatrix[0][j - 1] + distanceTo(T[j], C[0]);
				pathMatrix[0][j] = LEFT;
				optimalPathLength[0][j] = j;
			}
			if (j < seqLength)
				costMatrix[0][j] = Double.POSITIVE_INFINITY;

			for (i = 1; i < centerLength; i++) {
				jStart = Math.max(1, i-w);
				jStop = Math.min(seqLength, i+w+1);
				int indexInftyLeft = i-w-1;
				if(indexInftyLeft>=0)warpingMatrix[i][indexInftyLeft] = Double.POSITIVE_INFINITY;
				for (j = jStart; j < jStop; j++) {
					indiceRes = ArgMin3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j]);
					pathMatrix[i][j] = indiceRes;
					switch (indiceRes) {
						case DIAGONAL:
							res = costMatrix[i - 1][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j - 1] + 1;
							break;
						case LEFT:
							res = costMatrix[i][j - 1];
							optimalPathLength[i][j] = optimalPathLength[i][j - 1] + 1;
							break;
						case UP:
							res = costMatrix[i - 1][j];
							optimalPathLength[i][j] = optimalPathLength[i - 1][j] + 1;
							break;
					}
					costMatrix[i][j] = res + distanceTo(C[i], T[j]);
				}
				if (jStop < seqLength)
					costMatrix[i][jStop] = Double.POSITIVE_INFINITY;
			}

			nbTuplesAverageSeq = optimalPathLength[centerLength-1][seqLength-1] + 1;

			i = centerLength - 1;
			j = seqLength - 1;
			
			for (int t = nbTuplesAverageSeq - 1; t >= 0; t--) {
				tupleAssociation[i].add(T[j]);
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
		}
		
		for (int t = 0; t < centerLength; t++) {
			C[t] = barycenter((tupleAssociation[t].toArray()));
		}
		return C;
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

	public static double distanceTo(double a, double b) {
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
		int w = 1;
		double [][]sequences = new double[100][];
		for(int i=0;i<sequences.length;i++){
			sequences[i] = new double[20];
			for(int j=0;j<sequences[i].length;j++){
				sequences[i][j] = Math.cos(Math.random()*j/20.0*Math.PI) ;
			}
		}
		double [] averageSequence = new double[20];
		int choice = (int) Math.random()*100;
		for(int j=0;j<averageSequence.length;j++){
			averageSequence[j] = sequences[choice][j] ;
		}
		
		System.out.print("[");
		for(int j=0;j<averageSequence.length;j++){
			System.out.print(averageSequence[j]+" ");
		}
		System.out.println("]");
		
		for (int i=0; i<10; i++){
			DBA_update(averageSequence, sequences, w);
			
			System.out.print("[");
			for(int j=0;j<averageSequence.length;j++){
				System.out.print(averageSequence[j]+" ");
			}
			System.out.println("]");
		}
	}
}
