package beast.evolution.substitutionmodel;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;

import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import jdk.incubator.vector.*;
import java.util.concurrent.TimeUnit ;

public class BD extends SubstitutionModel.Base {
	public Input<RealParameter> nstate = new Input<RealParameter>("nstate", "same as what in BD model", Validate.REQUIRED);
	protected static int nrOfStates;
	protected static float[] binom_array;
	protected static float[] db_resultsi;
	protected static float[] db_resultsj;
	protected static float[] samplematrix;
	protected VectorSpecies<Float> SPECIES = FloatVector.SPECIES_512;

	protected Vector<Float> v0 = SPECIES.broadcast (0);
	protected Vector<Float> v1 = SPECIES.broadcast (1);
	protected static int upperBound;
	protected static int upperBoundr;
	public static int indexSeen = -5;
	public static FloatVector rangeSeen;
	public long totalTime1 = 0;
	public long totalTime2 = 0;

	public BD() {
		// this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
		frequenciesInput.setRule(Validate.OPTIONAL);
		try {
			// this call will be made twice when constructed from XML
			// but this ensures that the object is validly constructed for testing purposes.
			//System.out.println("Class Constructor");
			//System.out.println(nstate.get().getValue());
			//initAndValidate();
		} catch (Exception e) {
			e.printStackTrace();
			throw new RuntimeException("initAndValidate() call failed when constructing BD()");
		}
	}

	@Override
	public void initAndValidate() {
		int i, j, index;
		super.initAndValidate();
		//System.out.println("BD");
		//System.out.println(nstate.get().getValue());
		if (nstate.get() == null) {
			throw new IllegalArgumentException("number of states to consider is required");
		}
		nrOfStates = (int) Math.round(nstate.get().getValue());
		upperBound = SPECIES.loopBound(nrOfStates);
		//nrOfStates = 30;
		binom_array = new float[nrOfStates * nrOfStates];
		db_resultsi = new float[(nrOfStates+1) * (nrOfStates+1)];
		db_resultsj = new float[(nrOfStates+1) * (nrOfStates+1)];
		samplematrix = new float[(nrOfStates+1) * (nrOfStates+1)];
		for (i = 0; i < nrOfStates; ++i) {
			for (j = 0; j < nrOfStates; ++j) {
				index = i * nrOfStates + j;
				binom_array[index] = (float) binomialCoeff(i, j);
				db_resultsi[index] = i;
				db_resultsj[index] = j;
				samplematrix[index] = 0;
			}
		}
		//trMatrix = new double[(nrOfStates - 1) * (nrOfStates - 1)];
	}


	//public static final int nstates = 30;
	@Override
	public double[] getFrequencies() {
		return null;
	}

	public int getStateCount() {
		return nrOfStates;
	}
	//protected int nrOfStates = 30;


	public static float binomi(int n, int k) {
		return binom_array[n * nrOfStates + k];
	}

	public static double bd_prob(int child, int ancestor, double bd_rate, double distance) {
		double p = 0;
		int j;
		double p_j;
		int range = Math.min(child, ancestor) + 1;
		if (distance <=1) {
			for (j = 1; j < range; ++j ){
				p_j = (double) binomi(ancestor, j) * (double) binomi(child - 1, j - 1) * Math.pow(bd_rate * distance, child + ancestor -2*j);
				p += p_j;
			}
			p = p * Math.pow(bd_rate / (1 + distance * bd_rate), child + ancestor);
			//if (Math.pow(bd_rate / (1 + distance * bd_rate), child + ancestor) == 0) {
			//p = 0;
			//}
		}
		else {
			for (j = 1; j < range; ++j ){
				p_j = (double) binomi(ancestor, j) * (double) binomi(child - 1, j - 1) * Math.pow(bd_rate * distance, -2 * j);
				p += p_j;
			}
			p = p * Math.pow(bd_rate*distance / (1 + distance * bd_rate), child + ancestor);
			if (Math.pow(bd_rate*distance / (1 + distance * bd_rate), child + ancestor) == 0) {
				//p = 0;
			}
		}

		return p;
	}

	public FloatVector bd_probVectorized(FloatVector vbd_rate,FloatVector vdistance, FloatVector vprei, FloatVector vprej, FloatVector range,int index) {
		if(indexSeen==index) {
			return rangeSeen;
		}
		float p_j,p;
		int lj,r,index1,index2;
		FloatVector ve,b1,b2,sum,ve1,vbd,vd,vj,vi;
		for (r = 0; r < range.length(); ++r) {
			p=0;
			sum = FloatVector.zero(SPECIES);
			upperBound = SPECIES.loopBound((int) range.lane(r));
			//ve = FloatVector.fromArray(SPECIES, samplematrix, index)
			//float[] range = IntStream.rangeClosed(1, 10).toArray();;
			if (vdistance.lane(r) <= 1) {
				for (lj = 1; lj < upperBound; lj+=SPECIES.length()) {
					index1 = (int) vprei.lane(r) * nrOfStates +  lj;
					index2 = ((int) vprej.lane(r) - 1) * nrOfStates + ( lj - 1);
					ve1 = FloatVector.fromArray(SPECIES, samplematrix, index1);
					b1 = FloatVector.fromArray(SPECIES, binom_array, index1);
					b2 = FloatVector.fromArray(SPECIES,binom_array,index2);
					vbd = ve1.add(vbd_rate.lane(r));
					vd = ve1.add(vdistance.lane(r));
					vj = ve1.add(vprej.lane(r));
					vi = ve1.add(vprei.lane(r));
					sum = b1.fma(b2.mul(vbd.mul(vd).pow(vj.add(vi.add(-2*lj)))),sum);
				}
				p = sum.reduceLanes(VectorOperators.ADD);
				for (; lj < r; ++lj ){
					p_j = (float) binomi((int) vprei.lane(r), lj) * (float) binomi(((int) vprej.lane(r) - 1), lj - 1) * (float) Math.pow(vbd_rate.lane(r) * vdistance.lane(r), (int) vprej.lane(r) + (int) vprei.lane(r)-2*lj);
					p += p_j;
				}
				p = p * (float) Math.pow(vbd_rate.lane(r) / (1 + vdistance.lane(r) * vbd_rate.lane(r)), (int) vprej.lane(r) + (int) vprei.lane(r));
				range=range.withLane(r,p);


			}
			else {
				for (lj = 1; lj < upperBound; lj+=SPECIES.length()) {
					index1 = (int) vprei.lane(r) * nrOfStates +  lj;
					index2 = ((int) vprej.lane(r) - 1) * nrOfStates + ( lj - 1);
					ve1 = FloatVector.fromArray(SPECIES, samplematrix, index1);
					b1 = FloatVector.fromArray(SPECIES, binom_array, index1);
					b2 = FloatVector.fromArray(SPECIES,binom_array,index2);
					vbd = ve1.add(vbd_rate.lane(r));
					vd = ve1.add(vdistance.lane(r));
					vj = ve1.add(vprej.lane(r));
					vi = ve1.add(vprei.lane(r));
					sum = b1.fma(b2.mul(vbd.mul(vd).pow(-2*lj)),sum);
				}

				p = sum.reduceLanes(VectorOperators.ADD);
				for (; lj < r; ++lj ){
					p_j = (float) binomi((int) vprei.lane(r), lj) * (float) binomi(((int) vprej.lane(r) - 1), lj - 1) * (float) Math.pow(vbd_rate.lane(r) * vdistance.lane(r), -2*lj);
					p += p_j;
				}
				p = p * (float) Math.pow(vbd_rate.lane(r)*vdistance.lane(r) / (1 + vdistance.lane(r) * vbd_rate.lane(r)), (int) vprej.lane(r) + (int) vprei.lane(r));
				range=range.withLane(r,p);
			}
		}
		indexSeen = index;
		rangeSeen = range;
		return range;
	}

	protected boolean checkTransitionMatrix(double[] matrix) {
		double sum = 0;
		int i, j;
		int index;
		for (i = 0; i < nrOfStates; ++i) {
			for (j = 0; j < nrOfStates; ++j) {
				index = i * nrOfStates + j;
				sum = sum + matrix[index];
			}
			if (sum > 1.01 | sum < 0.95) {
				//System.out.println("current index:" + i);
				//System.out.println(sum);
				return true;
			}
			sum = 0;
		}
		return true;

	}

	public double getProbability(int i, int j,double distance, double bd_rate){
		if(i == 0){
			if (j == 0) {
				return 1;
			}
			else {
				return 0;
			}
		}
		else if(j == 0){
			return Math.pow((bd_rate * distance) / (1 + bd_rate * distance), i);
		}
		else if (i == 1){
			return Math.pow(distance, j - 1) / Math.pow((1 + distance), j + 1);
		}
		else{
			return bd_prob(j, i, bd_rate, distance);
		}
	}

	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		float bd_rate = 1;
		int index;
		int i=0,j=0;
		double prob;
		float distance = ((float) startTime - (float) endTime) * (float) rate;
		FloatVector vm,vprei,vprej,vdistance,vbd_rate,range;
		double[] dvm;
		long startTime1 = System.nanoTime();
		for(i =0; i<upperBound;i++){
			for(j=0;j<upperBound;j+=SPECIES.length()){
				index = i * nrOfStates + j;
				vm = FloatVector.fromArray(SPECIES, samplematrix, index);
				vprei = FloatVector.fromArray(SPECIES, db_resultsi, index);
				vprej = FloatVector.fromArray(SPECIES, db_resultsj, index);
				//vm=vm.sub(vm);
				vdistance = vm.add(distance);
				vbd_rate = vm.add(bd_rate);
				range = vprei.min(vprej).add(1);
				vm=vm.add(1,vprei.eq(0).and(vprej.eq(0)))
						.add(0,vprei.eq(0).and(vprej.eq(0).not()))
						.add(vbd_rate.mul(vdistance).div(vbd_rate.mul(vdistance).add(1))
								.pow(vprei),vprei.eq(0).not().and(vprej.eq(0)))
						.add(vdistance.pow(vprej.sub(1)).div(vdistance.add(1).pow(vprej.add(1))),
								vprei.eq(1).and(vprej.eq(0).not()))
						.add(bd_probVectorized(vbd_rate,vdistance, vprei, vprej, range,index).mul(vbd_rate.div(vbd_rate.mul(vdistance).add(1)).pow(vprej.add(vprei)))
								,vprei.eq(0).not().and(vprei.eq(1).
								not()).and(vprej.eq(0).not()).and(vdistance.eq(1).or(vdistance.lt(1))))
						.add(bd_probVectorized(vbd_rate,vdistance, vprei, vprej, range,index).mul(vbd_rate.mul(vdistance).div(vbd_rate.mul(vdistance).add(1)).pow(vprej.add(vprei)))
								,vprei.eq(0).not().and(vprei.eq(1).
								not()).and(vprej.eq(0).not()).and(vdistance.eq(1).not().and(vdistance.lt(1).not())));
				dvm = vm.toDoubleArray();
				for(int dmvindex = 0;dmvindex<dvm.length;dmvindex++){
					matrix[index + dmvindex] = dvm[dmvindex];
				}
			}

			for (; j < nrOfStates ; ++j) {
				index = i * nrOfStates + j;
				matrix[index] = getProbability(i,j,distance,bd_rate);
			}
		}
		System.out.println("SPECIES");
		System.out.println(SPECIES);
		long endTime1 = System.nanoTime();
		long duration1 = (endTime1 - startTime1);
		totalTime1 = duration1 + totalTime1;
		System.out.println("BD VECTORIZED TIME");
		System.out.println(totalTime1);


		long startTime2 = System.nanoTime();
		for (i = 0; i < nrOfStates ; ++i) {
			for (j = 0; j < nrOfStates ; ++j) {
				index = i * nrOfStates + j;
				matrix[index] = getProbability(i,j,distance,bd_rate);
			}
		}
		long endTime2 = System.nanoTime();
		long duration2 = (endTime2 - startTime2);
		totalTime2 = totalTime2 + duration2;
		System.out.println("BD NON VECTORIZED TIME");
		System.out.println(totalTime2);

	}

	public int binomialCoeff(int n, int k) {
		int res = 1;

		// Since C(n, k) = C(n, n-k)
		if (k > n - k)
			k = n - k;

		// Calculate value of
		// [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
		for (int i = 0; i < k; ++i) {
			res *= (n - i);
			res /= (i + 1);
		}

		return res;
	}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		// TODO Auto-generated method stub
		return dataType instanceof IntegerData;
	}
}

