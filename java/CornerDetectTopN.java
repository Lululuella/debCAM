
import edu.rit.numeric.NonNegativeLeastSquares;
import java.util.List;

class Convset implements Comparable<Convset> {
	private double normsqrsum;
	private int[] CornerIdx;
	
	public Convset (double norm, int[] Idx) {
		this.normsqrsum = norm;
		
		int len = Idx.length;
		this.CornerIdx = new int[len];
		System.arraycopy(Idx, 0, this.CornerIdx, 0, len);
	}
	
	@Override
	public int compareTo(Convset o) {
		double delta = this.normsqrsum - o.normsqrsum;
	    if(delta > 0) return 1;
	    if(delta < 0) return -1;
	    return 0;
	}
	
	public double getNorm() {
        return normsqrsum;
    }
	
	public int[] getCornerIdx() {
        return CornerIdx;
    }
}


public class CornerDetectTopN {
	
	private double[][] X;
	private int K;
	private int M;
	private int L;
	private int N;
	private NonNegativeLeastSquares nnls;
	private FixSizedPriorityQueue<Convset> pq;
	private List<Convset> topconv;


	public CornerDetectTopN (double[][] X, int K, int N) {
		this.K = K;
		this.M = X.length;
		this.L = X[0].length;
		this.N = N;
		this.X = new double[M][L];
		for (int i = 0; i < M; i ++) {
			System.arraycopy(X[i], 0, this.X[i], 0, L);
		}

		nnls = new NonNegativeLeastSquares(M+1,K);
		
		this.pq = new FixSizedPriorityQueue<Convset>(N);
		
		 
	}
	
	private void combine(int[] wholeset, int n) {  
        
        if(null == wholeset || wholeset.length == 0 || n <= 0 || n > wholeset.length)  
            return;  
              
        int[] combset = new int[n];  
        getCombination(wholeset, n , 0, combset, 0);  
    }  
  
    private void getCombination(int[] wholeset, int n, int begin, int[] combset, int index) {  
          
        if(n == 0){
        	if (K != index) {
        		throw new IllegalArgumentException
				("K != index");			
			}
        	
        	double normsqrsum_tmp = 0;
        	for (int i = 0; i < L; i++) {
        		boolean innerflag = true;
        		for (int j = 0; j < K; j++) {
        			if (i == combset[j]) {
        				innerflag = false;
        				break;
        			}  
        		}
        		
        		if (!innerflag) continue;
        		
        		for (int k = 0; k < K; k++) {        		
            		for(int j = 0; j < M; j++)
            			nnls.a[j][k] = 1E-5 * X[j][combset[k]];
            		nnls.a[M][k] = 1;
            	}
        		
        		for(int j = 0; j < M; j++)
        			nnls.b[j] = 1E-5 * X[j][i];
        		nnls.b[M] = 1;
        		
        		nnls.solve();
        		normsqrsum_tmp += nnls.normsqr;
        	}
        	
        	Convset conv = new Convset(normsqrsum_tmp, combset);
        	
        	pq.add(conv);
      	
            return;  
        }  
              
        for(int i = begin; i < wholeset.length; i++){  
              
        	combset[index] = wholeset[i];  
            getCombination(wholeset, n-1, i+1, combset, index+1);  
        }  
          
    } 
    
	public boolean search() {
		
		int[] wholeset = new int[L]; 
		for(int i = 0; i < L; i++){
			wholeset[i] = i;
		}
		combine(wholeset, K);
		
		topconv = pq.sortedList();
		
		return true;
	}
	
    public int[][] getTopNConv() {
    	int[][] idx = new int[topconv.size()][K];
    	int i = 0;
		
    	for (Convset item : topconv) {
    		System.arraycopy(item.getCornerIdx(), 0, idx[i], 0, K);
    		i++;
    	}
    	
    	for(i = 0; i < idx.length; i++) {
    		for(int j = 0; j < idx[0].length; j++) {
    			idx[i][j] = idx[i][j]+1;
        	}
    	}
    	
    	return idx;	
    }
    
    public double[] getTopNConvErr() {
    	double[] norm = new double[topconv.size()];
    	int i = 0;
		
    	for (Convset item : topconv) {
    		norm[i] = item.getNorm();
    		i++;
    	}
    	return norm;	
    }
    
        
    public static void main(String[] args) {
    	double[][] X = {{0.1,0.2,1,0,0,0.5,0.3},{0.1,0.7,0,1,0,0.5,0.3},{0.8,0.1,0,0,1,0,0.4}};
    	CornerDetectTopN ex = new CornerDetectTopN(X,3,4);
    	ex.search();
    	double[] norm = ex.getTopNConvErr();
    	for(int i = 0; i < norm.length; i++) {
    		System.out.print(norm[i] + ", ");
    	}
    	System.out.println();
    	
    	int[][] idx = ex.getTopNConv();
    	for(int i = 0; i < idx.length; i++) {
    		for(int j = 0; j < idx[0].length; j++) {
        		System.out.print(idx[i][j] + ", ");
        	}
    		System.out.println();
    	}
    	System.out.println();

	}
	

}
