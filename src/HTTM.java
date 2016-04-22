

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;



public class HTTM {

	String dataName = "";
	
	 /**
     * document data (term lists)
     */
   // int[][] documents;

    /**
     * vocabulary size
     */
    int V;

    /**
     * number of topics
     */
    int K;

    /**
     * Dirichlet parameter (document--topic associations)
     */
    double alpha=.1;

    /**
     * Dirichlet parameter (topic--term associations)
     */
    double beta=.1;
    
    /**
     * topic assignments for each document based on DMM.
     */
    int z[];
    
    /**
     * topic assignments for each word based on LDA.
     */
    int zw[][];
    
    
    /**
     * number of documents in cluster z.
     */
   // int m_z[];
    
    /**
     * number of words in cluster z.
     */
    int n_z[];
    
    /**
     * cwt[i][j] number of instances of word i (term?) assigned to topic j.
     */
    //int[][] nw;

    
    /**
     * number of occurrences of word w in cluster z.
     */
    int n_w_z[][];
    
    /**
     * na[i][j] number of words in document i assigned to topic j.
     */
    int[][] nd;
    
    
    /**
     * number of words in document d.
     */
    int N_d[];
    
    /**
     * number of occurrences of word w in document d
     */
   // int N_w_d[][];
    
    double[][] thetasum;
    
    /**
     * cumulative statistics of phi
     */
    double[][] phisum;
 
    /**
     * size of statistics
     */
    int numstats;
    
    /**
     * the result of clustering
     */
    
  //  int cluster_doc[][];
    
    private static int THIN_INTERVAL = 20;

    /**
     * burn-in period
     */
    private static int BURN_IN = 100;

    /**
     * max iterations
     */
    private static int ITERATIONS = 1000;

    /**
     * sample lag (if -1 only one sample taken)
     */
    private static int SAMPLE_LAG;
    /**
     * the number of clusters
    */
    
    private static int dispcol = 0;
    
    int clustering;
    
	ArrayList<Integer> lablesArr = new ArrayList<Integer>();
	//ArrayList<String> gene ;
	ArrayList<ArrayList<Integer>> sData ;
	//ArrayList<ArrayList<Integer>> sDataNum;
	ArrayList<String> wordsArr = new ArrayList<String>();
	ArrayList<String> nWordsArr = new ArrayList<String>();
	ArrayList<Integer> IdList = new ArrayList<Integer>();
	
	int clu[];
	
	int threshold;
    /**
     * Initialise the Gibbs sampler with data.
     * 
     * @param V
     *            vocabulary size
     * @param data
     */
    public HTTM(int[][] documents, int V) {

        //this.documents = documents;
        this.V = V;
    }

    public HTTM(int thres, String dn) {

    	threshold = thres;
    	dataName = dn;
    }
    
    /**
     * Initialisation: Must start with an assignment of observations to topics ?
     * Many alternatives are possible, I chose to perform random assignments
     * with equal probabilities
     * 
     * @param K
     *            number of clusters
     * @return z assignment of clusters to documents
     */
    public void initialState(int K) {
     

        int M = sData.size();

        // initialise count variables.
        //m_z = new int[K];
        n_z = new int[K];
        n_w_z = new int[K][V];
        N_d = new int[M];
        //N_w_d = new int[M][V];
        
        // The z_i are are initialised to values in [1,K] to determine the
        // initial state of the Markov chain.

        z = new int[M];
        zw = new int[M][];
        nd = new int[M][K];
       // nw = new int[V][K];
        clu = new int[M];
        for (int m = 0; m < M; m++) 
        {
        	int N = sData.get(m).size(); //documents[m].length;
        	
        	N_d[m] = N;
        	
        	zw[m] = new int[N];
        	
        	if(N < threshold)
        	{
        		int topic = (int) (Math.random() * K);
           	 	z[m] = topic;
           	 //	m_z[topic]++;
           	 	n_z[topic] += N;
           	 
           	 	for(int n=0; n<N; n++)
           	 	{
           	 		n_w_z[topic][sData.get(m).get(n)]++;
           		 
           	 		//N_w_d[m][sData.get(m).get(n)]++;
           	 		
           	 		zw[m][n] = topic;
           	 		
           	 		//nw[sData.get(m).get(n)][topic]++;
           	 	}
           	 	
           	 	nd[m][topic] = N;
           	 	
        	}else
        	{
        		
        		for (int n = 0; n < N; n++) {
                    int topic = (int) (Math.random() * K);
                    zw[m][n] = topic;
                    // number of instances of word i assigned to topic j
                    n_z[topic]++;
                    
                    n_w_z[topic][sData.get(m).get(n)]++;
                   // N_w_d[m][sData.get(m).get(n)]++;
                    
                   
                    
                   //nw[sData.get(m).get(n)][topic]++;
                    // number of words in document i assigned to topic j.
                    nd[m][topic]++;
                    // total number of words assigned to topic j.
                    //nwsum[topic]++;
                }
        		
        		
        		/*double text_topic[] = new double[K];
        		int topic = -1;
        		double max = 0;
        		for(int i=0; i<K; i++)
        		{
        			text_topic[i] = (nd[m][i] + alpha)/(N + K * alpha);
        			if(text_topic[i] > max)
        			{
        				max = text_topic[i];
        				topic = i;
        			}
        		}
        		
        		z[m] = topic;
        		m_z[topic]++;*/
        	}
        	 
        }
        
       // cluster_doc = new int[K][M];
        clustering = 0;
    }
    
    /**
     * Configure the gibbs sampler
     * 
     * @param iterations
     *            number of total iterations
     * @param burnIn
     *            number of burn-in iterations
     * @param thinInterval
     *            update statistics interval
     * @param sampleLag
     *            sample interval (-1 for just one sample at the end)
     */
    public void configure(int iterations, int burnIn, int thinInterval,
        int sampleLag) {
        ITERATIONS = iterations;
        BURN_IN = burnIn;
        THIN_INTERVAL = thinInterval;
        SAMPLE_LAG = sampleLag;
    }
    
    /**
     * Sample a topic z_i from the full conditional distribution: p(z_i = j |
     * z_-i, w) = (n_-i,j(w_i) + beta)/(n_-i,j(.) + W * beta) * (n_-i,j(d_i) +
     * alpha)/(n_-i,.(d_i) + K * alpha)
     * 
     * @param m
     *            document
     */
    private int sampleFullConditionalDMM(int m) {
        // remove z_i from the count variables
    	//System.out.println("DMM");
        int topic = z[m];
       // m_z[topic]--;
		n_z[topic] -= N_d[m];
		nd[m][topic] -= N_d[m];
		for(int n=0; n<sData.get(m).size(); n++)
		{
			n_w_z[topic][sData.get(m).get(n)] --;
			//nw[sData.get(m).get(n)][topic] --;
		}
		
        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
        	
        	double wordsS = 1;
    		double wordsT = 1;
    		for(int n=0; n<sData.get(m).size(); n++)
    		{
    			wordsS *= (n_w_z[k][sData.get(m).get(n)] + beta);
    			wordsT *= (n_z[k] + V*beta);
    		}
        	
        	/*double wordsS = 0;
    		double wordsT = 0;
    		for(int n=0; n<sData.get(m).size(); n++)
    		{
    			wordsS += (n_w_z[k][sData.get(m).get(n)] + beta);
    			wordsT += (n_z[k] + V*beta);
    		}*/
        	/*double wordsS = 0;
    		double wordsT = 0;
    		for(int n=0; n<sData.get(m).size(); n++)
    		{
    			wordsS += (n_w_z[k][sData.get(m).get(n)] + beta);
    			wordsT += (n_z[k] + V*beta);
    		}*/
    		
    		
        	//p[k] = (m_z[k] + alpha) / (sData.size()-1+K*alpha) * wordsS / wordsT;
    		p[k] =  wordsS / wordsT;
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        
        double u = Math.random() * p[K - 1];
       // System.out.println("u: " + u);
       // System.out.println("p : " + p[p.length-1]);
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }

        // add newly estimated z_i to count variables
        //z[m] = topic;
       /* if(topic == K)
        {
        	topic = z[m-1];
        }*/
       // System.out.println("topic: " + topic);
       // m_z[topic]++;
        n_z[topic] += N_d[m];
        nd[m][topic] += N_d[m];
   	 
   	 	for(int n=0; n<sData.get(m).size(); n++)
   	 	{
   	 		n_w_z[topic][sData.get(m).get(n)]++;
   	 		zw[m][n] = topic;
   	 		
   	 		//nw[sData.get(m).get(n)][topic]++;
   	 	}
   		 
   	 	
        return topic;
    }

    
    private int sampleFullConditionalLDA(int m, int n) {

    	//System.out.println("LDA");
        // remove z_i from the count variables
        int topic = zw[m][n];
        n_w_z[topic][sData.get(m).get(n)]--;
        
        nd[m][topic]--;
        n_z[topic]--;
        //N_d[m]--;

        // do multinomial sampling via cumulative method:
        double[] p = new double[K];
        for (int k = 0; k < K; k++) {
            p[k] = (n_w_z[k][sData.get(m).get(n)] + beta) / (n_z[k] + V * beta)
                * (nd[m][k] + alpha);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < p.length; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalised p[]
        double u = Math.random() * p[K - 1];
        for (topic = 0; topic < p.length; topic++) {
            if (u < p[topic])
                break;
        }

        // add newly estimated z_i to count variables
        n_w_z[topic][sData.get(m).get(n)]++;
        nd[m][topic]++;
        n_z[topic]++;
        //N_d[m]++;

        return topic;
    }
    
    /**
     
     * @param K
     *            number of clusters
     * @param alpha
     *            symmetric prior parameter on document--cluster associations
     * @param beta
     *            symmetric prior parameter on cluster--term associations
     */
    private void gibbs(int K, double alpha, double beta) {
        this.K = K;
        this.alpha = alpha;
        this.beta = beta;

        // initial state of the Markov chain:
     // init sampler statistics
        if (SAMPLE_LAG > 0) {
            thetasum = new double[sData.size()][K];
            phisum = new double[K][V];
            numstats = 0;
        }
        
        initialState(K);

        for(int i=0; i<ITERATIONS; i++)
        {
        	for (int m = 0; m < sData.size(); m++)
        	{
        		//System.out.println("m : " + m);
        		
        		int N = sData.get(m).size();
        		
        		if(N < threshold)
        		{
        			int topic = sampleFullConditionalDMM(m);
        			z[m] = topic;
        		}else
        		{
        			 //m_z[z[m]]--;
        			 for (int n = 0; n < N; n++) 
        			 {
        				 int topic = sampleFullConditionalLDA(m, n);
                         zw[m][n] = topic;
        			 }
        		    /*double text_topic[] = new double[K];
             		int topic = -1;
             		double max = 0;
             		for(int j=0; j<K; j++)
             		{
             			text_topic[j] = (nd[m][j] + alpha)/(N + K * alpha);
             			if(text_topic[j] > max)
             			{
             				max = text_topic[j];
             				topic = j;
             			}
             		}
             		
             		z[m] = topic;
             		m_z[topic]++;*/
        		}
        	}
        	
        	if ((i < BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("B");
                dispcol++;
            }
            // display progress
            if ((i > BURN_IN) && (i % THIN_INTERVAL == 0)) {
                System.out.print("S");
                dispcol++;
            }
            // get statistics after burn-in
            if ((i > BURN_IN) && (SAMPLE_LAG > 0) && (i % SAMPLE_LAG == 0)) {
                updateParams();
                System.out.print("|");
                if (i % THIN_INTERVAL != 0)
                    dispcol++;
            }
            if (dispcol >= 100) {
                System.out.println();
                dispcol = 0;
            }
        }
        
       /* for(int i=0; i<z.length; i++)
        {
        	//System.out.println(i + ":" + z[i] + " ");
        	cluster_doc[z[i]][i] = 1;
        }
       // System.out.println();
        for(int i=0; i<cluster_doc.length; i++)
   	 	{
   		 	for(int j=0; j<cluster_doc[i].length; j++)
   		 	{
   		 		if(cluster_doc[i][j]!=0)
   		 		{
   		 			clustering++;
   		 			break;
   		 		}
   		 		
   		 	}
   		 	
   	 	}*/
    }
    
    /**
     * Add to the statistics the values of theta and phi for the current state.
     */
    private void updateParams() {
        
        for (int k = 0; k < K; k++) {
            for (int w = 0; w < V; w++) {
                phisum[k][w] += (n_w_z[k][w] + beta) / (n_z[k] + V * beta);
            }
        }
        for (int m = 0; m < sData.size(); m++) {
            for (int k = 0; k < K; k++) {
                thetasum[m][k] += (nd[m][k] + alpha) / (N_d[m] + K * alpha);
            }
        }
        
        numstats++;
    }
    
    public double[][] getPhi() {
        double[][] phi = new double[K][V];
        if (SAMPLE_LAG > 0) {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    phi[k][w] = phisum[k][w] / numstats;
                }
            }
        } else {
            for (int k = 0; k < K; k++) {
                for (int w = 0; w < V; w++) {
                    phi[k][w] = (n_w_z[k][w] + beta) / (n_z[k] + V * beta);
                }
            }
        }
        return phi;
    }

    
    
    
    public void readText( String path)
    {
    	String csvFile = path;
    	BufferedReader br = null;
		String line = "";
		String cvsSplitBy = " ";
		
		//label = new ArrayList<Integer>();
		
		//gene = new ArrayList<String>();
		
		sData = new ArrayList<ArrayList<Integer>>();
		
	
		try {
	 
			br = new BufferedReader(new FileReader(csvFile));
			
			
			//V = ge.length-1;
			
			HashMap<Integer,Integer> vArr = new HashMap<Integer,Integer>();
			int id=0;
			
			while ((line = br.readLine()) != null) {
	 
			        // use comma as separator
				String[] num = line.split(cvsSplitBy);
				
				//System.out.println(num[0]);
				//int laberNum = Integer.parseInt(num[0]);
				
				//label.add(laberNum);
					
				ArrayList<Integer> sample = new ArrayList<Integer>();
		
				
				for(int i=0; i<num.length; i++)
				{
					if(num[i].equals(""))
						continue;
					int word = Integer.parseInt(num[i]);
					
					if(!vArr.containsKey(word))
					{
						vArr.put(word, id);
						sample.add(id);
						nWordsArr.add(wordsArr.get(word));
						id++;
						
					}else
					{
						sample.add(vArr.get(word));
					}
					
				
				}
				sData.add(sample);
				
			}
			V = vArr.size();
			br.close();
	 
		} catch(Exception e)
		{
			e.printStackTrace();
		}
	 
		
    }
    
    public int[]  getMaxIndex(double []array, int top_k)
    {
    	double[] max = new double[top_k];
        int[] maxIndex = new int[top_k];
        Arrays.fill(max, Double.NEGATIVE_INFINITY);
        Arrays.fill(maxIndex, -1);

        top: for(int i = 0; i < array.length; i++) {
            for(int j = 0; j < top_k; j++) {
                if(array[i] > max[j]) {
                    for(int x = top_k - 1; x > j; x--) {
                        maxIndex[x] = maxIndex[x-1]; max[x] = max[x-1];
                    }
                    maxIndex[j] = i; max[j] = array[i];
                    continue top;
                }
            }
        }
        return maxIndex;

    }
    
    public void printText()
    {
    	for(int i=0; i<sData.size(); i++)
    	{
    		System.out.println(i + " :");
    		for(int j=0; j<sData.get(i).size(); j++)
    		{
    			System.out.print(nWordsArr.get(sData.get(i).get(j)) + " ");
    		}
    		System.out.println();
    	}
    }
    
    public void readWord( String path)
    {
    	
    	BufferedReader br = null;
		String line = "";
		String cvsSplitBy = " ";
		
		try {
	 
			br = new BufferedReader(new FileReader(path));

			//ArrayList<Integer> vArr = new ArrayList<Integer>();
			
			while ((line = br.readLine()) != null) {
	 
			        // use comma as separator
				String[] num = line.split(cvsSplitBy);
				
				//System.out.println(num[0]);
				//int laberNum = Integer.parseInt(num[0]);
				
				//label.add(laberNum);
				//wordId.put(Integer.parseInt(num[0]), num[1]);
				wordsArr.add(num[0]);
			}
	 
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	 
		
    }
    
    public int commTwoString(ArrayList<String> list1, ArrayList<String> list2)
    {
    	int count = 0;
    	
    	for(int i=0; i<list1.size(); i++)
    		if(list2.contains(list1.get(i)))
    			count++;
    	return count;
    }
    
   public double computPurity(String refPath, String comPath)
   {
	   //	String refFile = refPath;
   		BufferedReader brRef = null;
   		BufferedReader brCom = null;
		String line = "";
		String splitBy = " ";
		
		ArrayList<ArrayList<String>> refList = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> comList = new ArrayList<ArrayList<String>>();
	
		try {
	 
			brRef = new BufferedReader(new FileReader(refPath));
			brCom = new BufferedReader(new FileReader(comPath));

			while ((line = brRef.readLine()) != null) {
				ArrayList<String> ref = new ArrayList<String>();
				String[] words = line.split(splitBy);
				for(int i=0; i<words.length; i++)
					ref.add(words[i]);
				refList.add(ref);
			}
			while ((line = brCom.readLine()) != null) {
				ArrayList<String> com = new ArrayList<String>();
				String[] words = line.split(splitBy);
				for(int i=0; i<words.length; i++)
					com.add(words[i]);
				comList.add(com);
			}
			
			
		}catch(Exception e)
		{
			e.printStackTrace();
		}
		
		//double res = 0;
		
		int count = 0;
		for(int i=0; i<refList.size(); i++)
		{
			int max = 0;
			for(int j=0; j<comList.size(); j++)
			{
				int comm = commTwoString(refList.get(i),comList.get(j));
				if(comm>max)
					max = comm;
			}
			
			count += max;
		}
		
		return 1.0/(refList.size()*refList.get(0).size())*count;
   }
   
   public void readLable( String path)
   {
   	
   	BufferedReader br = null;
		String line = "";
		String cvsSplitBy = " ";
		
		try {
	 
			br = new BufferedReader(new FileReader(path));

			//ArrayList<Integer> vArr = new ArrayList<Integer>();
			
			while ((line = br.readLine()) != null) {
	 
			        // use comma as separator
				String[] num = line.split(cvsSplitBy);
				
				//System.out.println(num[0]);
				//int laberNum = Integer.parseInt(num[0]);
				
				//label.add(laberNum);
				//wordId.put(Integer.parseInt(num[0]), num[1]);
				lablesArr.add(Integer.parseInt(num[0]));
			}
			//V = wordsArr.size();
	 
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
   }
   
 
    public void mainFun(String path, String wordPath, int topics, int topWords) throws Exception
    {
    	readWord(wordPath);
    	readText(path);
    	
    	
    	
    	int M = sData.size();
    	
    	 int K = topics;
    	// V = gene2Id.size();
    	 //System.out.println("v : " + V);
         // good values alpha = 2, beta = .5
         double alpha = .1;
         double beta = .1;

         System.out.println("Topic Model Unbalanced using Gibbs Sampling.");

         // LdaGibbsSampler lda = new LdaGibbsSampler(documents, V);
          configure(1000, 500, 100, 5);
          
          ArrayList<Double> nmiArr = new ArrayList<Double>();
          
          String resPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/20News/Test/Unbal_topics=" + topics + "topWords="+topWords+".txt";
          
          FileWriter fw = new FileWriter(resPath);
  		  BufferedWriter bw = new BufferedWriter(fw);
  		  
          
          gibbs(K, alpha, beta);

         
          double[][] phi = getPhi();
         
          for(int i=0; i<phi.length; i++)
          {
       	 //  ArrayList<String> wArr = new ArrayList<String>();
       	   
       	   int[] maxIndices = getMaxIndex(phi[i],topWords);
              
       	   for(int j=0; j<maxIndices.length; j++)
       	   {
       		  // wArr.add(wordId.get(maxIndices[j]));
       		   bw.write(wordsArr.get(maxIndices[j]));
       		   bw.write(" ");
       	   }
       	   bw.newLine();
       	  
          }
         
        // bw.write(res);
         bw.close();
         fw.close();
        // test();
    }
    
    public double mainFun2(int times, int topics, int topWords) throws Exception
    {
    
    	 int K = topics;
          
    	 String resPath = "";
         if(times>0)
       	  //resPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/googlenews/TopWords/LDA_topics=" + topics + "topWords="+topWords+"times="+times+".txt";
         	resPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/"+dataName +"/TopWords/HTTM_topics=" + topics + "topWords="+topWords+"times="+times+".txt";
         else
       	  resPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/" + dataName +"/TopWords/HTTM_topics=" + topics + "topWords="+topWords+".txt";

          FileWriter fw = new FileWriter(resPath);
  		  BufferedWriter bw = new BufferedWriter(fw);
  		  
  		  double alpha = 0.1;
		  double beta = 0.1;
          
          gibbs(K, alpha, beta);

         
          double[][] phi = getPhi();
          
          ArrayList<ArrayList<Integer>> topWordsArr = new ArrayList<ArrayList<Integer>>();
         
          for(int i=0; i<phi.length; i++)
          {
       	 //  ArrayList<String> wArr = new ArrayList<String>();
       	   
	       	   	int[] maxIndices = getMaxIndex(phi[i],topWords);
	       	   	ArrayList<Integer> wArr = new ArrayList<Integer>(); 
	       	   	for(int j=0; j<maxIndices.length; j++)
	       	   	{
	       		  // wArr.add(wordId.get(maxIndices[j]));
	       	   		wArr.add(maxIndices[j]);
	       		   bw.write(nWordsArr.get(maxIndices[j]));
	       		   bw.write(" ");
	       	   	}
	       	   	topWordsArr.add(wArr);
	       	   	bw.newLine();
       	  
          }
         
          getCluster();
          double nmi = computNMI();
          System.out.println(nmi);
        bw.close();
        fw.close();
        return nmi;
        
        // test();
    }
    
    public double mainFun(String path, String lablePath, int topics) throws Exception
    {
    	
    	readText(path);
    	
    	readLable(lablePath);
    	
    	
    //	int M = sData.size();
    	
    	 int K = topics;
    	// V = gene2Id.size();
    	 //System.out.println("v : " + V);
         // good values alpha = 2, beta = .5
         double alpha = .1;
         double beta = .1;

         System.out.println("Topic Model Unbalanced using Gibbs Sampling.");

         // LdaGibbsSampler lda = new LdaGibbsSampler(documents, V);
          configure(1000, 100, 100, 5);
          
      //    ArrayList<Double> nmiArr = new ArrayList<Double>();
          
          String resPath = "C:/Users/qjp/Desktop/TopicModel/dataset/googlenews/result/UNtopics=" + topics + "NMI.txt";
          
          FileWriter fw = new FileWriter(resPath);
  		  BufferedWriter bw = new BufferedWriter(fw);
  		  
          
          gibbs(K, alpha, beta);

         
         // double[][] theta = getTheta();
          //  double[][] phi = getPhi();
           
           getCluster();
           double nmi = computNMI();
           System.out.println(nmi);
         bw.write(String.valueOf(nmi));
         bw.close();
         fw.close();
         return nmi;
        // test();
    }
    
    public void getCluster()
    {
    	for (int m = 0; m < sData.size(); m++)
    	{
    		//System.out.println("m : " + m);
    		
    		int N = sData.get(m).size();
    		
    		if(N < threshold)
    		{
    			clu[m] = z[m];
    		}else
    		{
    			double theta[] = new double[K];
    			int maxI = 0;
    			double val = 0.0;
    			for (int k = 0; k < K; k++) {
                    theta[k] = (nd[m][k] + alpha) / (N_d[m] + K * alpha);
                    if(theta[k]>val)
                    {
                    	maxI = k;
                    	val = theta[k];
                    }
                }
    			clu[m] = maxI;
    		}
    	}
    	
    }
    
    
    public ArrayList<Integer> getCluster2()
    {
    	ArrayList<Integer> arr = new ArrayList<Integer>();
    	
    	for (int m = 0; m < sData.size(); m++)
    	{
    		//System.out.println("m : " + m);
    		
    		int N = sData.get(m).size();
    		
    		if(N < threshold)
    		{
    			clu[m] = z[m];
    		}else
    		{
    			double theta[] = new double[K];
    			int maxI = 0;
    			double val = 0.0;
    			for (int k = 0; k < K; k++) {
                    theta[k] = (nd[m][k] + alpha) / (N_d[m] + K * alpha);
                    if(theta[k]>val)
                    {
                    	maxI = k;
                    	val = theta[k];
                    }
                }
    			clu[m] = maxI;
    			if(val<0.4)
    				arr.add(m);
    		}
    	}
    	
    	return arr;
    	
    }
    
    
    public double computNMI()
    {
    	//double res = 0;
    	
    	ArrayList<ArrayList<Integer>> textLabel = new ArrayList<ArrayList<Integer>>();
    	
    	ArrayList<Integer> labelId = new ArrayList<Integer>();
    	
    	for(int i=0; i< lablesArr.size(); i++)
    	{
    		int id = lablesArr.get(i);
    		
    		if(labelId.contains(id))
    		{
    			int index = labelId.indexOf(id);
    			
    			textLabel.get(index).add(i);
    		}else
    		{
    			ArrayList<Integer> subLabel = new ArrayList<Integer>();
    			subLabel.add(i);
    			labelId.add(id);
    			textLabel.add(subLabel);
    		}
    	}
    	
    	ArrayList<ArrayList<Integer>> clusterLabel = new ArrayList<ArrayList<Integer>>();
    	
    	ArrayList<Integer> clusterlId = new ArrayList<Integer>();
    	
    	for(int i=0; i<clu.length ; i++)
    	{
    		int id = clu[i];
    		
    		if(clusterlId.contains(id))
    		{
    			int index = clusterlId.indexOf(id);
    			
    			clusterLabel.get(index).add(i);
    		}else
    		{
    			ArrayList<Integer> subLabel = new ArrayList<Integer>();
    			subLabel.add(i);
    			clusterlId.add(id);
    			clusterLabel.add(subLabel);
    		}
    	}
    	
    	System.out.println(" the cluster number : " + clusterLabel.size());
    	
    	double comRes = 0;
    	
    	for(int i=0; i<textLabel.size(); i++)
    	{
    		for(int j=0; j<clusterLabel.size(); j++)
    		{
    			int common = commonArray(textLabel.get(i),clusterLabel.get(j));
    			
    			if(common!=0)
    				comRes += (double)common*Math.log((double)clu.length*common/(textLabel.get(i).size()*clusterLabel.get(j).size()));
    		}	
    	}
    	
    	double comL = 0;
    	for(int i=0; i<textLabel.size(); i++)
    	{
    		comL += (double)textLabel.get(i).size()*Math.log((double)textLabel.get(i).size()/clu.length);
    	}
    	
    	double comC = 0;
    	for(int j=0; j<clusterLabel.size(); j++)
    		comC += (double)clusterLabel.get(j).size()*Math.log((double)clusterLabel.get(j).size()/clu.length);
    	
    	//System.out.println(comRes + " " + comL + " "+ comC);
    	
    	comRes /= Math.sqrt(comL*comC);
    	for(int i=0; i<clusterLabel.size(); i++)
    	{
    		System.out.println(i + " " +clusterLabel.get(i).toString());
    	}
    	
    	return comRes;
    }
    
    public int commonArray(ArrayList<Integer> arr1, ArrayList<Integer> arr2)
    {
    	int count = 0;
    	for(int i=0; i<arr1.size(); i++)
    		if(arr2.contains(arr1.get(i)))
    			count++;
    	
    	return count;
    }
    
    public void mainFunForpercent(int percent, int times) throws Exception
    {
    
    	 int K = 100;
    	 int topWords = 20;
    	 
    	 
          String resPath = "C:/Users/qjp/Desktop/TopicModel/dataset/NIPS2Result/Unbal" +"_Percent=" + percent + "Times="+times+".txt";
          
          FileWriter fw = new FileWriter(resPath);
  		  BufferedWriter bw = new BufferedWriter(fw);
  		  
          
          gibbs(K, alpha, beta);

         
          double[][] phi = getPhi();
         
          for(int i=0; i<phi.length; i++)
          {
       	 //  ArrayList<String> wArr = new ArrayList<String>();
       	   
       	   int[] maxIndices = getMaxIndex(phi[i],topWords);
              
       	   for(int j=0; j<maxIndices.length; j++)
       	   {
       		  // wArr.add(wordId.get(maxIndices[j]));
       		   int ind = maxIndices[j];
       		   bw.write(wordsArr.get(ind));
       		   bw.write(" ");
       	   }
       	   bw.newLine();
       	  
          }
         
        // bw.write(res);
         bw.close();
         fw.close();
        // test();
    }
    
    public void percentIterMain() throws Exception 
    {
    	String path = "C:/Users/jipeng/Desktop/TopicModel/dataset/NIPS/PercentText/both";
    	File folder = new File(path);
		File[] listOfFiles = folder.listFiles();

		//int topics = 100;
		
		//int topWords = 20;
		
		String wordPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/NIPS/PercentText/word.txt";
		readWord(wordPath);
		
		configure(1000, 500, 100, 5);
		
		int k = 20;
		int topWords = 20;
		
		for(int i=5; i<=60; i=i+5)
		{
			String path2 = path + String.valueOf(i) + ".txt";
			readText(path2);
			mainFun2(i, k, topWords);
		}
		
    }
    
    
    public void iterMain(String path, String wordPath) throws Exception
    {
    	readWord(wordPath);
    	readText(path);
    	
    	
    	readLable("C:/Users/jipeng/Desktop/TopicModel/dataset/" + dataName + "/label.txt");
    	
    	String nmiPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/" + dataName + "/NMI/HTTM_NMI.txt";
    	
    	FileWriter fw = new FileWriter(nmiPath);
		BufferedWriter bw = new BufferedWriter(fw);
   	// V = gene2Id.size();
   	 //System.out.println("v : " + V);
        // good values alpha = 2, beta = .5
     //  double alpha = .1;
        //double beta = .1;
//
        System.out.println("Topic Model Unbalanced using Gibbs Sampling.");

        // LdaGibbsSampler lda = new LdaGibbsSampler(documents, V);
         configure(1000, 500, 100, 5);
         for(int times=1; times<=5; times++ ) 
    	for(int k=8; k<=20; k=k+20)
    	{
    		for(int topWords=10; topWords<=10; topWords=topWords+20)
    		{
    			double nmi = mainFun2(times, k, topWords);
         		bw.write(String.valueOf(nmi) + " ");
    		}
    	}
         
         bw.close();
         fw.close();
    }
    
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		int thres = 10; // less than 100  is short text; 
		String dataName = "NewsTweet";
		
		HTTM lda = new HTTM(thres,dataName);
		
		
       // ArrayList<Integer> classes = new ArrayList<Integer>();
       // int a[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
        //classes.add(0);
		try
		{
			  //String path = "C:/Users/qjp/Desktop/TopicModel/dataset/title.txt";
			int k = 20;
		int topWords = 100; 
			  //String path = "C:/Users/qjp/Desktop/TopicModel/dataset/title.txt";
			 // int k = 135;
			//int topWords =20; 
		//String path = "C:/Users/jipeng/Desktop/TopicModel/dataset/NIPS/ShortText/sentences.txt";
		//String wordPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/NIPS/ShortText/word.txt";
		//String path = "C:/Users/jipeng/Desktop/TopicModel/dataset/20News/both10.txt";
		//String wordPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/20News/word.txt";
			String path = "C:/Users/jipeng/Desktop/TopicModel/dataset/" + lda.dataName +"/text.txt";
			String wordPath = "C:/Users/jipeng/Desktop/TopicModel/dataset/" + lda.dataName +"/word.txt";
			 //lda.percentIterMain();
	          //lda.mainFun(path, lablePath, k);
		   lda.iterMain(path, wordPath);
	          //lda.mainFun(path, wordPath, k, topWords);
		}catch(Exception e)
		{
			e.printStackTrace();
		}
       
	}

}


