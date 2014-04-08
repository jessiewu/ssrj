package beast.app.seqgen;

import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.alignment.Alignment;
import beast.util.MersenneTwisterFast;
import beast.util.Randomizer;
import beast.util.TreeParser;
import beast.util.XMLProducer;
import com.sun.org.apache.bcel.internal.generic.NEW;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 29/10/2011
 * Time: 4:42:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class HeterogeneousSequenceSimulator extends SequenceSimulator{

    SiteModel[] siteModels;
    int[] siteModelAssign;
    protected double[][][] m_probabilities;
    Alignment data;
    public HeterogeneousSequenceSimulator(
            Alignment data,
            Tree tree,
            SiteModel[] siteModels,
            BranchRateModel branchRateModel,
            int siteModelAssign[],
            int sequenceLength,
            int stateCount,
            int categoryCount) {
        this.data = data;
    	m_tree = tree;
    	this.siteModels = siteModels;
    	m_branchRateModel = branchRateModel;
        this.siteModelAssign = siteModelAssign;
    	m_sequenceLength = sequenceLength;
    	m_stateCount = stateCount;
        m_categoryCount = categoryCount;
        m_probabilities = new double[siteModels.length][m_categoryCount][m_stateCount * m_stateCount];
    }

    public Alignment simulate() throws Exception {
    	Node root =  m_tree.getRoot();


    	double [][] categoryProbs = new double[siteModels.length][];
        for(int i = 0; i < categoryProbs.length;i++){
            categoryProbs[i]=siteModels[i].getCategoryProportions(root);
            /*for(double prop:categoryProbs[i]){
                System.out.println(prop);
            }*/
        }
    	int [] category  = new int[m_sequenceLength];
    	for (int i  = 0; i < m_sequenceLength; i++) {

    		category[i] = Randomizer.randomChoicePDF(categoryProbs[siteModelAssign[i]]);
            System.out.println(categoryProbs[siteModelAssign[i]][0]+" "+category[i]);
    	}

        //Specifies the root states
        double [][] frequencies = new double[siteModels.length][];
        for(int i = 0; i < siteModels.length;i++){
            frequencies[i] = siteModels[i].getSubstitutionModel().getFrequencies();
        }

    	int [] seq = new int[m_sequenceLength];
    	for (int i  = 0; i < m_sequenceLength; i++) {
        	seq[i] = Randomizer.randomChoicePDF(frequencies[siteModelAssign[i]]);
    	}


    	Alignment alignment = new Alignment();
          //alignment.setDataType(m_siteModel.getFrequencyModel().getDataType());

    	traverse(root, seq, category, alignment);



    	return alignment;
    }

    protected void traverse(Node node, int [] parentSequence, int [] category, Alignment alignment) throws Exception {
		for (int iChild = 0; iChild < 2; iChild++) {
			Node child = (iChild == 0 ? node.getLeft() : node.getRight());
            for(int i = 0; i < siteModels.length;i++){
                for (int j = 0; j < m_categoryCount; j++) {
            	    getTransitionProbabilities(m_tree, child,i, j, m_probabilities[i][j]);
                }
            }

        	int [] seq = new int[m_sequenceLength];
    		double [] cProb = new double[m_stateCount];
        	for (int i  = 0; i < m_sequenceLength; i++) {
                if(category[i] == 0){
                    throw new RuntimeException("");
                }
                System.arraycopy(m_probabilities[siteModelAssign[i]][category[i]], parentSequence[i]*m_stateCount, cProb, 0, m_stateCount);

            	seq[i] = Randomizer.randomChoicePDF(cProb);
        	}
            //System.out.println();
            if (child.isLeaf()) {
            	alignment.m_pSequences.setValue(intArray2Sequence(seq, child), alignment);
            } else {
            	traverse(child, seq, category, alignment);
            }
		}
	}

    /** get transition probability matrix for particular rate category **/
    protected void getTransitionProbabilities(Tree tree, Node node, int siteModelIndex,int rateCategory, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (m_branchRateModel == null ? 1.0 : m_branchRateModel.getRateForBranch(node));

        if(node.isRoot()){
            branchRate = 1.0;
        }else{
            String s = node.m_sMetaData;
            if(s != null){
                branchRate = (Double.parseDouble(s.substring(s.indexOf("rate")).split("=")[1]));
            //System.out.println(branchRate);
            }
        }

        branchRate *= siteModels[siteModelIndex].getRateForCategory(rateCategory, node);
        //System.out.println(branchRate);
        //System.out.println(siteModels[siteModelIndex].getRateForCategory(rateCategory, node));
        siteModels[siteModelIndex].getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);


    }

    public static int[][] getCombinations(
            int[] indicies1,
            int[] indicies2,
            HashMap<Integer,HashMap<Integer,Integer>> comboMap){
        int length;
        if(indicies1.length != indicies2.length){
            throw new RuntimeException("The lengths of the two int arrays must be equal.");
        }else{
            length = indicies1.length;

        }
        int comboCount = 0;
        //HashMap<Integer, HashMap<Integer,Integer>> comboMap = new HashMap<Integer, HashMap<Integer,Integer>>();
        //HashMap<Integer, ArrayList<Integer>> comboOrderMap = new HashMap<Integer, ArrayList<Integer>>();
        for(int i = 0; i < length;i++){
            if(comboMap.containsKey(indicies1[i])){
                if(!comboMap.get(indicies1[i]).containsKey(indicies2[i])){
                    //System.out.println("hello?");
                    comboMap.get(indicies1[i]).put(indicies2[i], comboCount++);

                }
            }else{
                HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
                map.put(indicies2[i],comboCount++);
                comboMap.put(indicies1[i],map);
            }
        }

        int[][] combo = new int[comboCount][2];
        int k  = 0;
        Integer[] uniqueIndices1 = comboMap.keySet().toArray(new Integer[comboMap.size()]);
        for(int i = 0; i < uniqueIndices1.length;i++){
            Set keys = comboMap.get(uniqueIndices1[i]).keySet();
            Integer[] keyArray = new Integer[keys.size()];
            keys.toArray(keyArray);
            for(int key: keyArray){
                combo[k][0] = uniqueIndices1[i];
                combo[k++][1] = key;
            }

        }
        return combo;

    }

    Sequence intArray2Sequence(int [] seq, Node node) throws Exception {
		DataType dataType = data.getDataType();
		String sSeq = dataType.state2string(seq);
//    	}
    	List<Sequence> taxa = data.m_pSequences.get();
    	String sTaxon = taxa.get(node.getNr()).m_sTaxon.get();
        //System.out.println("seq: "+ sSeq.toString());
		return new Sequence(sTaxon, sSeq.toString());
    } // intArray2Sequence

    public static int[] processPointers(String line) throws Exception{


        String[] substStr = line.split("\\t");
        int[] subst = new int[substStr.length-1];
        for(int i = 0;i < subst.length;i++){
            subst[i] = (int)Double.parseDouble(substStr[i+1]);
        }

        return subst;

    }

    public static double[] processPointersDouble(String line) throws Exception{
        String[] substStr = line.split("\\t");
        double[] subst = new double[substStr.length-1];
        for(int i = 0;i < subst.length;i++){
            subst[i] = Double.parseDouble(substStr[i + 1]);
        }

        return subst;

    }


    public static void main(String[] args){
        String[] taxonNames = new String[]{"Canis","Felis","Homo","Macaca","Microcebus","Mus","Ochotona","Oryctolagus","Otolemur","Pan","Pongo","Rattus"};


        Randomizer.setSeed(127);
        int counter = 0;
        try{

            BufferedReader treeReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_1_re.trees"));

            BufferedReader paramListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_paramList_1_re.log"));
            BufferedReader modelListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_modelList_1_re.log"));
            BufferedReader freqsListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_freqsList_1_re.log"));
            BufferedReader alphaListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_ratesList_1_re.log"));
            BufferedReader invPrListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_ratesList_1_re.log"));
            BufferedReader siteModelListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_ratesList_1_re.log"));
            BufferedReader ratesListReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_ratesList_1_re.log"));
            BufferedReader boundaryReader = new BufferedReader(new FileReader("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_ratesList_1_re.log"));

             Alignment data = getAlignment(taxonNames);

            String treeStr = "";

            String paramListStr = "";
            String modelListStr = "";
            String freqsListStr = "";
            String alphaListStr = "";
            String invPrListStr = "";
            String siteModelListStr = "";
            String ratesListStr = "";
            String boundaryStr = "";
            while((treeStr = treeReader.readLine())!= null){

                paramListStr = paramListReader.readLine();
                modelListStr = modelListReader.readLine();
                freqsListStr = freqsListReader.readLine();
                alphaListStr = alphaListReader.readLine();
                invPrListStr = invPrListReader.readLine();
                siteModelListStr = siteModelListReader.readLine();
                ratesListStr = ratesListReader.readLine();
                boundaryStr = boundaryReader.readLine();



                treeStr = labelTree(taxonNames,treeStr);
                TreeParser tree = new TreeParser();
                tree.initByName("taxa", data, "newick", treeStr);


                double[] paramList = processPointersDouble(paramListStr);
                double[] modelList = processPointersDouble(modelListStr);
                double[] freqsList = processPointersDouble(freqsListStr);
                double[] alphaList = processPointersDouble(alphaListStr);
                double[] invPrList = processPointersDouble(invPrListStr);
                double[] siteModelList = processPointersDouble(siteModelListStr);
                double[] ratesList = processPointersDouble(ratesListStr);
                double[] boundaryList = processPointersDouble(boundaryStr);

                SiteModel[] siteModels = new SiteModel[modelList.length];


                for(int i = 0; i < siteModels.length;i++){
                    QuietRealParameter logKappa = new QuietRealParameter(new Double[]{paramList[i*5]});
                    QuietRealParameter logTN = new QuietRealParameter(new Double[]{paramList[i*5+1]});
                    QuietRealParameter logAC = new QuietRealParameter(new Double[]{paramList[i*5+2]});
                    QuietRealParameter logAT = new QuietRealParameter(new Double[]{paramList[i*5+3]});
                    QuietRealParameter logGC = new QuietRealParameter(new Double[]{paramList[i*5+4]});
                    QuietRealParameter model = new QuietRealParameter(new Double[]{modelList[i]});
                    QuietRealParameter freqs = new QuietRealParameter(new Double[]{freqsList[i*4],freqsList[i*4+1],freqsList[i*4+2],freqsList[i*4+3]});
                    SwitchingNtdBMA substModel = new SwitchingNtdBMA(logKappa,logTN,logAC,logAT,logGC,model,freqs);
                    RealParameter alpha = new RealParameter(new Double[]{alphaList[i]});
                    RealParameter invPr = new RealParameter(new Double[]{invPrList[i]});
                    RealParameter rate = new RealParameter(new Double[]{ratesList[i]});
                    siteModels[i] = new SiteModel();
                    siteModels[i].initByName(
                            "substModel", substModel,
                            "shape", alpha,
                            "proportionInvariant", invPr,
                            "mutationRate", rate,
                            "gammaCategoryCount", 4);
                }


                int[] siteModelIndex = new int[data.getSiteCount()];
                int start, end;
                for(int i = 0; i< boundaryList.length;i++){
                    start = (int)boundaryList[i];
                    end = (int)boundaryList[i++];
                    i = i + 1;

                    for(int j = start; j <= end; j++){
                        siteModelIndex[j] = i;
                    }
                }




                HeterogeneousSequenceSimulator simulator = new HeterogeneousSequenceSimulator(
                        data,
                        tree,
                        siteModels,
                        null,
                        siteModelIndex,
                        468,
                        4,
                        1);

                Alignment alignment = simulator.simulate();
                alignment.initAndValidate();
                System.out.println(alignment.getNrTaxa());
                List<Sequence> seqs = alignment.m_pSequences.get();
                Sequence seq;

                /*PrintWriter pw = new PrintWriter("C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\SubstBMA\\xml\\data\\mammal\\model\\mammal_468_ntd_rate_sep_sim_"+(counter+1)+".txt");
                for(int i = 0; i < alignment.getNrTaxa();i++){
                    seq = seqs.get(i);
                    pw.println(seq.m_sTaxon.get()+"\t"+seq.m_sData.get());
                }
                pw.close();
                */
                produceXML(
                        alignment,
                        "C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\SubstBMA\\xml\\data\\mammal\\model\\template.xml",
                        "C:\\Users\\Jessie Wu\\Documents\\Jessie\\Research\\BEAST\\BEAST_logs\\model\\mammal_468_ntd_rate_sep_sim_"+(counter+1)+".xml",
                        (counter+1),
                        Randomizer.getSeed()
                );
                counter++;

                //Randomizer.setSeed((new MersenneTwisterFast()).getSeed());
                Randomizer.setSeed(Math.round(Randomizer.nextDouble()));
                System.out.println(Randomizer.getSeed());

            }

            treeReader.close();
            paramListReader.close();
            modelListReader.close();
            freqsListReader.close();
            alphaListReader.close();
            invPrListReader.close();
            ratesListReader.close();
            siteModelListReader.close();
        }catch(Exception e ){
            throw new RuntimeException(e);

        }
    }

    public static Alignment getAlignment(String[] taxonNames) throws Exception {
        Alignment data = new Alignment();
        for(int i = 0; i < taxonNames.length;i++){
            Sequence seq = new Sequence(taxonNames[i],"A");
            data.m_pSequences.setValue(seq, data);
        }

        data.initByName(
                "dataType", "nucleotide"
        );


        return data;
    }


    public static Alignment getAlignment(List<String> taxonNames) throws Exception {
        Alignment data = new Alignment();
        for(int i = 0; i < taxonNames.size();i++){
            Sequence seq = new Sequence(taxonNames.get(i),"A");
            data.m_pSequences.setValue(seq, data);
        }

        data.initByName(
                "dataType", "nucleotide"
        );


        return data;
    }

    public static String labelTree(String[] taxonNames, String tree){
        for(int i = 0;i < taxonNames.length;i++){
            tree = tree.replaceFirst(i + ":", taxonNames[i] + ":");
        }
        System.out.println(tree);
        return tree;

    }

    public static String labelTree2(String[] taxonNames, String tree){
        for(int i = 0;i < taxonNames.length;i++){

            String s = "("+(i+1)+"#";
                String s1 = "("+taxonNames[i]+"#";
                tree = tree.replace(s, s1);

        }
        tree = tree.replaceAll("#","");
        tree = tree.substring(tree.indexOf("("));
        System.out.println(tree);
        return tree;

    }

    public static String labelTree3(String[] taxonNames, String tree){
        tree = tree.replaceAll(":","#");
        for(int i = 0;i < taxonNames.length;i++){

            if(tree.contains("("+(i+1)+"#")){
                String s = "("+(i+1)+"#";
                String s1 = "("+taxonNames[i]+"#";
                tree = tree.replace(s, s1);
            }else if (tree.contains(","+(i+1)+"#")){
                tree = tree.replaceFirst(","+(i+1)+"#", ","+taxonNames[i]+"#");
            }
        }
        tree = tree.replaceAll("#",":");
        tree = tree.substring(tree.indexOf("("));
        System.out.println(tree);
        return tree;

    }

    public static String labelTree2(List<String> taxonNames, String tree){
        for(int i = 0;i < taxonNames.size();i++){


            tree = tree.replaceFirst((i+1)+"#", taxonNames.get(i));
        }
        tree = tree.replaceAll("#","");
        tree = tree.substring(tree.indexOf("("));
        //System.out.println(tree);
        return tree;

    }

    public static void produceXML(
            Alignment alignment,
            String templateFile,
            String outputFile,
            int xmlID,
            long seed) throws Exception{
        BufferedReader templateReader = new BufferedReader(new FileReader(templateFile));
        PrintWriter pw = new PrintWriter(outputFile);
        String line = "";
        while((line = templateReader.readLine()) != null){
            if(line.contains("[Alignment]")){
                pw.println("<!-- The sequence alignment -->");
                pw.println("<!-- Generated with seed: "+ seed +" -->");
                pw.println("<!-- ntaxa="+alignment.getNrTaxa()+" nchar="+alignment.getSiteCount()+" -->");
                pw.println("<!-- npat="+alignment.getPatternCount()+" -->");

                List<Sequence> seqs = alignment.m_pSequences.get();
                Sequence seq;
                pw.println("\t<data id=\"alignment\" dataType=\"nucleotide\">");
                for(int i = 0; i < alignment.getNrTaxa();i++){
                    seq = seqs.get(i);
                    pw.println("\t\t<sequence id=\"seq_"+seq.m_sTaxon.get()+"\" taxon=\""+seq.m_sTaxon.get()+"\" totalcount=\"4\" value=\""+seq.m_sData.get()+"\"/>");
                }
                pw.println("\t</data>");

            }else if(line.contains("<logger")){
                pw.println(line.replace("_x.","_"+xmlID+"."));

            }else {
                pw.println(line);

            }

        }
        pw.close();

    }



    public static void produceXML(
            Alignment alignment,
            String templateFile,
            String outputFile,
            String xmlID,
            long seed,
            String treeStr) throws Exception{
        BufferedReader templateReader = new BufferedReader(new FileReader(templateFile));
        PrintWriter pw = new PrintWriter(outputFile);
        String line = "";
        while((line = templateReader.readLine()) != null){
            if(line.contains("[Alignment]")){
                pw.println("<!-- The sequence alignment -->");
                pw.println("<!-- Generated with seed: "+ seed +" -->");
                pw.println("<!-- ntaxa="+alignment.getNrTaxa()+" nchar="+alignment.getSiteCount()+" -->");
                pw.println("<!-- npat="+alignment.getPatternCount()+" -->");

                List<Sequence> seqs = alignment.m_pSequences.get();
                Sequence seq;
                pw.println("\t<data id=\"alignment\" dataType=\"nucleotide\">");
                for(int i = 0; i < alignment.getNrTaxa();i++){
                    seq = seqs.get(i);
                    pw.println("\t\t<sequence id=\"seq_"+seq.m_sTaxon.get()+"\" taxon=\""+seq.m_sTaxon.get()+"\" totalcount=\"4\" value=\""+seq.m_sData.get()+"\"/>");
                }
                pw.println("\t</data>");

            }else if(line.contains("<logger")){
                pw.println(line.replace("_x.","_"+xmlID+"."));
            }else if(line.contains("[tree]")){

                pw.println(line.replace("[tree]",treeStr));

            }else {
                pw.println(line);

            }

        }
        pw.close();

    }


    public static void produceXML(
            Alignment alignment,
            String templateFile,
            String outputFile,
            String xmlID,
            long seed) throws Exception{
        BufferedReader templateReader = new BufferedReader(new FileReader(templateFile));
        PrintWriter pw = new PrintWriter(outputFile);
        String line = "";
        while((line = templateReader.readLine()) != null){
            if(line.contains("[Alignment]")){
                pw.println("<!-- The sequence alignment -->");
                pw.println("<!-- Generated with seed: "+ seed +" -->");
                pw.println("<!-- ntaxa="+alignment.getNrTaxa()+" nchar="+alignment.getSiteCount()+" -->");
                pw.println("<!-- npat="+alignment.getPatternCount()+" -->");

                List<Sequence> seqs = alignment.m_pSequences.get();
                Sequence seq;
                pw.println("\t<data id=\"alignment\" dataType=\"nucleotide\">");
                for(int i = 0; i < alignment.getNrTaxa();i++){
                    seq = seqs.get(i);
                    pw.println("\t\t<sequence id=\"seq_"+seq.m_sTaxon.get()+"\" taxon=\""+seq.m_sTaxon.get()+"\" totalcount=\"4\" value=\""+seq.m_sData.get()+"\"/>");
                }
                pw.println("\t</data>");

            }else if(line.contains("<logger")){
                pw.println(line.replace("_x.","_"+xmlID+"."));

            }else {
                pw.println(line);

            }

        }
        pw.close();

    }


}
