package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.evolution.sitemodel.DPSiteModel;
import beast.evolution.sitemodel.SiteModel;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 16/07/13
 * Time: 3:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class MixtureNtdSubstLikelihood extends DPTreeLikelihood{
    public void initAndValidate() throws Exception{
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);
        dpVal = dpValInput.get();
        if(!(m_pSiteModel.get() instanceof DPSiteModel)){
            throw new RuntimeException("DPSiteModel required for site model.");
        }
        dpSiteModel = (DPSiteModel) m_pSiteModel.get();


        alignment = m_data.get();
        int patternCount = alignment.getPatternCount();



        int[][] clusterWeights = new int[dpSiteModel.getSiteModelCount()][patternCount];

        int siteModelCount = dpSiteModel.getSiteModelCount();

        int siteCount = alignment.getSiteCount();
        for(int i = 0; i < siteCount; i++){
            //System.err.println("substModelIndices[i]: "+dpNtdSiteModel.getSubstCurrCluster(i));
            //System.out.println(dpSiteModel.getCurrCluster(i)+" "+alignment.getPatternIndex(i));
            //clusterWeights[dpSiteModel.getCurrCluster(i)][alignment.getPatternIndex(i)]++;
            clusterWeights[dpVal.getCurrCategory(i)][alignment.getPatternIndex(i)]++;
        }

        for(int i = 0; i < siteModelCount;i++){
            NewWVTreeLikelihood2 treeLik = new NewWVTreeLikelihood2(
                    clusterWeights[i],
                    alignment,
                    m_tree.get(),
                    useAmbiguitiesInput.get(),
                    dpSiteModel.getSiteModel(i),
                    m_pBranchRateModel.get());
            treeLiks.add(treeLik);


        }


    }

    public void addTreeLikelihood(){
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] patternWeights = new int[alignment.getPatternCount()];
        patternWeights[alignment.getPatternIndex(dpSiteModel.getLastDirtySite())] +=1;

        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood2 treeLik = new NewWVTreeLikelihood2(
                patternWeights,
                alignment,
                m_tree.get(),
                useAmbiguitiesInput.get(),
                siteModel,
                m_pBranchRateModel.get());
        try{
            /*treeLik.initByName(
                    "data", alignment,
                    "tree", treeInput.get(),
                    "siteModel", siteModel,
                    "branchRateModel", branchRateModelInput.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );*/

            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(dpSiteModel.getLastAddedIndex(),treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public void splitTreeLikelihood(){


        //alignment.printPattern();
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] clusterSites = dpVal.getClusterSites(dpSiteModel.getLastAddedIndex());

        //System.out.println("Last: "+dpSiteModel.getLastAddedIndex());

        int[] patternWeights = new int[alignment.getPatternCount()];

        int prevCluster = dpSiteModel.getDirtySiteModelIndex();
        //System.out.println("prevCluster: "+prevCluster);
        NewWVTreeLikelihood prevTreeLikelihood = treeLiks.get(prevCluster);
        for(int i = 0; i < clusterSites.length;i++){

            int patternIndex = alignment.getPatternIndex(clusterSites[i]);
            //System.out.println("clusterSites: "+clusterSites[i]+" PatIndex: "+patternIndex);
            patternWeights[patternIndex]++;
            //System.out.println(i+" splitting: "+patternWeights[patternIndex]);
            prevTreeLikelihood.removeWeight(patternIndex,1);
        }
        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood2 treeLik = new NewWVTreeLikelihood2(
                patternWeights,
                alignment,
                m_tree.get(),
                useAmbiguitiesInput.get(),
                siteModel,
                m_pBranchRateModel.get());
        try{

            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(dpSiteModel.getLastAddedIndex(),treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }



    }
}
