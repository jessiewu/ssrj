package beast.evolution.likelihood;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 16/07/13
 * Time: 3:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class NewWVTreeLikelihood2 extends NewWVTreeLikelihood{
    public NewWVTreeLikelihood2(int[] patternWeights,
                               Alignment data,
                               Tree tree,
                               boolean useAmbiguities,
                               SiteModel siteModel,
                               BranchRateModel.Base branchRateModel){
        super(patternWeights, data, tree, useAmbiguities, siteModel, branchRateModel);

    }
    @Override
    protected void setup() {
        int nodeCount = tree.getNodeCount();
        //m_siteModel = m_pSiteModel.get();
        m_siteModel.setDataType(data.getDataType());
        //m_substitutionModel = m_siteModel.m_pSubstModel.get();
        //m_substitutionModel = ((QuietSiteModel)m_siteModel).getSubstitutionModel();

        /*if (m_pBranchRateModel.get() != null) {
        	m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }*/
    	m_branchLengths = new double[nodeCount];
    	m_StoredBranchLengths = new double[nodeCount];

        int nStateCount = data.getMaxStateCount();
        //System.out.println("nStateCount: "+nStateCount);
        int nPatterns = data.getPatternCount();
        //System.out.println("patternWeights.length:"+patternWeights.length);

        boolean[] unmasked = new boolean[patternWeights.length];
        for(int i = 0; i < unmasked.length;i++){
                ///System.out.println(patternWeights[i]);
                unmasked[i] = patternWeights[i] > 0;
            }
        if (nStateCount == 4) {

            m_likelihoodCore = new WVLikelihoodCore4(unmasked);
        } else {
            m_likelihoodCore = new WVLikelihoodCore(nStateCount,unmasked);
        }
        //System.err.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        m_fProportionInvariant = m_siteModel.getProportianInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if(m_siteModel instanceof QuietGammaSiteBMA){
            calcConstantPatternIndices(nPatterns, nStateCount);

        }else if (m_fProportionInvariant > 0) {
        	calcConstantPatternIndices(nPatterns, nStateCount);
        }
        addedPatternIds = new int[nPatterns];
        initCore();

        m_fPatternLogLikelihoods = new double[nPatterns];
        storedPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        m_nMatrixSize = (nStateCount +1)* (nStateCount+1);
        m_fProbabilities = new double[(nStateCount +1)* (nStateCount+1)];
        //System.out.println("m_fProbabilities: "+m_fProbabilities.length+" "+m_data.get());
        Arrays.fill(m_fProbabilities, 1.0);

        if (data instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }
    }
}
