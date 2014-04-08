package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 3/07/13
 * Time: 5:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class MixedNtdBMAGammaBMASiteModel extends MixedNtdBMAGISiteModel{
    public Input<ParameterList> gammaSiteModelIndicatorListInput = new Input<ParameterList>(
            "gammaSiteModelIndicatorList",
            "A list of unique indicator values that determines whether the gamma site model includes the alpha or proportion invariant parameter.",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> invPrLogitInput = new Input<Boolean>(
            "invPrLogit",
            "Is transforming the invPr to logit space.",
            false
    );



    private ParameterList gammaSiteModelIndicatorList;
    private boolean invPrLogit;
    public MixedNtdBMAGammaBMASiteModel(){
        super();

    }


    public void initAndValidate() throws Exception{
        //super.initAndValidate();

        dpNtdBMA = dpNtdBMAInput.get();
        int ntdBMACount = dpNtdBMA.getDimension();
        ratePointers = ratePointersInput.get();
        rateList = rateListInput.get();
        alphaList = alphaListInput.get();
        invPrList = invPrListInput.get();


        if(rateList.getDimension() != dpNtdBMA.getDimension()){
            throw new RuntimeException("The number of clusters for rates and substitution models must be the same.");
        }

        //Setting up site models
        siteModels = new ArrayList<QuietSiteModel>();
        gammaCategoryCount =  gammaCategoryCountInput.get();
        invPrLogit = invPrLogitInput.get();

        gammaSiteModelIndicatorList = gammaSiteModelIndicatorListInput.get();
        for(int i = 0;i < ntdBMACount; i++){
            QuietRealParameter muParameter = rateList.getParameter(i);
            QuietRealParameter shapeParameter = alphaList.getParameter(i);
            QuietRealParameter invPrParameter = invPrList.getParameter(i);
            QuietRealParameter modelChoice = gammaSiteModelIndicatorList.getParameter(i);
            //System.out.println("gamma: "+gammaCategoryCount);
             QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(i),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    modelChoice,
                     invPrLogit
            );

           /*QuietSiteModel siteModel = new QuietSiteModel(
                    dpNtdBMA.getModel(i),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount
            );*/
            siteModels.add(siteModel);
        }
    }


    protected void addSiteModel(){
        try{

            QuietRealParameter muParameter = rateList.getParameter(rateList.getLastAddedIndex());
            QuietRealParameter shapeParameter = alphaList.getParameter(alphaList.getLastAddedIndex());
            QuietRealParameter invPrParameter = invPrList.getParameter(invPrList.getLastAddedIndex());
            QuietRealParameter modelChoice = gammaSiteModelIndicatorList.getParameter(gammaSiteModelIndicatorList.getLastAddedIndex());
            //System.out.println("Added: "+gammaSiteModelIndicatorList.getDimension());
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getLastAddedIndex()),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    modelChoice,
                    invPrLogit
                    );

            /*QuietSiteModel siteModel = new QuietSiteModel(
                    dpNtdBMA.getModel(dpNtdBMA.getLastAddedIndex()),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount
            );*/
            siteModels.add(dpNtdBMA.getLastAddedIndex(),siteModel);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }


    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        ChangeType substModelChangeType = dpNtdBMA.getChangeType();
        //System.out.println(rateList.somethingIsDirty() +" "+ dpNtdBMA.isDirtyCalculation());

        if(rateList.somethingIsDirty() && dpNtdBMA.isDirtyCalculation()){
            changeType = rateList.getChangeType();
            if(changeType != substModelChangeType){
                System.out.println(substModelChangeType+" "+changeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }


            if(changeType == ChangeType.ADDED||changeType == ChangeType.SPLIT){
                addSiteModel();

                //setupPointerIndices();
            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){
                //System.out.println(ChangeType.SPLIT_AND_VALUE_CHANGE);
                addSiteModel();
                //System.out.println("dirtyIndex: "+dpNtdBMA.getDirtyModelIndex());
                MCMCNodeFactory.checkDirtiness(siteModels.get(dpNtdBMA.getDirtyModelIndex()));

            }else if(changeType == ChangeType.REMOVED || changeType == ChangeType.MERGE){
                removeSiteModel(rateList.getRemovedIndex());
                //setupPointerIndices();
            }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){
                removeSiteModel(rateList.getRemovedIndex());
                MCMCNodeFactory.checkDirtiness(siteModels.get(dpNtdBMA.getDirtyModelIndex()));
            }else if(changeType == ChangeType.VALUE_CHANGED){

                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }else {
                this.changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }
                //setupPointerIndices();
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (ratePointers.somethingIsDirty()){

            changeType = ratePointers.getChangeType();
            if(changeType != substModelChangeType){
                System.out.println(changeType+" "+substModelChangeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }
            recalculate = true;

            //this.changeType = ChangeType.POINTER_CHANGED;


            //setupPointerIndices();

        }else{
            if(rateList.somethingIsDirty()){
                this.changeType = rateList.getChangeType();
                /*for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }*/
                recalculate = true;

            }else if(alphaList.somethingIsDirty()){
                this.changeType = alphaList.getChangeType();
                recalculate = true;

            }else if(invPrList.somethingIsDirty()){
                this.changeType = invPrList.getChangeType();
                recalculate = true;

            }else if(gammaSiteModelIndicatorList.somethingIsDirty()){
                this.changeType = gammaSiteModelIndicatorList.getChangeType();
                recalculate = true;

            }else if(dpNtdBMA.isDirtyCalculation() && dpNtdBMA.getChangeType() == ChangeType.VALUE_CHANGED){
                this.changeType = ChangeType.VALUE_CHANGED;
                /*for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }*/
                recalculate = true;
            }
            for(SiteModel siteModel:siteModels){

                MCMCNodeFactory.checkDirtiness(siteModel);
                //System.out.println("siteModel: "+siteModel.isDirtyCalculation());
            }
        }

        return recalculate;
    }


}
