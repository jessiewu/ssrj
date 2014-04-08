package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class MixedNtdBMAGISiteModel extends DPNtdRateSiteModel {
    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
            "alphaList",
            "A list of unique alpha values of Gamma distribution used to model sites rates.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPrListInput = new Input<ParameterList>(
        "invPrList",
        "a list of unique invariant proportion values used to model sites rates.",
            Input.Validate.REQUIRED
    );

    public Input<Integer> gammaCategoryCountInput =
            new Input<Integer>("gammaCategoryCount", "gamma category count (default=zero for no gamma)", 4);

    protected ParameterList alphaList;
    protected ParameterList invPrList;

    public MixedNtdBMAGISiteModel(){
        super();

    }
    protected int gammaCategoryCount;
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
        for(int i = 0;i < ntdBMACount; i++){
            QuietRealParameter muParameter = rateList.getParameter(i);
            QuietRealParameter shapeParameter = alphaList.getParameter(i);
            QuietRealParameter invPrParameter = invPrList.getParameter(i);
            QuietSiteModel siteModel = new QuietSiteModel(
                    dpNtdBMA.getModel(i),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount
            );
            siteModels.add(siteModel);
        }
    }

    protected void addSiteModel(){
        try{

            QuietRealParameter muParameter = rateList.getParameter(rateList.getLastAddedIndex());
            QuietRealParameter shapeParameter = alphaList.getParameter(alphaList.getLastAddedIndex());
            QuietRealParameter invPrParameter = invPrList.getParameter(invPrList.getLastAddedIndex());

            QuietSiteModel siteModel = new QuietSiteModel(
                    dpNtdBMA.getModel(dpNtdBMA.getLastAddedIndex()),
                    muParameter,
                    shapeParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount);
            siteModels.add(dpNtdBMA.getLastAddedIndex(),siteModel);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }


    public boolean requiresRecalculation(){

        boolean recalculate = false;
        ChangeType substModelChangeType = dpNtdBMA.getChangeType();
        if(rateList.somethingIsDirty() && dpNtdBMA.isDirtyCalculation()){
            changeType = rateList.getChangeType();
            if(changeType != substModelChangeType){
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }

            if(changeType == ChangeType.ADDED||changeType == ChangeType.SPLIT){
                addSiteModel();
            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){
                addSiteModel();
                MCMCNodeFactory.checkDirtiness(siteModels.get(dpNtdBMA.getDirtyModelIndex()));
            }else if(changeType == ChangeType.REMOVED || changeType == ChangeType.MERGE){
                removeSiteModel(rateList.getRemovedIndex());
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
            }
            recalculate = true;

        }else if (ratePointers.somethingIsDirty()){

            changeType = ratePointers.getChangeType();
            if(changeType != substModelChangeType){
                System.out.println(changeType+" "+substModelChangeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }
            recalculate = true;


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

    public void store(){
        //System.out.println("store");

        for(QuietSiteModel siteModel:siteModels){
            siteModel.store();
        }
        super.store();
    }

    public void restore(){
        //System.out.println("restore");
        super.restore();
        for(QuietSiteModel siteModel:siteModels){
            siteModel.restore();
        }
    }


}
