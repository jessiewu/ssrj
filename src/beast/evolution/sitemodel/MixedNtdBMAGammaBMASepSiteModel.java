package beast.evolution.sitemodel;

import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class MixedNtdBMAGammaBMASepSiteModel extends DPNtdBMAGIBMASepSiteModel{
    public Input<ParameterList> siteModelBoundaryListInput = new Input<ParameterList>(
            "siteModelBoundaryList",
            "A list of the boundaries of the site rate models",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> siteModelBoundaryPointersInput = new Input<DPPointer>(
            "siteModelBoundaryPointers",
            "A list of pointers referring to a set of unique parameters.",
            Input.Validate.REQUIRED
    );


    private ParameterList siteModelBoundaryList;
    private DPPointer siteModelBoundaryPointers;

    public MixedNtdBMAGammaBMASepSiteModel(){
        ratesPointersInput.setRule(Input.Validate.OPTIONAL);

    }


    public void initAndValidate() throws Exception{
        siteModelBoundaryList = siteModelBoundaryListInput.get();
        siteModelBoundaryPointers = siteModelBoundaryPointersInput.get();
        siteModels = new ArrayList<QuietSiteModel>();

        dpNtdBMA = dpNtdBMAInput.get();

        ratesList = ratesListInput.get();
        eltCount = siteModelBoundaryPointers.getDimension();
        if(dpNtdBMA.getCount() != eltCount){
            throw new RuntimeException("The number of model elements is different to the number of rate elements.");
        }
        dpValSubst = dpValSubstInput.get();
        dpValRate = dpValRateInput.get();

        substClusterLimit = substClusterLimitInput.get();
        ratesClusterLimit = ratesClusterLimitInput.get();
        gammaCategoryCount =  gammaCategoryCountInput.get();
        siteModelChoiceList = siteModelChoiceListInput.get();
        alphaList = alphaListInput.get();
        invPrList = invPrListInput.get();
        invPrLogit = invPrLogitInput.get();

        setup();
    }

    /*
     * Setup the site model and weight matrices.
     */
    public void setup(){
        siteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        storedSiteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        siteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        storedSiteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        //substWeights = new int[substClusterLimit];
        //ratesWeights = new int[ratesClusterLimit];
        clusterMap = new int[2][siteModelBoundaryPointers.getDimension()];
        storedClusterMap = new int[2][siteModelBoundaryPointers.getDimension()];

        //int[] substModelPointerIndicies = dpNtdBMA.getPointerIndices();
        int substModelIndex = -1;
        int rateIndex = -1;

        for(int i = 0; i < eltCount; i++){
            substModelIndex = dpNtdBMA.getModelBySiteIndexIDNumber(i);
            //System.out.println("eltCount: "+substModelIndex);
            rateIndex = siteModelBoundaryPointers.indexInList(i,siteModelBoundaryList);
            //System.out.print(rateIndex);
            if(siteModelsMatrix[substModelIndex][rateIndex] == null){

                try{

                    QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                            dpNtdBMA.getModel(substModelIndex),
                            ratesList.getParameter(rateIndex),
                            alphaList.getParameter(rateIndex),
                            invPrList.getParameter(rateIndex),
                            true,
                            gammaCategoryCount,
                            siteModelChoiceList.getParameter(rateIndex),
                            invPrLogit
                    );
                    siteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    storedSiteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    siteModels.add(siteModel);

                }catch(Exception e){
                    throw new RuntimeException(e);
                }

            }
            siteModelWeights[substModelIndex][rateIndex]++;
            storedSiteModelWeights[substModelIndex][rateIndex]++;
            //substWeights[substModelIndex]++;
            //ratesWeights[rateIndex]++;
            clusterMap[NTDBMA][i] = substModelIndex;
            clusterMap[RATES][i] = rateIndex;
            storedClusterMap[NTDBMA][i] = substModelIndex;
            storedClusterMap[RATES][i] = rateIndex;
        }

    }

    protected void handleSplit(int changedInput) throws Exception{

        int ntdBMAIDNum;
        int muIDNum;
        if(changedInput == NTDBMA){
            //Retreive the new NtdBMA model
            ntdBMAIDNum = dpNtdBMA.getModelByListIndexIDNumber(dpNtdBMA.getLastAddedIndex());
            //System.out.println("ntdBMAIDNum: "+ntdBMAIDNum+ " "+dpNtdBMA.getModelByListIndexIDNumber(dpNtdBMA.getLastAddedIndex()-1));
            //Retrieve all the (dirty) sites that uses the newly added model
            clusterSites = dpValSubst.getClusterSites(dpNtdBMA.getLastAddedIndex());
            //System.out.println(dpNtdBMA.getLastAddedIndex());

            for(int dirtySite: clusterSites){
                //System.out.println("dirtysite: "+dirtySite);

                //Get the rate of the each of the dirty sites
                muIDNum = siteModelBoundaryPointers.getParameter(dirtySite).getIDNumber();
                //Update the substModel and rate combination
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }

        }else if(changedInput == RATES){

            muIDNum = ratesList.getParameterIDNumber(ratesList.getLastAddedIndex());
            clusterSites = dpValRate.getClusterSites(ratesList.getLastAddedIndex());
            for(int dirtySite:clusterSites){
                //System.out.println("dirtysite: "+dirtySite);
                ntdBMAIDNum = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite)).getIDNumber();
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

    }


    public void addSiteModel(int changedInput, int lastDirtySite) throws Exception{
        int ntdBMACluster;
        int rateCluster;
        if(changedInput == NTDBMA){
            ntdBMACluster = dpNtdBMA.getLastModel().getIDNumber();
            rateCluster = siteModelBoundaryPointers.getParameter(lastDirtySite).getIDNumber();

        }else if(changedInput == RATES){
            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = ratesList.getLastParameter().getIDNumber();
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

        updateMap(lastDirtySite, ntdBMACluster, rateCluster);

    }


    public void handleMerge(int changedInput) throws Exception{

        int ntdBMAIDNum;
        int muIDNum;
        if(changedInput == NTDBMA){
            int removedSubstModelIndex = dpNtdBMA.getRemovedIndex();
            //System.out.println(removedSubstModelIndex);
            clusterSites = dpValSubst.getStoredClusterSites(removedSubstModelIndex);
            //This the substModel cluster that takes up the members of the removed substModel cluster
            ntdBMAIDNum = dpNtdBMA.getDirtyModelIDNumber();
            //System.out.println("ntdBMAIDNum: "+dpNtdBMA.getDimension());
            for(int dirtySite: clusterSites){
                muIDNum = siteModelBoundaryPointers.getParameterIDNumber(dirtySite);
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }

        }else if(changedInput == RATES){
            //The rate cluster that takes up the member of the removed rate cluster
            muIDNum = ratesList.getDirtyParameterIDNumber();

            int removedRateIndex = ratesList.getRemovedIndex();
            clusterSites = dpValRate.getStoredClusterSites(removedRateIndex);

            for(int dirtySite:clusterSites){
                ntdBMAIDNum = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite)).getIDNumber();
                //System.out.println(ntdBMAIDNum +" " +muIDNum);
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);

            }
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

    }


    public void removeSiteModel(int changedInput, int lastDirtySite) throws Exception {
        //int[] siteClusterMap = new int[2];
        int ntdBMACluster;
        int rateCluster;
        if(changedInput == NTDBMA){

            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = siteModelBoundaryPointers.getParameter(lastDirtySite).getIDNumber();

        }else if(changedInput == RATES){

            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = siteModelBoundaryPointers.getParameterIDNumber(lastDirtySite);

        }else{
            throw new RuntimeException("Can only remove clusters for either ntdBMA model or rates.");
        }


        //System.out.println("changedInput: "+changedInput+" "+ntdBMAId+" "+rateId);
        //siteModels.remove(siteModels.indexOf(siteModelsMatrix[ntdBMAId][rateId]));
        //siteModelsMatrix[ntdBMAId][rateId] = null;

        updateMap(lastDirtySite, ntdBMACluster, rateCluster);

    }


    protected void handleMultiPointerChanges(int changedInput) throws Exception{

        if(changedInput == NTDBMA){
            clusterSites = dpValSubst.getLastDirtySites();
        }else if(changedInput == RATES){
            clusterSites = dpValRate.getLastDirtySites();
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

        for(int dirtySite: clusterSites){
            SwitchingNtdBMA ntdBMA = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite));
            QuietRealParameter muParameter = siteModelBoundaryPointers.getParameter(dirtySite);

            int ntdBMACluster = ntdBMA.getIDNumber();
            int rateCluster = muParameter.getIDNumber();

            //update mapping and weights
            updateWeights(dirtySite, ntdBMACluster, rateCluster);



            //updateModelMatrix(dirtySite);



        }
        for(int dirtySite: clusterSites){
            updateModelMatrix(dirtySite);

        }

    }

    /*
     * Handles changes in pointers that does not involve changes in the number of clusters
     */
    protected void handlePointerChange(int lastDirtySite) throws Exception{

        //Previous cluster ids of the last dirty site
        int prevNtdBMAIdNum = clusterMap[NTDBMA][lastDirtySite];
        int prevRateIdNum = clusterMap[RATES][lastDirtySite];

        //Current cluster ids of the last dirty site
        SwitchingNtdBMA ntdBMA = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)); //todo
        int listIndex = siteModelBoundaryPointers.indexInList(lastDirtySite,siteModelBoundaryList);
        QuietRealParameter muParameter = ratesList.getParameter(listIndex);
        QuietRealParameter alphaParameter = alphaList.getParameter(listIndex);
        QuietRealParameter invPrParameter = invPrList.getParameter(listIndex);
        QuietRealParameter siteModelChoice = siteModelChoiceList.getParameter(listIndex);

        //int[] siteClusterMap = new int[2];
        int ntdBMACluster = ntdBMA.getIDNumber();
        int rateCluster = muParameter.getIDNumber();



        //update mapping and weights
        updateMap(lastDirtySite, ntdBMACluster, rateCluster);



        //If the propose combination is new then create a new site model
        if(siteModelsMatrix[ntdBMACluster][rateCluster] == null){
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    ntdBMA,
                    muParameter,
                    alphaParameter,
                    invPrParameter,
                    true,
                    gammaCategoryCount,
                    siteModelChoice,
                    invPrLogit
            );
            siteModelsMatrix[ntdBMACluster][rateCluster] = siteModel;
            siteModels.add(siteModel);
        }

        //If the previous combination no longer has any weight then remove
        //System.out.println(prevNtdBMAIdNum+" "+prevRateIdNum);

    }


    public void updateModelMatrix(int siteIndex) throws Exception{
        //Remove site model if it has zero pattern weight
        if(siteModelWeights[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] == 0
                && siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] != null){
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]]));
            siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            int listIndex = siteModelBoundaryPointers.indexInList(siteIndex,siteModelBoundaryList);
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),
                    ratesList.getParameter(listIndex),
                    alphaList.getParameter(listIndex),
                    invPrList.getParameter(listIndex),
                    true,
                    gammaCategoryCount,
                    siteModelChoiceList.getParameter(listIndex),
                    invPrLogit
            );
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

    }


    /*
     * Update cluster mapping
     */
    public void updateMap(int siteIndex, int ntdBMACluster, int rateCluster) throws Exception{
        //System.out.println("siteModel dim: "+siteModels.size());
        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("NTDBMA: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }

        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored NTDBMA: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/

        mappingChanged = true;

        //Update the weight matrix
        moveWeight(
                clusterMap[NTDBMA][siteIndex],  //current substModel
                clusterMap[RATES][siteIndex],   //current rate
                ntdBMACluster,                  //potentially new substModel
                rateCluster,                    //potentially new rate
                1
        );
        //System.out.println(clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
        //System.out.println(ntdBMACluster+" "+rateCluster);

        int prevNtdBMAId = clusterMap[NTDBMA][siteIndex];
        int prevRateId = clusterMap[RATES][siteIndex];

        //Update mapping
        clusterMap[NTDBMA][siteIndex] = ntdBMACluster;
        clusterMap[RATES][siteIndex] = rateCluster;

        //Stored map
        //storedClusterMap[NTDBMA][siteIndex] = prevNtdBMAId;
        //storedClusterMap[RATES][siteIndex] = prevRateId;

        //System.out.println("flag1: "+siteModelWeights[prevNtdBMAId][prevRateId]);


        //Remove site model if it has zero pattern weight
        if(siteModelWeights[prevNtdBMAId][prevRateId] == 0){
            //System.out.println("remove");
            //System.out.println("siteModels: "+siteModels.size());
            //System.out.println(siteModelsMatrix[prevNtdBMAId][prevRateId]);
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[prevNtdBMAId][prevRateId]));
            siteModelsMatrix[prevNtdBMAId][prevRateId] = null;
        }
        //System.out.println("siteModel dim: "+siteModels.size());

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);


            int listIndex = siteModelBoundaryPointers.indexInList(siteIndex,siteModelBoundaryList);
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);

            QuietGammaSiteBMA siteModel = new QuietGammaSiteBMA(
                    dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),
                    ratesList.getParameter(listIndex),
                    alphaList.getParameter(listIndex),
                    invPrList.getParameter(listIndex),
                    true,
                    gammaCategoryCount,
                    siteModelChoiceList.getParameter(listIndex),
                    invPrLogit
            );
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;

            siteModels.add(siteModel);
        }

        //System.out.println("siteModel dim: "+siteModels.size());

        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("RATE: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }
        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored RATE: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/




    }






    public void handlePointerSwap(int siteIndex1, int siteIndex2) throws Exception{
        throw new RuntimeException("Not applicable.");

    }





    public boolean requiresRecalculation(){

        boolean recalculate = false;
        mappingChanged = false;
        try{

            //ChangeType substModelChangeType = dpNtdBMA.getChangeType();
            if(ratesList.somethingIsDirty() ||
                    alphaList.somethingIsDirty() ||
                    invPrList.somethingIsDirty() ||
                    siteModelChoiceList.somethingIsDirty() ||
                    siteModelBoundaryPointers.somethingIsDirty()||
                    dpNtdBMA.isDirtyCalculation()){
                //ChangeType inputChangeType = null;


                //Check whether subst model or rates have changed and
                //retrieve change type.
                if(dpNtdBMA.isDirtyCalculation()){

                    changeType = dpNtdBMA.getChangeType();
                    changedInput = NTDBMA;
                    lastDirtySite = dpNtdBMA.getLastDirtySite();

                }else{
                    if(ratesList.somethingIsDirty()){
                        changeType = ratesList.getChangeType();
                    }else if(alphaList.somethingIsDirty()){
                        changeType = alphaList.getChangeType();
                    }else if(invPrList.somethingIsDirty()){
                        changeType = invPrList.getChangeType();
                    }else if(siteModelChoiceList.somethingIsDirty()){
                        changeType = siteModelChoiceList.getChangeType();
                    }else{
                        changeType = siteModelBoundaryPointers.getChangeType();
                    }
                    lastDirtySite = siteModelBoundaryPointers.getLastDirty();
                    changedInput = RATES;
                }


                //System.out.println("changeType: "+this.changeType);
                if(changeType == ChangeType.ADDED){
                    //System.out.println("ADD");
                    addSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.SPLIT){

                    handleSplit(changedInput);
                }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){
                    handleSplit(changedInput);
                    checkSiteModelsDirtiness(changedInput);



                }else if(changeType == ChangeType.REMOVED){
                    //System.out.println("REMOVED");
                    removeSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.MERGE){

                    handleMerge(changedInput);
                }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){

                    handleMerge(changedInput);
                    checkSiteModelsDirtiness(changedInput);

                }else if(changeType == ChangeType.POINTER_CHANGED){
                    //System.out.println("POINTER_CHANGED");

                    handlePointerChange(lastDirtySite);
                }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){

                    handleMultiPointerChanges(changedInput);


                }else if(changeType == ChangeType.POINTERS_SWAPPED){
                    if(changedInput == NTDBMA){
                        clusterSites = dpNtdBMA.getSwappedSites();

                    }else{
                        clusterSites = siteModelBoundaryPointers.getSwappedSites();
                    }

                    if(clusterSites[0] != clusterSites[1]){
                        handlePointerSwap(clusterSites[0], clusterSites[1]);
                        //handlePointerChange(clusterSites[0]);
                        //handlePointerChange(clusterSites[1]);

                    }

                }else if(changeType == ChangeType.VALUE_CHANGED){
                    //System.out.println("VALUE_CHANGED");
                    //When there's a value change, the cluster assignment doesn't change.
                    checkSiteModelsDirtiness(changedInput);

                }else {
                    this.changeType = ChangeType.ALL;
                    for(int i = 0; i < siteModelsMatrix.length; i++){
                        for(int j = 0; j < siteModelsMatrix[i].length; j++){
                            if(siteModelsMatrix[i][j] != null){
                                MCMCNodeFactory.checkDirtiness(siteModelsMatrix[i][j]);

                            }

                        }
                    }
                    //setupPointerIndices();
                }


                recalculate = true;
                //System.err.println("dirty1");

            }

            /*if(siteModels.size()< siteModelBoundaryList.getDimension() ||
                    siteModels.size() < dpNtdBMA.getDimension()){
                System.out.println("problem: "+siteModels.size()+" "+
                dpNtdBMA.getDimension()
                +" "+siteModelBoundaryList.getDimension());
                System.out.println(changedInput+" "+changeType);
                throw new RuntimeException("");
            }*/



        }catch(Exception e){
            throw new RuntimeException(e);
        }
        //System.out.println("siteModel changeType: "+changeType);
        return recalculate;
    }






}
