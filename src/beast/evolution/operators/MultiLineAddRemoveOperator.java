package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ConditionalCategoricalDistribution;
import beast.math.distributions.DirichletDistribution;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 23/10/13
 * Time: 2:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiLineAddRemoveOperator extends Operator {
    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );



    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
                "alphaList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPListInput = new Input<ParameterList>(
                "invPrList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );

    public Input<Boolean> invPrLogitInput = new Input<Boolean>(
                "invPrLogit",
                "An ParameterList object that contains a list of parameters.",
                false
        );

    public Input<ParameterList> modelListInput = new Input<ParameterList>(
                "modelList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );

    public Input<ParameterList> siteModelListInput = new Input<ParameterList>(
                "gammaSiteModelIndicatorList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );

    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "Reference pointers that refers to the category assigned",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> boundaryPointsInput = new Input<ParameterList>(
            "boundaryPoints",
            "A list of points that splits the line/queue.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> parametricDistributionInput = new Input<ParametricDistribution>(
            "distr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<DirichletDistribution> freqsDistributionInput = new Input<DirichletDistribution>(
            "freqsDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> ratesDistributionInput = new Input<ParametricDistribution>(
            "rateDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> alphaDistributionInput = new Input<ParametricDistribution>(
            "alphaDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> invPrDistributionInput = new Input<ParametricDistribution>(
            "invPrDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );





    public Input<ConditionalCategoricalDistribution> substModelIndicatorDistrInput = new Input<ConditionalCategoricalDistribution>(
            "substModelIndicatorDistr",
            "The distribution to be sampled from given the current value.",
            Input.Validate.REQUIRED
    );

    public Input<ConditionalCategoricalDistribution> siteModelIndicatorDistrInput = new Input<ConditionalCategoricalDistribution>(
            "siteModelIndicatorDistr",
            "The distribution to be sampled from given the current value.",
            Input.Validate.REQUIRED
    );

    //private double logSplitJacobian;
    //private double logMergeJacobian;
    //public static final double LOG_MERGE_JACOBIAN = Math.log(0.5);

    private boolean invPrLogit;

    private ParametricDistribution parametricDistribution;
    private DirichletDistribution freqsDistribution;
    private ParametricDistribution alphaDistribution;
    private ParametricDistribution invPrDistribution;
    private ParametricDistribution ratesDistribution;
    private ConditionalCategoricalDistribution substModelIndicatorDistr;
    private ConditionalCategoricalDistribution siteModelIndicatorDistr;
    private ParameterList boundaryPoints;
    private DPPointer pointers;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
        //logSplitJacobian = parameterListInput.get().getParameterDimension()*Math.log(2);
        //logMergeJacobian = -logSplitJacobian;
        freqsDistribution = freqsDistributionInput.get();
        ratesDistribution = ratesDistributionInput.get();
        alphaDistribution = alphaDistributionInput.get();
        invPrDistribution = invPrDistributionInput.get();
        substModelIndicatorDistr = substModelIndicatorDistrInput.get();
        siteModelIndicatorDistr = siteModelIndicatorDistrInput.get();
        invPrLogit = invPrLogitInput.get();
    }

    public double proposal() {
        double logq = 0.0;


        try{
            ParameterList parameterList = parameterListInput.get();
            pointers = pointersInput.get();
            boundaryPoints = boundaryPointsInput.get();
            boolean splitMove;
            splitMove = Randomizer.nextBoolean();
            ///System.out.println(splitMove+" "+parameterList.getDimension());

            if(!splitMove && boundaryPoints.getDimension() == 1){

                 return Double.NEGATIVE_INFINITY;


            }else if (splitMove && boundaryPoints.getDimension() == pointers.getDimension()){

                return Double.NEGATIVE_INFINITY;

            }
            //System.out.println("Hello!: "+splitMove);

            int categoryIndex = -1;
            int currCategoryCount = boundaryPoints.getDimension();
            int l = pointers.getDimension();


            int[] suitableCategories = new int[boundaryPoints.getDimension()];
            int suitableCategoryCount = 0;
            for(int i = 0; i< suitableCategories.length; i++){
                //System.out.println(boundaryPoints.getParameter(i));
                if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                    suitableCategories[suitableCategoryCount++] = i;
                }

            }

            int[] potentialSplitIndexCounts = new int[suitableCategoryCount];
            int potentialSplitIndexTotal = 0;
            for(int i = 0;i < potentialSplitIndexCounts.length; i++){
                potentialSplitIndexCounts[i] = (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
                potentialSplitIndexTotal+= potentialSplitIndexCounts[i];
            }


            int index = Randomizer.nextInt(potentialSplitIndexTotal)+ 1;

            for(int i = 0; i< potentialSplitIndexCounts.length; i++){
                if(index <= potentialSplitIndexCounts[i]){
                    categoryIndex = suitableCategories[i];
                    break;
                }
                index -= potentialSplitIndexCounts[i];
            }
            if(categoryIndex < 0){
                throw new RuntimeException("Category index error!");
            }


            int splitIndex = (int)boundaryPoints.getValue(categoryIndex,0)+ index;


            int mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);

            int categoryIndex2 = categoryIndex + 1;

            while(mergeIndex == categoryIndex || mergeIndex == categoryIndex2 ){
                mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
            }


            logq += split(
                    categoryIndex,
                    splitIndex,
                    potentialSplitIndexTotal,
                    currCategoryCount
            ) ;
            logq += merge(mergeIndex) ;




        }catch(Exception e){
            throw new RuntimeException(e);
        }


        //System.out.println("logq: "+logq);
        return logq;

    }

    private void temp(){

        int categoryIndex = -1;
        int currCategoryCount = boundaryPoints.getDimension();
        int l = pointers.getDimension();


        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            //System.out.println(boundaryPoints.getParameter(i));
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }

        int[] potentialSplitIndexCounts = new int[suitableCategoryCount];
        int potentialSplitIndexTotal = 0;
        for(int i = 0;i < potentialSplitIndexCounts.length; i++){
            potentialSplitIndexCounts[i] = (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
            potentialSplitIndexTotal+= potentialSplitIndexCounts[i];
        }


        int index = Randomizer.nextInt(potentialSplitIndexTotal)+ 1;

        for(int i = 0; i< potentialSplitIndexCounts.length; i++){
            if(index <= potentialSplitIndexCounts[i]){
                categoryIndex = suitableCategories[i];
                break;
            }
            index -= potentialSplitIndexCounts[i];
        }
        if(categoryIndex < 0){
            throw new RuntimeException("Category index error!");
        }


        int splitIndex = (int)boundaryPoints.getValue(categoryIndex,0)+ index;


        int mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);

        int categoryIndex2 = categoryIndex + 1;

        while(mergeIndex == categoryIndex || mergeIndex == categoryIndex2 ){
            mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
        }

    }

    private double split(
            int categoryIndex,
            int splitIndex,
            int potentialSplitIndexTotal,
            int currCategoryCount) throws Exception{




        double logq = 0.0;
        ParameterList parameterList = parameterListInput.get();
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        logq += proposeNewValue(parameterList, categoryIndex, parametricDistribution);
        //System.out.println("a1: "+logq);
        logq += proposeNewFreqValues(freqsList,categoryIndex);
        //System.out.println("a2: "+logq);
        logq += proposeNewAlphaValue(alphaList, categoryIndex, alphaDistribution);
        //System.out.println("a3: "+logq);
        logq += proposeNewAlphaValue(ratesList, categoryIndex, ratesDistribution);
        //System.out.println("a4: "+logq);
        if(invPrLogit){
            logq += proposeNewValue(invPList, categoryIndex, invPrDistribution);
        }else{
            logq += proposeNewInvPValue(invPList, categoryIndex);
        }

        //System.out.println("a5: "+logq);
        logq += proposeNewDiscreteValue(modelList, categoryIndex, substModelIndicatorDistr);
        //System.out.println("a6: "+logq);

        logq += proposeNewDiscreteValue(siteModelList, categoryIndex, siteModelIndicatorDistr);

        double oldCategoryUpperBound = boundaryPoints.getValue(categoryIndex,0);
        double newCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);

        double prop = ((double)splitIndex - oldCategoryUpperBound+1.0)/((double)newCategoryUpperBound - oldCategoryUpperBound + 1.0);

        boundaryPoints.splitParameter(categoryIndex,1, (double)(splitIndex - 1),categoryIndex+1, new Double[]{(double)splitIndex,newCategoryUpperBound});
        int[] newCateogrySites = new int[(int)newCategoryUpperBound - splitIndex + 1];
        for(int i = 0; i < newCateogrySites.length;i++ ){
            newCateogrySites[i] = i + splitIndex;
        }
        pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(categoryIndex + 1));




        logq += Math.log(potentialSplitIndexTotal);//-Math.log(currCategoryCount);

        return logq;
    }

    private double proposeNewValue(ParameterList parameterList, int categoryIndex, ParametricDistribution distr) throws Exception{
        double[] oldValue = new double[parameterList.getParameterDimension()];
        Double[] newValues = new Double[oldValue.length];
        Double[] sampleVal = distr.sample(1)[0];

        for(int i = 0; i < oldValue.length; i++){
            oldValue[i] = parameterList.getValue(categoryIndex, i);
            newValues[i] = oldValue[i] + sampleVal[i];
        }
        parameterList.splitParameter(categoryIndex, categoryIndex+1, newValues);
        //Jacobian is 0.



        return -distr.calcLogP(new QuietRealParameter(sampleVal));

    }

    private double proposeNewAlphaValue(ParameterList parameterList, int categoryIndex, ParametricDistribution distr) throws Exception{

        /*double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double newValue = oldValue * Math.exp(sampleVal);
        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);


        return -distr.logDensity(sampleVal)+Math.log(newValue);*/

        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double newValue = Math.exp(sampleVal+Math.log(oldValue));
        //Double[] newValues = new Double[]{newValue};
        //Double[] newValue = distr.sample(1)[0];
        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);
        //System.out.println((distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)}))));
        //System.out.println(tmp+" "+(distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)}))-Math.log(tmp)));

        //System.out.println(tmp+" "+oldValue+" "+(distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)-Math.log(oldValue)}))-Math.log(tmp)));

        return -(distr.calcLogP(new QuietRealParameter(new Double[]{sampleVal}))-Math.log(newValue));

    }

    private double proposeNewInvPValue(ParameterList parameterList, int categoryIndex) throws Exception{
        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = invPrDistribution.sample(1)[0][0];
        double oldLogitValue = logit(oldValue);
        double newValue = invLogit(oldLogitValue + sampleVal);
        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);
        double logJacobian = Math.log(newValue)+Math.log(1.0-newValue);

        //System.out.println(-invPrDistribution.logDensity(sampleVal)+logJacobian);

        return -invPrDistribution.logDensity(sampleVal)+logJacobian;
    }

    private double proposeNewFreqValues(ParameterList parameterList, int categoryIndex) throws Exception{
        //System.out.println(parameterList.getParameter(categoryIndex));
        Double[] oldValues = parameterList.getValues(categoryIndex);
        Double[] newValues =  freqsDistribution.nextDirichletScale(oldValues,freqsDistribution.getScaleValue());

        int validCount = 0;
        //while(validCount < 4){
            validCount = 0;
            for(double newVal: newValues){
                if(newVal == 0.0){
                   // newValues =  freqsDistribution.nextDirichletScale(oldValues,freqsDistribution.getScaleValue());
                    /*for(int i = 0; i < newValues.length;i++){
                        System.out.print(newValues[i]+" ");
                    }
                    System.out.println(); */
                    break;
                }else{
                    validCount++;
                }
            }
            if(validCount < 4){
                return Double.NEGATIVE_INFINITY;
            }

        //}

        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValues);
        return -freqsDistribution.logPDF(newValues,oldValues,freqsDistribution.getScaleValue());
    }


    private double proposeNewDiscreteValue(ParameterList parameterList, int categoryIndex, ConditionalCategoricalDistribution conditionalDistr) throws Exception{
        double oldValue = parameterList.getValue(categoryIndex, 0);
        int sampleVal = Randomizer.randomChoicePDF(conditionalDistr.conditionalDensities((int) oldValue));
        double newValue = sampleVal;

        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);


        return -conditionalDistr.logConditionalDensity((int)oldValue,sampleVal);

    }



    private double merge(int mergeIndex) throws Exception{
        //int l = pointers.getDimension();

        //int mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
        //int newCategoryCount = boundaryPoints.getDimension() - 1 ;
        double logq = 0;
        ParameterList parameterList = parameterListInput.get();
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        double temp = boundaryPoints.getValue(mergeIndex+1,1) - boundaryPoints.getValue(mergeIndex,0);
        double temp1 = boundaryPoints.getValue(mergeIndex,0);
        double temp2 = boundaryPoints.getValue(mergeIndex,1);
        double temp3 = boundaryPoints.getValue(mergeIndex+1,0);
        double temp4 = boundaryPoints.getValue(mergeIndex+1,1);
        logq += mergeParameter(parameterList,mergeIndex, parametricDistribution);

        logq += mergeFreqValues(freqsList,mergeIndex);


        logq += mergeAlphaParameter(ratesList,mergeIndex,ratesDistribution);

        logq += mergeAlphaParameter(alphaList,mergeIndex, alphaDistribution);

        if(invPrLogit){
            logq += mergeParameter(invPList,mergeIndex, invPrDistribution);
        }else{
            logq += mergeInvPrParameter(invPList,mergeIndex);
        }

        logq += mergeDiscreteValue(modelList,mergeIndex,substModelIndicatorDistr);
        logq += mergeDiscreteValue(siteModelList,mergeIndex,siteModelIndicatorDistr);

        int removedCatStart = (int)boundaryPoints.getValue(mergeIndex + 1,0);
        int removedCatEnd = (int)boundaryPoints.getValue(mergeIndex + 1,1);
        int[] sitesToMerge = new int[removedCatEnd - removedCatStart + 1];
        for(int i = 0; i < sitesToMerge.length; i++){
            sitesToMerge[i] = removedCatStart + i;
        }

        //System.out.println("sites changed: "+sitesToMerge.length+" "+((double)sitesToMerge.length/(temp)) + " " +temp1 + " "+temp2+" "+temp3+" " +temp4);
        pointers.multiPointerChanges(sitesToMerge,removedCatStart -1);

        double mergedCategoryUpper =  boundaryPoints.getValue(mergeIndex + 1,1);
        //boundaryPoints.setValue(mergeIndex, 1, mergedCategoryUpper);
        //boundaryPoints.removeParameter(mergeIndex + 1);

        boundaryPoints.mergeParameter(mergeIndex+1,mergeIndex,1,mergedCategoryUpper);

        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){

            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }

        int suitableCategoryCountTotal = 0;
        for(int i = 0; i < suitableCategoryCount;i++){
            suitableCategoryCountTotal += (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
        }

        logq += -Math.log(suitableCategoryCountTotal);//+Math.log(newCategoryCount);
        return logq;
    }

    public double mergeParameter(ParameterList parameterList, int mergeIndex, ParametricDistribution distr) throws Exception{
        Double[] oldValues1 = parameterList.getValues(mergeIndex);
        Double[] oldValues2 = parameterList.getValues(mergeIndex + 1);
        Double[] sampleVal = new Double[oldValues1.length];

        for(int i = 0; i < oldValues1.length;i++){
            sampleVal[i] = (oldValues1[i] - oldValues2[i]);
        }

        parameterList.mergeParameter(mergeIndex + 1, mergeIndex);
        //return parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
        return distr.calcLogP(new RealParameter(sampleVal));

    }

    private double mergeAlphaParameter(ParameterList parameterList,int mergeIndex, ParametricDistribution distr) throws Exception{
        /*double oldValue1 = parameterList.getValue(mergeIndex, 0);
        double oldValue2 = parameterList.getValue(mergeIndex+1,0);
        double sample = Math.log(oldValue2) - Math.log(oldValue1);
        parameterList.mergeParameter(mergeIndex + 1,mergeIndex);
        return distr.logDensity(sample)-Math.log(oldValue2); */




        RealParameter oldParameter = parameterList.getParameter(mergeIndex + 1);
        parameterList.mergeParameter(mergeIndex + 1,mergeIndex);
        return distr.calcLogP(new RealParameter(new Double[]{Math.log(oldParameter.getValue()/parameterList.getValue(mergeIndex, 0))}))-Math.log(oldParameter.getValue());
    }

    private double mergeInvPrParameter(ParameterList parameterList, int mergeIndex){
        double oldValue1 = parameterList.getValue(mergeIndex,0);
        double oldValue2 = parameterList.getValue(mergeIndex + 1,0);
        double oldLogitValue1 = logit(oldValue1);
        double oldLogitValue2 = logit(oldValue2);

        double sample = (oldLogitValue2 - oldLogitValue1);
        double logJacobian = -Math.log(oldValue2) - Math.log(1.0 - oldValue2);
        parameterList.mergeParameter(mergeIndex+1,mergeIndex);
        //System.out.println(invPrDistribution.logDensity(sample)+" "+sample);
        return invPrDistribution.logDensity(sample)+logJacobian;
    }

    private double mergeDiscreteValue(ParameterList parameterList, int mergeIndex, ConditionalCategoricalDistribution conditionalDistr) throws Exception{
        double newValue1 = parameterList.getValue(mergeIndex,0);
        double newValue2 = parameterList.getValue(mergeIndex+1,0);

        parameterList.mergeParameter(mergeIndex + 1,mergeIndex);
        return conditionalDistr.logConditionalDensity((int)newValue1,(int)newValue2);


    }

    private double mergeFreqValues(ParameterList parameterList, int mergeIndex) throws Exception{
        Double[] newValue1 = parameterList.getValues(mergeIndex);
        Double[] newValue2 = parameterList.getValues(mergeIndex + 1);

        parameterList.mergeParameter(mergeIndex + 1,mergeIndex);
        return freqsDistribution.logPDF(newValue2,newValue1,freqsDistribution.getScaleValue());


    }

    public double logit(double p){
        return Math.log(p/(1.0-p));
    }
    public double invLogit(double logit){
        return 1.0/(1.0+Math.exp(-logit));
    }
}

