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
 * Date: 1/07/13
 * Time: 3:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiLineSplitMergeRetainPropRJOperator extends Operator {
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

    public Input<ParameterList> modelListInput = new Input<ParameterList>(
                "modelList",
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

    private double logSplitJacobian;
    private double logMergeJacobian;
    //public static final double LOG_MERGE_JACOBIAN = Math.log(0.5);

    private ParametricDistribution parametricDistribution;
    private DirichletDistribution freqsDistribution;
    private ParametricDistribution alphaDistribution;
    private ParametricDistribution invPrDistribution;
    private ParametricDistribution ratesDistribution;
    private ConditionalCategoricalDistribution substModelIndicatorDistr;
    private ParameterList boundaryPoints;
    private DPPointer pointers;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
        logSplitJacobian = parameterListInput.get().getParameterDimension()*Math.log(2);
        logMergeJacobian = -logSplitJacobian;
        freqsDistribution = freqsDistributionInput.get();
        ratesDistribution = ratesDistributionInput.get();
        alphaDistribution = alphaDistributionInput.get();
        invPrDistribution = invPrDistributionInput.get();
        substModelIndicatorDistr = substModelIndicatorDistrInput.get();
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

            if(splitMove){
                //System.out.println("Split");
                logq += split(parameterList) ;
                //System.out.println(logq);

            }else{
                //System.out.println("Merge");
                logq += merge(parameterList) ;



            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }


        //System.out.println("logq: "+logq);
        return logq;

    }

    private double split(ParameterList parameterList) throws Exception{

        //System.out.println(splitIndex);
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
        //System.out.println(currCategoryCount+" "+suitableCategoryCount);
        //categoryIndex = suitableCategories[Randomizer.nextInt(suitableCategoryCount)];
        //System.out.println( categoryIndex);
        //int potentialSplitIndexCount = (int)(boundaryPoints.getValue(categoryIndex,1) - (int)boundaryPoints.getValue(categoryIndex,0));
        //int splitIndex = (int)boundaryPoints.getValue(categoryIndex,0)+Randomizer.nextInt(potentialSplitIndexCount)+1;
        //System.out.println(categoryIndex+" "+splitIndex);

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


        int oldCategoryLowerBound = (int)boundaryPoints.getValue(categoryIndex,0);
        double oldCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);

        double prop = ((double)splitIndex - oldCategoryLowerBound+1.0)/((double)oldCategoryUpperBound - oldCategoryLowerBound + 1.0);

        int oldIndex = categoryIndex;
        int newIndex = categoryIndex + 1;
        double logq = prop;

        if(Randomizer.nextDouble() > prop){
            oldIndex = categoryIndex;
            newIndex  = categoryIndex;
            logq = 1.0 - logq;

        }


        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        logq += proposeNewValue(parameterList, oldIndex, newIndex);
        //System.out.println("a1: "+logq);
        logq += proposeNewFreqValues(freqsList,oldIndex, newIndex);
        //System.out.println("a2: "+logq);
        logq += proposeNewAlphaValue(alphaList,  alphaDistribution, oldIndex, newIndex);
        //System.out.println("a3: "+logq);
        logq += proposeNewAlphaValue(ratesList,  ratesDistribution, oldIndex, newIndex);
        //System.out.println("a4: "+logq);
        logq += proposeNewInvPValue(invPList,  oldIndex, newIndex);
        //System.out.println("a5: "+logq);
        logq += proposeNewDiscreteValue(modelList, substModelIndicatorDistr, oldIndex, newIndex);
        //System.out.println("a6: "+logq);

        if(newIndex != oldIndex){

            boundaryPoints.splitParameter(oldIndex, 1, (double)(splitIndex - 1),newIndex, new Double[]{(double)splitIndex,oldCategoryUpperBound});
            int[] newCateogrySites = new int[(int)oldCategoryUpperBound - splitIndex + 1];
            for(int i = 0; i < newCateogrySites.length;i++ ){
                newCateogrySites[i] = i + splitIndex;
            }
            pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(newIndex));
        }else{
            /*System.out.println("Oh?");
            System.out.println(prop);
            System.out.println(splitIndex+" "+oldCategoryLowerBound+" "+oldCategoryUpperBound );   */
            boundaryPoints.splitParameter(newIndex, 0, (double)splitIndex, oldIndex, new Double[]{(double)oldCategoryLowerBound,(double)(splitIndex - 1)});
            int[] newCateogrySites = new int[splitIndex - oldCategoryLowerBound];
            for(int i = 0; i < newCateogrySites.length;i++){
                newCateogrySites[i] = i + oldCategoryLowerBound;
            }
             pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(newIndex));

        }


        logq += Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);

        return logq;
    }

    private double proposeNewValue(ParameterList parameterList, int oldIndex, int newIndex) throws Exception{
        double[] oldValue = new double[parameterList.getParameterDimension()];
        Double[] newValues = new Double[oldValue.length];
        Double[] sampleVal = parametricDistribution.sample(1)[0];

        for(int i = 0; i < oldValue.length; i++){
            oldValue[i] = parameterList.getValue(oldIndex, i);
            newValues[i] = oldValue[i] + sampleVal[i];
        }
        parameterList.splitParameter(oldIndex, newIndex, newValues);
        //Jacobian is 0.

        /*for(int i = 0; i < oldValue.length; i++){
            System.out.println("relativeRates: "+oldValue[i] + " "+newValues[i]);
        }*/

        return -parametricDistribution.calcLogP(new QuietRealParameter(sampleVal));

    }

    private double proposeNewAlphaValue(
            ParameterList parameterList,
            ParametricDistribution distr,
            int oldIndex,
            int newIndex) throws Exception{

        /*double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double newValue = oldValue * Math.exp(sampleVal);
        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);


        return -distr.logDensity(sampleVal)+Math.log(newValue);*/

        double oldValue = parameterList.getValue(oldIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double newValue = Math.exp(sampleVal+Math.log(oldValue));
        //Double[] newValues = new Double[]{newValue};
        //Double[] newValue = distr.sample(1)[0];
        parameterList.splitParameter(oldIndex,newIndex,newValue);
        //System.out.println((distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)}))));
        //System.out.println(tmp+" "+(distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)}))-Math.log(tmp)));

        //System.out.println(tmp+" "+oldValue+" "+(distr.calcLogP(new QuietRealParameter(new Double[]{Math.log(tmp)-Math.log(oldValue)}))-Math.log(tmp)));
        //System.out.println("alpha: "+newValue);

        return -(distr.calcLogP(new QuietRealParameter(new Double[]{sampleVal}))-Math.log(newValue));

    }

    private double proposeNewInvPValue(
            ParameterList parameterList,
            int oldIndex,
            int newIndex) throws Exception{
        double oldValue = parameterList.getValue(oldIndex, 0);
        double sampleVal = invPrDistribution.sample(1)[0][0];
        double oldLogitValue = logit(oldValue);
        double newValue = invLogit(oldLogitValue + sampleVal);
        parameterList.splitParameter(oldIndex,newIndex,newValue);
        double logJacobian = Math.log(newValue)+Math.log(1.0-newValue);

        //System.out.println(-invPrDistribution.logDensity(sampleVal)+logJacobian);

        return -invPrDistribution.logDensity(sampleVal)+logJacobian;
    }

    private double proposeNewFreqValues(ParameterList parameterList, int oldIndex, int newIndex) throws Exception{
        //System.out.println(parameterList.getParameter(categoryIndex));
        Double[] oldValues = parameterList.getValues(oldIndex);
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

        /*for(int i = 0; i < oldValues.length; i++){
            System.out.println( oldValues[i]+ " "+newValues[i]);
        }*/

        parameterList.splitParameter(oldIndex,newIndex,newValues);
        return -freqsDistribution.logPDF(newValues,oldValues,freqsDistribution.getScaleValue());
    }


    private double proposeNewDiscreteValue(
            ParameterList parameterList,
            ConditionalCategoricalDistribution conditionalDistr,
            int oldIndex,
            int newIndex) throws Exception{
        double oldValue = parameterList.getValue(oldIndex, 0);
        int sampleVal = Randomizer.randomChoicePDF(conditionalDistr.conditionalDensities((int) oldValue));
        double newValue = sampleVal;
        parameterList.splitParameter(oldIndex,newIndex,newValue);


        return -conditionalDistr.logConditionalDensity((int)oldValue,sampleVal);

    }



    private double merge(ParameterList parameterList) throws Exception{
        int l = pointers.getDimension();

        int category1Index = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
        int newCategoryCount = boundaryPoints.getDimension() - 1 ;



        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        //double temp = boundaryPoints.getValue(category1Index+1,1) - boundaryPoints.getValue(category1Index,0);
        double category1Lower = boundaryPoints.getValue(category1Index,0);
        double category1Upper = boundaryPoints.getValue(category1Index,1);
        double category2Lower = boundaryPoints.getValue(category1Index+1,0);
        double category2Upper = boundaryPoints.getValue(category1Index+1,1);
        double prop = ((double)(category1Upper) -(category1Lower) + 1.0)/((double)category2Upper - category1Lower + 1.0);
        double logq = prop;
        int mergeToIndex = category1Index;
        int mergeFromIndex = category1Index + 1;

        if(Randomizer.nextDouble() > prop){
            mergeToIndex = category1Index + 1;
            mergeFromIndex = category1Index;
            logq  = 1.0 - prop;
        }



        logq += mergeParameter(parameterList, mergeToIndex, mergeFromIndex);
        logq += mergeFreqValues(freqsList, mergeToIndex, mergeFromIndex);
        logq += mergeAlphaParameter(ratesList, ratesDistribution, mergeToIndex, mergeFromIndex);
        logq += mergeAlphaParameter(alphaList, alphaDistribution, mergeToIndex, mergeFromIndex);
        logq += mergeInvPrParameter(invPList, mergeToIndex, mergeFromIndex);
        logq += mergeDiscreteValue(modelList, substModelIndicatorDistr, mergeToIndex, mergeFromIndex);

            int removedCatStart = (int)boundaryPoints.getValue(mergeFromIndex,0);
            int removedCatEnd = (int)boundaryPoints.getValue(mergeFromIndex,1);
            int[] sitesToMerge = new int[removedCatEnd - removedCatStart + 1];
            for(int i = 0; i < sitesToMerge.length; i++){
                sitesToMerge[i] = removedCatStart + i;
            }
        if(mergeToIndex == category1Index){

            //System.out.println("sites changed: "+sitesToMerge.length+" "+((double)sitesToMerge.length/(temp)) + " " +temp1 + " "+temp2+" "+temp3+" " +temp4);
            pointers.multiPointerChanges(sitesToMerge,removedCatStart -1);
            double mergedCategoryUpper =  boundaryPoints.getValue(mergeFromIndex,1);
            boundaryPoints.mergeParameter(mergeFromIndex,mergeToIndex,1,mergedCategoryUpper);
        }  else {

            //System.out.println("sites changed: "+sitesToMerge.length+" "+((double)sitesToMerge.length/(temp)) + " " +temp1 + " "+temp2+" "+temp3+" " +temp4);
            pointers.multiPointerChanges(sitesToMerge,removedCatEnd + 1);
            double mergedCategoryLower =  boundaryPoints.getValue(mergeFromIndex,0);
            boundaryPoints.mergeParameter(mergeFromIndex,mergeToIndex,0,mergedCategoryLower);

        }


        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            //System.out.println(boundaryPoints.getParameter(i));
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }
        //double potentialSplitIndexCount = boundaryPoints.getValue(mergeIndex, 1) - boundaryPoints.getValue(mergeIndex, 0);

        int suitableCategoryCountTotal = 0;
        for(int i = 0; i < suitableCategoryCount;i++){
            suitableCategoryCountTotal += (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
        }
        //System.out.println(parameterList);
        //logq +=Math.log(l -  newCategoryCount) - Math.log(newCategoryCount) + Math.log(newCategoryCount)-Math.log(pointers.getDimension()-1) ;


        //logq += Math.log(l -  newCategoryCount)-Math.log(suitableCategoryCount) - Math.log(potentialSplitIndexCount);
        //logq += - Math.log(newCategoryCount)+Math.log(l -  newCategoryCount) - Math.log(suitableCategoryCountTotal)+Math.log(newCategoryCount);
        logq += -Math.log(suitableCategoryCountTotal)+Math.log(newCategoryCount);
        return logq;
    }

    public double mergeParameter(
            ParameterList parameterList,
            int mergeToIndex,
            int mergeFromIndex) throws Exception{
        Double[] oldValues1 = parameterList.getValues(mergeToIndex);
        Double[] oldValues2 = parameterList.getValues(mergeFromIndex);
        Double[] sampleVal = new Double[oldValues1.length];

        for(int i = 0; i < oldValues1.length;i++){
            sampleVal[i] = (oldValues1[i] - oldValues2[i]);
        }

        parameterList.mergeParameter(mergeFromIndex, mergeToIndex);
        //return parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
        return parametricDistribution.calcLogP(new RealParameter(sampleVal));

    }

    private double mergeAlphaParameter(
            ParameterList parameterList,
            ParametricDistribution distr,
            int mergeToIndex,
            int mergeFromIndex) throws Exception{

        RealParameter oldParameter = parameterList.getParameter(mergeFromIndex);
        double temp = Math.log(oldParameter.getValue()/parameterList.getValue(mergeToIndex, 0));
        parameterList.mergeParameter(mergeFromIndex, mergeToIndex);

        return distr.calcLogP(new RealParameter(new Double[]{temp}))-Math.log(oldParameter.getValue());
    }

    private double mergeInvPrParameter(
            ParameterList parameterList,
            int mergeToIndex,
            int mergeFromIndex){
        double oldValue1 = parameterList.getValue(mergeToIndex,0);
        double oldValue2 = parameterList.getValue(mergeFromIndex,0);
        double oldLogitValue1 = logit(oldValue1);
        double oldLogitValue2 = logit(oldValue2);

        double sample = (oldLogitValue2 - oldLogitValue1);
        double logJacobian = -Math.log(oldValue2) - Math.log(1.0 - oldValue2);
        parameterList.mergeParameter(mergeFromIndex,mergeToIndex);
        //System.out.println(invPrDistribution.logDensity(sample)+" "+sample);
        return invPrDistribution.logDensity(sample)+logJacobian;
    }

    private double mergeDiscreteValue(
            ParameterList parameterList,
            ConditionalCategoricalDistribution conditionalDistr,
            int mergeToIndex,
            int mergeFromIndex) throws Exception{
        double newValue1 = parameterList.getValue(mergeToIndex,0);
        double newValue2 = parameterList.getValue(mergeFromIndex,0);

        parameterList.mergeParameter(mergeFromIndex,mergeToIndex);
        return conditionalDistr.logConditionalDensity((int)newValue1,(int)newValue2);


    }

    private double mergeFreqValues(
            ParameterList parameterList,
            int mergeToIndex,
            int mergeFromIndex) throws Exception{
        Double[] newValue1 = parameterList.getValues(mergeToIndex);
        Double[] newValue2 = parameterList.getValues(mergeFromIndex);

        parameterList.mergeParameter(mergeFromIndex, mergeToIndex);
        return freqsDistribution.logPDF(newValue2,newValue1,freqsDistribution.getScaleValue());


    }

    public double logit(double p){
        return Math.log(p/(1.0-p));
    }
    public double invLogit(double logit){
        return 1.0/(1.0+Math.exp(-logit));
    }

}
