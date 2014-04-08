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
 * Date: 29/06/13
 * Time: 4:32 PM
 * To change this template use File | Settings | File Templates.
 */
public class GIBMAMultiLineSplitMergeRetainRJOperator extends MultiLineSplitMergeRetainRJOperator {

    public GIBMAMultiLineSplitMergeRetainRJOperator(){
        parametricDistributionInput.setRule(Input.Validate.OPTIONAL);
        parameterListInput.setRule(Input.Validate.OPTIONAL);
        modelListInput.setRule(Input.Validate.OPTIONAL);
        substModelIndicatorDistrInput.setRule(Input.Validate.OPTIONAL);
        freqsListInput.setRule(Input.Validate.OPTIONAL);
        freqsDistributionInput.setRule(Input.Validate.OPTIONAL);

    }
    public void initAndValidate(){
        ratesDistribution = ratesDistributionInput.get();
        alphaDistribution = alphaDistributionInput.get();
        invPrDistribution = invPrDistributionInput.get();
        siteModelIndicatorDistr = siteModelIndicatorDistrInput.get();

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
                logq += split() ;
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

    private double split() throws Exception{

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




        double logq = 0.0;
        //ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        //ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        //logq += proposeNewValue(parameterList, categoryIndex, parametricDistribution);
        //System.out.println("a1: "+logq);
        //logq += proposeNewFreqValues(freqsList,categoryIndex);
        //System.out.println("a2: "+logq);
        logq += proposeNewAlphaValue(alphaList, categoryIndex, alphaDistribution);
        //System.out.println("a3: "+logq);
        logq += proposeNewAlphaValue(ratesList, categoryIndex, ratesDistribution);
        //System.out.println("a4: "+logq);
        if(invPrLogit){
            logq += proposeNewValue(invPList, categoryIndex, invPrDistribution);
        }else{
            logq += proposeNewInvPValue(invPList, categoryIndex,invPrDistribution);
        }

        //System.out.println("a5: "+logq);
        //logq += proposeNewDiscreteValue(modelList, categoryIndex, substModelIndicatorDistr);
        //System.out.println("a6: "+logq);

        logq += proposeNewDiscreteValue(siteModelList, categoryIndex, siteModelIndicatorDistr);

        double oldCategoryUpperBound = boundaryPoints.getValue(categoryIndex,0);
        double newCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);

        double prop = ((double)splitIndex - oldCategoryUpperBound+1.0)/((double)newCategoryUpperBound - oldCategoryUpperBound + 1.0);

        //boundaryPoints.setValue(categoryIndex,1,splitIndex - 1);
        //boundaryPoints.addParameter(categoryIndex + 1, new QuietRealParameter(new Double[]{(double)splitIndex,newCategoryUpperBound}));
        boundaryPoints.splitParameter(categoryIndex,1, (double)(splitIndex - 1),categoryIndex+1, new Double[]{(double)splitIndex,newCategoryUpperBound});
        int[] newCateogrySites = new int[(int)newCategoryUpperBound - splitIndex + 1];
        for(int i = 0; i < newCateogrySites.length;i++ ){
            newCateogrySites[i] = i + splitIndex;
        }
        //System.out.println("sites changed: "+((double)newCateogrySites.length/(newCategoryUpperBound - oldCategoryUpperBound)));
        pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(categoryIndex + 1));
        //System.out.println("Sites split: "+newCateogrySites.length);
        //System.out.println(Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount) ;

        //logq += -Math.log(l - currCategoryCount)+Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+ Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);

        logq += Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);
        //System.out.println("a7: "+logq);

        //System.out.println("category: "+categoryIndex);
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

    private double proposeNewInvPValue(ParameterList parameterList, int categoryIndex,ParametricDistribution distr) throws Exception{
        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double oldLogitValue = logit(oldValue);
        double newValue = invLogit(oldLogitValue + sampleVal);
        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);
        double logJacobian = Math.log(newValue)+Math.log(1.0-newValue);

        if(oldValue == 0 && newValue== 0){
            return 0.0;

        }

        return -distr.logDensity(sampleVal)+logJacobian;
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
        double newValue = sampleVal+conditionalDistr.getOffset();

        parameterList.splitParameter(categoryIndex,categoryIndex+1,newValue);


        return -conditionalDistr.logConditionalDensity((int)oldValue,(int)newValue);

    }



    private double merge(ParameterList parameterList) throws Exception{
        int l = pointers.getDimension();

        int mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
        int newCategoryCount = boundaryPoints.getDimension() - 1 ;
        double logq = 0;
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
        //logq += mergeParameter(parameterList,mergeIndex, parametricDistribution);

        //logq += mergeFreqValues(freqsList,mergeIndex);


        logq += mergeAlphaParameter(ratesList,mergeIndex,ratesDistribution);

        logq += mergeAlphaParameter(alphaList,mergeIndex, alphaDistribution);

        if(invPrLogit){
            logq += mergeParameter(invPList,mergeIndex, invPrDistribution);
        }else{
            logq += mergeInvPrParameter(invPList,mergeIndex);
        }

        //logq += mergeDiscreteValue(modelList,mergeIndex,substModelIndicatorDistr);
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
        if(oldValue1 ==0 || oldValue2 == 0){
            return 0.0;
        }
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
