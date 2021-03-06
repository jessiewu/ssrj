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

import java.text.DecimalFormat;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 17/06/13
 * Time: 1:36 AM
 * To change this template use File | Settings | File Templates.
 */
public class MultiLineSplitMergeTestRJOperator extends Operator {
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

    public Input<ConditionalCategoricalDistribution> siteModelIndicatorDistributionInput = new Input<ConditionalCategoricalDistribution>(
            "siteModelIndicatorDistr",
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

    private double windowSize = 0.007;

    protected double logSplitJacobian;
    protected double logMergeJacobian;
    //public static final double LOG_MERGE_JACOBIAN = Math.log(0.5);

    protected ParametricDistribution parametricDistribution;
    protected DirichletDistribution freqsDistribution;
    protected ParametricDistribution alphaDistribution;
    protected ParametricDistribution invPrDistribution;
    protected ParametricDistribution ratesDistribution;
    protected ConditionalCategoricalDistribution siteModelIndicatorDistr;
    protected ConditionalCategoricalDistribution substModelIndicatorDistr;
    protected ParameterList boundaryPoints;
    protected DPPointer pointers;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
        logSplitJacobian = parameterListInput.get().getParameterDimension()*Math.log(2);
        logMergeJacobian = -logSplitJacobian;
        freqsDistribution = freqsDistributionInput.get();
        ratesDistribution = ratesDistributionInput.get();
        alphaDistribution = alphaDistributionInput.get();
        invPrDistribution = invPrDistributionInput.get();
        substModelIndicatorDistr = substModelIndicatorDistrInput.get();
        siteModelIndicatorDistr = siteModelIndicatorDistributionInput.get();
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

        return logq;

    }

    protected double split(ParameterList parameterList) throws Exception{

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
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        logq += proposeNewValue(parameterList, categoryIndex);
        //System.out.println("a1: "+logq);
        logq += proposeNewFreqValues(freqsList,categoryIndex);
        //System.out.println("a2: "+logq);
        logq += proposeNewAlphaValue(alphaList, categoryIndex, alphaDistribution);
        //System.out.println("a3: "+logq);
        logq += proposeNewAlphaValue(ratesList, categoryIndex, ratesDistribution);
        //System.out.println("a4: "+logq);
        logq += proposeNewInvPValue(invPList, categoryIndex);
        //System.out.println("a5: "+logq);
        logq += proposeNewDiscreteValue(modelList, categoryIndex, substModelIndicatorDistr);
        logq += proposeNewDiscreteValue(siteModelList, categoryIndex, siteModelIndicatorDistr);
        //System.out.println("a6: "+logq);
        double newCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);
        //boundaryPoints.setValue(categoryIndex,1,splitIndex - 1);
        //boundaryPoints.addParameter(categoryIndex + 1, new QuietRealParameter(new Double[]{(double)splitIndex,newCategoryUpperBound}));
        boundaryPoints.splitParameter(categoryIndex,1, (double)(splitIndex - 1),categoryIndex+1, new Double[]{(double)splitIndex,newCategoryUpperBound});
        int[] newCateogrySites = new int[(int)newCategoryUpperBound - splitIndex + 1];
        for(int i = 0; i < newCateogrySites.length;i++ ){
            newCateogrySites[i] = i + splitIndex;
        }
        pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(categoryIndex + 1));
        //System.out.println(Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount) ;

        //logq += -Math.log(l - currCategoryCount)+Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+ Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);

        logq += Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);
        //System.out.println("a7: "+logq);

        //System.out.println("category: "+categoryIndex);
        return logq;
    }

    protected double proposeNewValue(ParameterList parameterList, int categoryIndex) throws Exception{
        double[] oldValue = new double[parameterList.getParameterDimension()];
        Double[] newValues1 = new Double[oldValue.length];
        Double[] newValues2 = new Double[oldValue.length];
        Double[] sampleVal = parametricDistribution.sample(1)[0];

        for(int i = 0; i < oldValue.length; i++){
            oldValue[i] = parameterList.getValue(categoryIndex, i);
            newValues1[i] = oldValue[i] + sampleVal[i];
            newValues2[i] = oldValue[i] - sampleVal[i];
        }

        parameterList.splitParameter(categoryIndex, newValues1, categoryIndex+1, newValues2);

        //return -parametricDistribution.calcLogP(new QuietRealParameter(sampleVal))+ LOG_SPLIT_JACOBIAN;

        return -parametricDistribution.calcLogP(new QuietRealParameter(sampleVal))+ logSplitJacobian;

    }

    /*protected double proposeNewAlphaValue(ParameterList parameterList, int categoryIndex, ParametricDistribution distr) throws Exception{

        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = distr.sample(1)[0][0];
        double newValue1 = oldValue * Math.exp(sampleVal);
        double newValue2 = oldValue * Math.exp(-sampleVal);
        parameterList.splitParameter(categoryIndex,newValue1,categoryIndex+1,newValue2);
        //System.out.println(oldValue+" "+sampleVal+" "+newValue1+" "+newValue2);

        return -distr.logDensity(sampleVal)+Math.log(2.0*oldValue);

    } */

    protected double proposeNewAlphaValue(ParameterList parameterList, int categoryIndex, ParametricDistribution distr) throws Exception{

        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = Randomizer.nextGaussian()*windowSize;
        double newValue1 = oldValue * Math.exp(sampleVal);
        double newValue2 = oldValue * Math.exp(-sampleVal);
        parameterList.splitParameter(categoryIndex,newValue1,categoryIndex+1,newValue2);
        //System.out.println(oldValue+" "+sampleVal+" "+newValue1+" "+newValue2);
        double temp = 1.0/(Math.sqrt(2.0*Math.PI*windowSize))*Math.exp(-sampleVal*sampleVal/2.0/windowSize);
        return -Math.log(temp)+Math.log(2.0*oldValue);

    }

    protected double proposeNewInvPValue(ParameterList parameterList, int categoryIndex) throws Exception{
        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = invPrDistribution.sample(1)[0][0];
        double oldLogitValue = logit(oldValue);
        double newValue1 = invLogit(oldLogitValue + sampleVal);
        double newValue2 = invLogit(oldLogitValue - sampleVal);
        parameterList.splitParameter(categoryIndex,newValue1,categoryIndex+1,newValue2);

        double expSampleVal = Math.exp(sampleVal);
        double jnum = 2.0*oldValue*(1 - oldValue)*Math.exp(2.0*sampleVal);
        double jdenumPt1 = (oldValue*expSampleVal + (1 - oldValue));
        double jdenumPt2 = (oldValue + (1 - oldValue)*expSampleVal);
        double jacobian = jnum/(jdenumPt1*jdenumPt1*jdenumPt2*jdenumPt2);
        //System.out.println(oldValue+" "+newValue2);
        return -invPrDistribution.logDensity(sampleVal)+Math.log(jacobian);
    }

    protected double proposeNewFreqValues(ParameterList parameterList, int categoryIndex) throws Exception{
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

        if(Randomizer.nextDouble() < 0.5){
            parameterList.splitParameter(categoryIndex,oldValues,categoryIndex+1,newValues);
        }else{
            parameterList.splitParameter(categoryIndex,newValues,categoryIndex+1,oldValues);
        }
        if(Double.isNaN(freqsDistribution.logPDF(newValues,oldValues,freqsDistribution.getScaleValue()))||
                Double.NEGATIVE_INFINITY == (freqsDistribution.logPDF(newValues,oldValues,freqsDistribution.getScaleValue()))){
            throw new RuntimeException("Crap");
        }


        return -freqsDistribution.logPDF(newValues,oldValues,freqsDistribution.getScaleValue());
    }


    protected double proposeNewDiscreteValue(ParameterList parameterList, int categoryIndex, ConditionalCategoricalDistribution conditionalDistr) throws Exception{
        double oldValue = parameterList.getValue(categoryIndex, 0);
        int sampleVal = Randomizer.randomChoicePDF(conditionalDistr.conditionalDensities((int) oldValue));
        double newValue = sampleVal;

        if(Randomizer.nextDouble() < 0.5){
            parameterList.splitParameter(categoryIndex,oldValue,categoryIndex+1,newValue);
        }else{
            parameterList.splitParameter(categoryIndex,newValue,categoryIndex+1,oldValue);
        }


        return -conditionalDistr.logConditionalDensity((int)oldValue,sampleVal);

    }



    protected double merge(ParameterList parameterList) throws Exception{
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
        logq += mergeParameter(parameterList,mergeIndex);

        logq += mergeFreqValues(freqsList,mergeIndex);

        logq += mergeAlphaParameter(ratesList,mergeIndex, ratesDistribution);

        logq += mergeAlphaParameter(alphaList,mergeIndex, alphaDistribution);
        logq += mergeInvPParameter(invPList,mergeIndex);
        logq += mergeDiscreteValue(modelList,mergeIndex,substModelIndicatorDistr);
        logq += mergeDiscreteValue(siteModelList,mergeIndex,siteModelIndicatorDistr);

        int removedCatStart = (int)boundaryPoints.getValue(mergeIndex + 1,0);
        int removedCatEnd = (int)boundaryPoints.getValue(mergeIndex + 1,1);
        int[] sitesToMerge = new int[removedCatEnd - removedCatStart + 1];
        for(int i = 0; i < sitesToMerge.length; i++){
            sitesToMerge[i] = removedCatStart + i;
        }
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

    protected double mergeParameter(ParameterList parameterList, int mergeIndex) throws Exception{
        Double[] oldValues1 = parameterList.getValues(mergeIndex);
        Double[] oldValues2 = parameterList.getValues(mergeIndex + 1);
        Double[] newValues1 = new Double[oldValues1.length];
        Double[] newValues2 = new Double[oldValues2.length];

        for(int i = 0; i < oldValues1.length;i++){
            newValues1[i] = (oldValues1[i] + oldValues2[i])/2.0;
            newValues2[i] = (oldValues1[i] - oldValues2[i])/2.0;
        }

        parameterList.mergeParameter(mergeIndex + 1, mergeIndex, newValues1);
        //return parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
        return parametricDistribution.calcLogP(new RealParameter(newValues2)) + logMergeJacobian;

    }

    protected double mergeAlphaParameter(ParameterList parameterList,int mergeIndex, ParametricDistribution distr){
        double oldValue1 = parameterList.getValue(mergeIndex,0);
        double oldValue2 = parameterList.getValue(mergeIndex+1,0);
        double newValue1 = Math.sqrt(oldValue1*oldValue2);
        double sample = Math.log(Math.sqrt(oldValue1/oldValue2));
        //System.out.println(oldValue1+" "+oldValue2);
        parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue1);

        double temp = 1.0/(Math.sqrt(2.0*Math.PI*windowSize))*Math.exp(-sample*sample/2.0/windowSize);
        return temp+Math.log(1.0/2.0/Math.sqrt(oldValue1*oldValue2));
    }

    protected double mergeInvPParameter(ParameterList parameterList, int mergeIndex){
        double oldValue1 = parameterList.getValue(mergeIndex,0);
        double oldValue2 = parameterList.getValue(mergeIndex + 1,0);
        double oldLogitValue1 = logit(oldValue1);
        double oldLogitValue2 = logit(oldValue2);

        double newLogitValue = (oldLogitValue1+oldLogitValue2)/2.0;
        double newValue1 = invLogit(newLogitValue);
        double sample = (oldLogitValue1 - oldLogitValue2)/2.0;
        parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue1);

        double anum = Math.exp(-newLogitValue);
        double jacobian = anum/((1.0 + anum)*(1.0 + anum))/2*(1.0/oldValue1 + 1.0/(1 - oldValue1))*(1.0/oldValue2 + 1.0/(1 - oldValue2));

        //System.out.println(invPrDistribution.logDensity(sample)+" "+sample);
        return invPrDistribution.logDensity(sample)+Math.log(jacobian);
    }

    protected double mergeDiscreteValue(ParameterList parameterList, int mergeIndex, ConditionalCategoricalDistribution conditionalDistr) throws Exception{
        double newValue1 = parameterList.getValue(mergeIndex,0);
        double newValue2 = parameterList.getValue(mergeIndex+1,0);

        if(Randomizer.nextDouble() > 0.5){
            parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue1);
            return conditionalDistr.logConditionalDensity((int)newValue1,(int)newValue2);
        }else{
            parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue2);
            return conditionalDistr.logConditionalDensity((int)newValue2,(int)newValue1);
        }

    }

    protected double mergeFreqValues(ParameterList parameterList, int mergeIndex) throws Exception{
        Double[] newValue1 = parameterList.getValues(mergeIndex);
        Double[] newValue2 = parameterList.getValues(mergeIndex + 1);

        if(Randomizer.nextDouble() > 0.5){
            parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue1);
            return freqsDistribution.logPDF(newValue2,newValue1,freqsDistribution.getScaleValue());
        }else{
            parameterList.mergeParameter(mergeIndex + 1,mergeIndex,newValue2);
            return freqsDistribution.logPDF(newValue1,newValue2,freqsDistribution.getScaleValue());
        }

    }

    protected double logit(double p){
        return Math.log(p/(1.0-p));
    }
    protected double invLogit(double logit){
        return 1.0/(1.0+Math.exp(-logit));
    }

    @Override
    public double getCoercableParameterValue() {
        return windowSize;
    }

    @Override
    public void setCoercableParameterValue(double fValue) {
        windowSize = fValue;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double fDelta = calcDelta(logAlpha);

        fDelta += Math.log(windowSize);
        windowSize = Math.exp(fDelta);
        System.out.println(logAlpha+" "+fDelta+" "+windowSize);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = windowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }
    
   
}