package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi
 */
public class SimpleLineSplitMergeRJOperator extends Operator {
    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
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



    public static final double LOG_SPLIT_JACOBIAN = Math.log(2.0);
    public static final double LOG_MERGE_JACOBIAN = Math.log(0.5);

    private ParametricDistribution parametricDistribution;
    private ParameterList boundaryPoints;
    private DPPointer pointers;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
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
            if(!splitMove && parameterList.getDimension() == 1){
                 return Double.NEGATIVE_INFINITY;


            }else if (splitMove && parameterList.getDimension() == pointers.getDimension()){
                return Double.NEGATIVE_INFINITY;

            }

            if(splitMove){
                logq += split(parameterList) ;
                //System.out.println(logq);
                //System.out.println("Split");
            }else{
                logq += merge(parameterList) ;
                //System.out.println("Merge");
            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return logq;

    }

    private double split(ParameterList parameterList) throws Exception{

        //System.out.println(splitIndex);
        int categoryIndex = -1;
        int currCategoryCount = parameterList.getDimension();
        int l = pointers.getDimension();
        /*
        int l = pointers.getDimension();
        int splitIndex = Randomizer.nextInt(pointers.getDimension()-1)+1;

        for(int i = 0; i< boundaryPoints.getDimension(); i++){
                //System.out.println(splitIndex >= boundaryPoints.getValue(i,0));
                //System.out.println(splitIndex <= boundaryPoints.getValue(i,1));
                //System.out.println(splitIndex +" "+ boundaryPoints.getValue(i,1));
                if(splitIndex >= (int)boundaryPoints.getValue(i,0) && splitIndex <= (int)boundaryPoints.getValue(i,1)){
                    categoryIndex = i;
                    break;
                }

        }

        if(categoryIndex == -1){
            throw new RuntimeException("There is something wrong with the boundaries.");
        }

        if(((int)(boundaryPoints.getValue(categoryIndex,1) - boundaryPoints.getValue(categoryIndex,0)) + 1) == 1){
            //System.out.println();
            return Double.NEGATIVE_INFINITY;
        }

        if(splitIndex == boundaryPoints.getValue(categoryIndex,0)){
            //System.out.println();
            return Double.NEGATIVE_INFINITY;
        }
        */



        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            //System.out.println(boundaryPoints.getParameter(i));
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }
        //System.out.println(currCategoryCount+" "+suitableCategoryCount);
        categoryIndex = suitableCategories[Randomizer.nextInt(suitableCategoryCount)];
        //System.out.println( categoryIndex);
        int potentialSplitIndexCount = (int)(boundaryPoints.getValue(categoryIndex,1) - (int)boundaryPoints.getValue(categoryIndex,0));
        int splitIndex = (int)boundaryPoints.getValue(categoryIndex,0)+ Randomizer.nextInt(potentialSplitIndexCount)+1;







        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = parametricDistribution.sample(1)[0][0];
        double newValue1 = oldValue + sampleVal;
        double newValue2 = oldValue - sampleVal;
        parameterList.setValue(categoryIndex,0,newValue1);
        parameterList.addParameter(categoryIndex+1,new QuietRealParameter(newValue2));

        double newCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);
        boundaryPoints.setValue(categoryIndex,1,splitIndex - 1);
        boundaryPoints.addParameter(categoryIndex + 1, new QuietRealParameter(new Double[]{(double)splitIndex,newCategoryUpperBound}));
        //System.out.println(Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN);
        /*
         * The Hasting's ratio:
         * Prior Pr(K + 1)/Pr(K) * Pr(Select a random category)/Pr(select random index) * 1/g(u) + |J_s|
         * K/(L - K) * (L - 1)/K * 1/g(u) + |J|
         * (L - 1)/(L - K) * 1/g(u) + |J|
         * K = curr number of category
         * L = number of units in line
         * return Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;
         * -Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;
         */


        //return Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;
        return Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;
        //System.out.println(currCategoryCount+" "+suitableCategoryCount+" "+potentialSplitIndexCount+boundaryPoints.getParameter(categoryIndex));
        //return -Math.log(l - currCategoryCount)+Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;

    }

    private double merge(ParameterList parameterList) throws Exception{
        int l = pointers.getDimension();
        int mergeIndex = Randomizer.nextInt(parameterList.getDimension() - 1);
        int newCategoryCount = parameterList.getDimension() - 1 ;
        double oldValue1 = parameterList.getValue(mergeIndex,0);
        double oldValue2 = parameterList.getValue(mergeIndex + 1,0);
        double newValue1 = (oldValue1 + oldValue2)/2.0;
        double newValue2 = (oldValue1 - oldValue2)/2.0;


        parameterList.setValue(mergeIndex, 0, newValue1);
        parameterList.removeParameter(mergeIndex + 1);

        double mergedCategoryUpper =  boundaryPoints.getValue(mergeIndex + 1,1);
        boundaryPoints.setValue(mergeIndex, 1, mergedCategoryUpper);
        boundaryPoints.removeParameter(mergeIndex + 1);
        //System.out.println(parameterList);

        /*
         * The Hasting's ratio:
         * Prior Pr(K)/Pr(K - 1) *Pr(select random index)/ Pr(Select a random category) * g(u) + |J_m|
         * K/(L - K) * (L - 1)/K * 1/g(u) + |J|
         * (L - 1)/(L - K) * 1/g(u) + |J|
         * K = curr number of category
         * L = number of units in line
         * Math.log(l -  newCategoryCount) - Math.log(newCategoryCount) + Math.log(newCategoryCount)-Math.log(pointers.getDimension()-1) + parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
         */
        //return Math.log(l -  newCategoryCount)-Math.log(pointers.getDimension()-1) + parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;

        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            //System.out.println(boundaryPoints.getParameter(i));
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }
        double potentialSplitIndexCount = boundaryPoints.getValue(mergeIndex, 1) - boundaryPoints.getValue(mergeIndex, 0);

        //return -Math.log(newCategoryCount)+Math.log(l -  newCategoryCount)-Math.log(suitableCategoryCount) - Math.log(potentialSplitIndexCount)+Math.log(newCategoryCount) + parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
        return -Math.log(suitableCategoryCount) - Math.log(potentialSplitIndexCount)+Math.log(newCategoryCount) + parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
    }
}
