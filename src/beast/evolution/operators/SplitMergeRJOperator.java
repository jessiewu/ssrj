package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 14/06/13
 * Time: 6:30 PM
 * To change this template use File | Settings | File Templates.
 */
public class SplitMergeRJOperator extends Operator {
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
            //System.out.println(splitMove+" "+parameterList.getDimension());
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

        int currCategoryCount = parameterList.getDimension();
        int categoryIndex = Randomizer.nextInt(currCategoryCount);

        double oldValue = parameterList.getValue(categoryIndex, 0);
        double sampleVal = parametricDistribution.sample(1)[0][0];
        double newValue1 = oldValue + sampleVal;
        double newValue2 = oldValue - sampleVal;
        parameterList.setValue(categoryIndex,0,newValue1);
        parameterList.addParameter(categoryIndex+1,new QuietRealParameter(newValue2));

       // return Math.log(currCategoryCount)+Math.log(1.0/2.0*currCategoryCount*(currCategoryCount+1))-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;
         return -parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN;

    }

    private double merge(ParameterList parameterList) throws Exception{
        int currCategoryCount = parameterList.getDimension();
        int categoryIndex1 = Randomizer.nextInt(currCategoryCount);
        int categoryIndex2 = Randomizer.nextInt(currCategoryCount);
        while(categoryIndex1 == categoryIndex2){
            categoryIndex2 = Randomizer.nextInt(currCategoryCount);
        }

        int newCategoryCount = parameterList.getDimension() - 1 ;
        double oldValue1 = parameterList.getValue(categoryIndex1,0);
        double oldValue2 = parameterList.getValue(categoryIndex2,0);
        double newValue1 = (oldValue1 + oldValue2)/2.0;
        double newValue2 = (oldValue1 - oldValue2)/2.0;


        parameterList.setValue(categoryIndex1, 0, newValue1);
        parameterList.removeParameter(categoryIndex2);


        //System.out.println(parameterList);
        //System.out.println(-Math.log(currCategoryCount - 1) + Math.log(1.0 / 2.0 * currCategoryCount * (currCategoryCount - 1)));
        //return -Math.log(currCategoryCount - 1)+Math.log(1.0/2.0*currCategoryCount*(currCategoryCount - 1)) + parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
        return parametricDistribution.logDensity(newValue2) + LOG_MERGE_JACOBIAN;
    }
}
