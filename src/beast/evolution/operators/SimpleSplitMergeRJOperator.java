package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;


/**
 * @author Chieh-Hsi Wu
 */
public class SimpleSplitMergeRJOperator extends Operator {
    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
            "An ParameterList object that contains a list of parameters.",
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

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
    }

    public double proposal() {
        double logq = 0.0;
        try{
            ParameterList parameterList = parameterListInput.get(this);
            if(parameterList.getDimension() == 1){
                logq = split(parameterList) + LOG_SPLIT_JACOBIAN;
                //System.out.println("Split");
            }else{
                logq = merge(parameterList) + LOG_MERGE_JACOBIAN;
                //System.out.println("Merge");
            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return logq;

    }

    private double split(ParameterList parameterList) throws Exception{
        double oldValue = parameterList.getValue(0, 0);
        double sampleVal = parametricDistribution.sample(1)[0][0];
        double newValue1 = oldValue + sampleVal;
        double newValue2 = oldValue - sampleVal;
        parameterList.setValue(0,0,newValue1);
        parameterList.addParameter(new QuietRealParameter(newValue2));
        return -parametricDistribution.logDensity(sampleVal);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double oldValue1 = parameterList.getValue(0,0);
        double oldValue2 = parameterList.getValue(1,0);
        double newValue1 = (oldValue1 + oldValue2)/2.0;
        double newValue2 = (oldValue1 - oldValue2)/2.0;
        parameterList.setValue(0, 0, newValue1);
        parameterList.removeParameter(1);
        return parametricDistribution.logDensity(newValue2);
    }
}
