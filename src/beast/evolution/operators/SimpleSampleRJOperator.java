package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
public class SimpleSampleRJOperator extends Operator {
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



    private ParametricDistribution parametricDistribution;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
    }

    public double proposal() {
        double logq = 0.0;
        try{
            ParameterList parameterList = parameterListInput.get(this);
            if(parameterList.getDimension() == 1){
                logq = split(parameterList);
                //System.out.println("Split");
            }else{
                logq = merge(parameterList);
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
        double newValue = sampleVal;
        if(Randomizer.nextDouble() < 0.5){
            parameterList.addParameter(new QuietRealParameter(newValue));

        }else{
            parameterList.setValue(0,0,newValue);
            parameterList.addParameter(new QuietRealParameter(oldValue));
        }


        return -parametricDistribution.logDensity(sampleVal);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double newValue1 = parameterList.getValue(0,0);
        double newValue2 = parameterList.getValue(1,0);

        if(Randomizer.nextDouble() > 0.5){
            parameterList.setValue(0, 0, newValue1);
            parameterList.removeParameter(1);
            return parametricDistribution.logDensity(newValue2);
        }else{
            parameterList.setValue(0, 0, newValue2);
            parameterList.removeParameter(1);
            return parametricDistribution.logDensity(newValue1);
        }

    }
}
