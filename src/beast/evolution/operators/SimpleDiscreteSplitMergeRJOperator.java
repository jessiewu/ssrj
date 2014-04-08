package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ConditionalCategoricalDistribution;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
public class SimpleDiscreteSplitMergeRJOperator extends Operator {
    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );

    public Input<ConditionalCategoricalDistribution> conditionalDistrInput = new Input<ConditionalCategoricalDistribution>(
            "conditionalDistr",
            "The distribution to be sampled from given the current value.",
            Input.Validate.REQUIRED
    );



    private ConditionalCategoricalDistribution conditionalDistr;

    public void initAndValidate(){
        conditionalDistr = conditionalDistrInput.get();
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
        int sampleVal = Randomizer.randomChoicePDF(conditionalDistr.conditionalDensities((int) oldValue));
        double newValue = sampleVal;

        if(Randomizer.nextDouble() < 0.5){
            QuietRealParameter newParameter = new QuietRealParameter(newValue);
            newParameter.setBounds(parameterList.getParameter(0).getLower(),parameterList.getParameter(0).getUpper());
            parameterList.addParameter(newParameter);

        }else{
            parameterList.setValue(0,0,newValue);
            QuietRealParameter newParameter = new QuietRealParameter(oldValue);
            newParameter.setBounds(parameterList.getParameter(0).getLower(),parameterList.getParameter(0).getUpper());
            parameterList.addParameter(newParameter);
        }


        return -conditionalDistr.logConditionalDensity((int)oldValue,sampleVal);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double newValue1 = parameterList.getValue(0,0);
        double newValue2 = parameterList.getValue(1,0);

        if(Randomizer.nextDouble() > 0.5){
            parameterList.setValue(0, 0, newValue1);
            parameterList.removeParameter(1);
            return conditionalDistr.logConditionalDensity((int)newValue1,(int)newValue2);
        }else{
            parameterList.setValue(0, 0, newValue2);
            parameterList.removeParameter(1);
            return conditionalDistr.logConditionalDensity((int)newValue2,(int)newValue1);
        }

    }
}
