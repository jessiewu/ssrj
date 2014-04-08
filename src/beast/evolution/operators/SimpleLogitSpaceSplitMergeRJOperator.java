package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;

/**
 * @author Chieh-Hsi Wu
 */
public class SimpleLogitSpaceSplitMergeRJOperator extends Operator {
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

        double oldLogitValue = logit(oldValue);


        double newValue1 = invLogit(oldLogitValue + sampleVal);
        double newValue2 = invLogit(oldLogitValue - sampleVal);
        parameterList.setValue(0,0,newValue1);
        QuietRealParameter newParameter = new QuietRealParameter(newValue2);
        newParameter.setBounds(0.0, 1.0);
        parameterList.addParameter(newParameter);

        double expSampleVal = Math.exp(sampleVal);
        double jnum = 2.0*oldValue*(1 - oldValue)*Math.exp(2.0*sampleVal);
        double jdenumPt1 = (oldValue*expSampleVal + (1 - oldValue));
        double jdenumPt2 = (oldValue + (1 - oldValue)*expSampleVal);
        double jacobian = jnum/(jdenumPt1*jdenumPt1*jdenumPt2*jdenumPt2);

        return -parametricDistribution.logDensity(sampleVal)+Math.log(jacobian);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double oldValue1 = parameterList.getValue(0,0);
        double oldValue2 = parameterList.getValue(1,0);
        double oldLogitValue1 = logit(oldValue1);
        double oldLogitValue2 = logit(oldValue2);

        double newLogitValue = (oldLogitValue1+oldLogitValue2)/2.0;
        double newValue1 = invLogit(newLogitValue);
        double sample = (oldLogitValue1 - oldLogitValue2)/2.0;
        parameterList.setValue(0, 0, newValue1);
        parameterList.removeParameter(1);

        double anum = Math.exp(-newLogitValue);
        double jacobian = anum/((1.0 + anum)*(1.0 + anum))/2*(1.0/oldValue1 + 1.0/(1 - oldValue1))*(1.0/oldValue2 + 1.0/(1 - oldValue2));

        return parametricDistribution.logDensity(sample)+Math.log(jacobian);
    }

    public double logit(double p){
        return Math.log(p/(1.0-p));
    }
    public double invLogit(double logit){
        return 1.0/(1.0+Math.exp(-logit));
    }

}
