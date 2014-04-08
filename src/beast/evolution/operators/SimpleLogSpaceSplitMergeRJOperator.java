package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 11/06/13
 * Time: 2:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimpleLogSpaceSplitMergeRJOperator extends Operator {
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
            /*if(parameterList.getDimension() == 1){
                logq = split(parameterList);
                //System.out.println("Split");
            }else{
                logq = merge(parameterList);
                //System.out.println("Merge");
            }*/

            boolean splitMove;
            splitMove = Randomizer.nextBoolean();
            ///System.out.println(splitMove+" "+parameterList.getDimension());

            if(!splitMove && parameterList.getDimension() == 1){

                 return Double.NEGATIVE_INFINITY;


            }else if (splitMove &&parameterList.getDimension() == 3){

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

    /*private double split(ParameterList parameterList) throws Exception{
        double oldValue = parameterList.getValue(0, 0);
        double sampleVal = parametricDistribution.sample(1)[0][0];
        double newValue1 = oldValue * Math.exp(sampleVal);
        double newValue2 = oldValue * Math.exp(-sampleVal);
        parameterList.setValue(0,0,newValue1);
        parameterList.addParameter(new QuietRealParameter(newValue2));
        return -parametricDistribution.logDensity(sampleVal)+Math.log(2.0*oldValue);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double oldValue1 = parameterList.getValue(0,0);
        double oldValue2 = parameterList.getValue(1,0);
        double newValue1 = Math.sqrt(oldValue1*oldValue2);
        double sample = Math.log(Math.sqrt(oldValue1/oldValue2));
        parameterList.setValue(0, 0, newValue1);
        parameterList.removeParameter(1);
        //System.out.println(oldValue1+" "+oldValue2);
        return parametricDistribution.logDensity(sample)+Math.log(1.0/2.0/Math.sqrt(oldValue1*oldValue2));
    }*/


    private double split(ParameterList parameterList) throws Exception{
        double oldValue = parameterList.getValue(0, 0);
        double sampleVal = parametricDistribution.sample(1)[0][0];
        double newValue = oldValue * Math.exp(sampleVal);
        parameterList.addParameter(new QuietRealParameter(newValue));
        return -parametricDistribution.logDensity(sampleVal)+Math.log(newValue);

    }

    private double merge(ParameterList parameterList) throws Exception{
        double oldValue1 = parameterList.getValue(0,0);
        double oldValue2 = parameterList.getValue(1,0);
        double sample = Math.log(oldValue2/oldValue1);
        parameterList.removeParameter(1);
        //System.out.println(oldValue1+" "+oldValue2);
        return parametricDistribution.logDensity(sample)-Math.log(oldValue2);
    }
}
