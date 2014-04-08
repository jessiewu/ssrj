package beast.evolution.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ConditionalParametricDistribution;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 9/09/13
 * Time: 10:48 AM
 * To change this template use File | Settings | File Templates.
 */
public class SiteModelRandomWalkOperator extends NetworkIntRandomWalkOperator {

    public Input<List<RealParameter>> associatedParametersInput = new Input<List<RealParameter>>(
            "associatedParameter",
            "Real parameters to be operated on with site model random walk.",
            new ArrayList<RealParameter>()
    );

    public Input<List<ConditionalParametricDistribution>> conditionalParametricDistributionListInput = new Input<List<ConditionalParametricDistribution>>(
            "proposalDistr",
            "A list of proposal distributions on the parameters given the site model proposal.",
            new ArrayList<ConditionalParametricDistribution>()
    );

    public double proposal() {

        double logq = 0.0;
        try{
            List<RealParameter> parameters = associatedParametersInput.get();
            List<ConditionalParametricDistribution> conditionalParametricDistributions = conditionalParametricDistributionListInput.get();

            RealParameter parameter = parameterInput.get(this);
            int currVertex = (int)(parameter.getValue() - offset);
            int nextVertex = neighbours[currVertex][Randomizer.nextInt(neighbours[currVertex].length)];
            parameter.setValue(0,(double)nextVertex+offset);

            //System.out.println("1: "+currVertex+" "+nextVertex +" "+parameter);
            for(int i = 0; i < parameters.size(); i++){
                logq += conditionalParametricDistributions.get(i).calcLogP(parameters.get(i),currVertex);

                double proposalValue = conditionalParametricDistributions.get(i).sample(1,nextVertex)[0][0];
                parameters.get(i).setValue(proposalValue);
                //System.out.println("g: "+parameters.get(i));
                logq -= conditionalParametricDistributions.get(i).calcLogP(parameters.get(i),nextVertex);
            }

            //System.out.println("2: "+currVertex+" "+nextVertex +" "+parameter);
            logq+= logHastingsRatios[currVertex][nextVertex];
            return logq;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }
}
