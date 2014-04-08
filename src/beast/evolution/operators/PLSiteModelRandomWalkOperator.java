package beast.evolution.operators;

import beast.core.Input;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ConditionalParametricDistribution;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
public class PLSiteModelRandomWalkOperator extends PLNetworkIntRandomWalkOperator{
    public Input<List<ParameterList>> associatedParametersInput = new Input<List<ParameterList>>(
            "associatedParameters",
            "Real parameters to be operated on with site model random walk.",
            new ArrayList<ParameterList>()
    );

    public Input<List<ConditionalParametricDistribution>> conditionalParametricDistributionListInput = new Input<List<ConditionalParametricDistribution>>(
            "proposalDistr",
            "A list of proposal distributions on the parameters given the site model proposal.",
            new ArrayList<ConditionalParametricDistribution>()
    );
    public double proposal(){
        try{
            double logq = 0.0;
            List<ParameterList> associatedParameterLists = associatedParametersInput.get();
            List<ConditionalParametricDistribution> conditionalParametricDistributions = conditionalParametricDistributionListInput.get();

            ParameterList paramList = parameterListInput.get(this);
            int iParam = Randomizer.nextInt(paramList.getDimension());
            int currVertex = (int)paramList.getValue(iParam) - offset ;
            int nextVertex = neighbours[currVertex][Randomizer.nextInt(neighbours[currVertex].length)];
            paramList.setValue(iParam,0,nextVertex+ offset);

            for(int i = 0; i < associatedParameterLists.size(); i++){
                logq += conditionalParametricDistributions.get(i).calcLogP(
                        associatedParameterLists.get(i).getParameter(iParam),
                        currVertex
                );


                double proposalValue = conditionalParametricDistributions.get(i).sample(1,nextVertex)[0][0];
                while(proposalValue > associatedParameterLists.get(i).getUpper() ||
                        proposalValue < associatedParameterLists.get(i).getLower()){
                    proposalValue = conditionalParametricDistributions.get(i).sample(1,nextVertex)[0][0];

                }
                //System.out.println(proposalValue+" "+nextVertex);
                associatedParameterLists.get(i).setValue(iParam,0,proposalValue);

                logq -= conditionalParametricDistributions.get(i).calcLogP(
                        associatedParameterLists.get(i).getParameter(iParam),
                        nextVertex
                );
            }

            logq+= logHastingsRatios[currVertex][nextVertex];
            //System.out.println("logq: "+logq);
            return logq;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

}
