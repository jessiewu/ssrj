package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.util.Randomizer;

import java.text.DecimalFormat;

/**
 * @author
 */
public class BoundaryRandomWalkOperator extends Operator {
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

    public Input<Double> windowSizeInput = new Input<Double>(
            "windowSize",
            "The maximum step size of the random walk on the boundary.",
            Input.Validate.REQUIRED
    );

    public final Input<Boolean> autoOptimizeInput =
            new Input<Boolean>(
                    "autoOptimize",
                    "if true, window size will be adjusted during the MCMC run to improve mixing.",
                    true
            );


    private double windowSize;
    private boolean autoOptimize;
    public void initAndValidate(){
        windowSize = windowSizeInput.get();
        autoOptimize = autoOptimizeInput.get();
    }

    public double proposal(){
        double logq = 0;
       ParameterList boundaryList = boundaryPointsInput.get();
        DPPointer pointers = pointersInput.get();
        if(boundaryList.getDimension() == 1 || pointers.getDimension() == boundaryList.getDimension()){
            return Double.NEGATIVE_INFINITY;
        }

         int lastSplitPointIndex = boundaryList.getDimension() - 1;
        /*int[] potentialMovablePoints = new int[lastSplitPointIndex];
        int k = 0;
        for(int i = 0; i < lastSplitPointIndex; i++){
            if((boundaryList.getValue(i+1,1) - boundaryList.getValue(i,0))> 1){
                potentialMovablePoints[k++] = i;
            }
        }


        int movePointIndex = potentialMovablePoints[Randomizer.nextInt(k)];*/
        int movePointIndex = Randomizer.nextInt(lastSplitPointIndex);
        int lowerBound = (int) boundaryList.getValue(movePointIndex,0);
        int upperBound = (int) boundaryList.getValue(movePointIndex+1,1) - 1;
        int step = Randomizer.nextInt((int) Math.ceil(windowSize))+1;
        step = Randomizer.nextBoolean()? step:-step;
        //System.out.println(step);
        int nextPos = (int)boundaryList.getValue(movePointIndex,1)+step;
        if(nextPos < lowerBound || nextPos > upperBound){
            return Double.NEGATIVE_INFINITY;
        }
        if(step != 0){

            if(step < 0){
                //nextPos = (int)boundaryList.getValue(movePointIndex,1) - step;
                int temp = (int)boundaryList.getValue(movePointIndex,1);
                int[] sitesAffected = new int[temp - nextPos];
                for(int i = 0; i < sitesAffected.length; i++){
                    sitesAffected[i] = nextPos + 1 + i;
                }
                pointers.multiPointerChanges(sitesAffected,temp+1);

            }else{
                //nextPos = (int)boundaryList.getValue(movePointIndex,1)+ step;
                int temp = (int)boundaryList.getValue(movePointIndex,1);
                int[] sitesAffected = new int[nextPos-temp];
                for(int i = 0; i < sitesAffected.length; i++){
                    sitesAffected[i] = temp+i + 1;

                }
                pointers.multiPointerChanges(sitesAffected,temp);

            }

            boundaryList = boundaryPointsInput.get(this);
            boundaryList.setValue(movePointIndex,1,nextPos);
            boundaryList.setValue(movePointIndex+1,0,nextPos+1);
        }


        /*int moveToLeftCount = (int)(boundaryList.getValue(movePointIndex,1) - boundaryList.getValue(movePointIndex,0));
        int moveToRightCount = (int)(boundaryList.getValue(movePointIndex+1,1) - boundaryList.getValue(movePointIndex+1,0));
        int potentialPosCount =  moveToLeftCount + moveToRightCount;   */


        /*int step = Randomizer.nextInt((int)Math.ceil(windowSize));
        int potentialDirections = 0;
        int nextPos;
        boolean moveToLeft;
        if(moveToLeftCount >= step && moveToRightCount >= step){
            potentialDirections = 2;
            moveToLeft = Randomizer.nextBoolean();

        }else if(moveToLeftCount >= step){
            potentialDirections = 1;
            moveToLeft = true;

        }else if(moveToRightCount >= step){
            potentialDirections = 1;
            moveToLeft = false;
        }else{
            return Double.NEGATIVE_INFINITY;
        }


        int step = Randomizer.nextInt((int)Math.ceil(windowSize))+1;
        boolean moveToLeft = Randomizer.nextBoolean();
        int nextPos;
        if(moveToLeft && (((int)boundaryList.getValue(movePointIndex,1) - step) < (int)boundaryList.getValue(movePointIndex,0))){
            return Double.NEGATIVE_INFINITY;
        }else if((((int)boundaryList.getValue(movePointIndex,1) + step) >= (int)boundaryList.getValue(movePointIndex+1,1))){
            return Double.NEGATIVE_INFINITY;
        }

        if(moveToLeft){nextPos = (int)boundaryList.getValue(movePointIndex,1) - step;
            int temp = (int)boundaryList.getValue(movePointIndex,1);
            int[] sitesAffected = new int[temp - nextPos];
            for(int i = 0; i < sitesAffected.length; i++){
                sitesAffected[i] = nextPos + 1 + i;
            }
            pointers.multiPointerChanges(sitesAffected,temp+1);

        }else{
            nextPos = (int)boundaryList.getValue(movePointIndex,1)+ step;
            int temp = (int)boundaryList.getValue(movePointIndex,1);
            int[] sitesAffected = new int[nextPos-temp];
            for(int i = 0; i < sitesAffected.length; i++){
                sitesAffected[i] = temp+i + 1;

            }
            pointers.multiPointerChanges(sitesAffected,temp);

        }





        boundaryList = boundaryPointsInput.get(this);
        boundaryList.setValue(movePointIndex,1,nextPos);
        boundaryList.setValue(movePointIndex+1,0,nextPos+1);       */

        /*int newMoveToLeftCount = (int)(boundaryList.getValue(movePointIndex,1) - boundaryList.getValue(movePointIndex,0)) - 1;
        int newMoveToRightCount = (int)(boundaryList.getValue(movePointIndex+1,1) - boundaryList.getValue(movePointIndex+1,0)) - 1;

        int newPotentialDirections = 0;
        if(moveToLeftCount >= step && moveToRightCount >= step){
            newPotentialDirections = 2;

        }else if(moveToLeftCount >= step){
            newPotentialDirections = 1;

        }else if(moveToRightCount >= step){
            newPotentialDirections = 1;
        }    */


        //Pr(step)/Pr(step) * Pr(newPotentialDirections)/Pr(potentialDirections)
        //return ((double)potentialDirections/(double)newPotentialDirections);
        return 0.0;
    }
    @Override
    public double getCoercableParameterValue() {
        return windowSize;
    }

    @Override
    public void setCoercableParameterValue(final double fValue) {
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
    public void optimize(final double logAlpha) {
        // must be overridden by operator implementation to have an effect
        if (autoOptimize) {
            double fDelta = calcDelta(logAlpha);
            fDelta += Math.log(windowSize);
            windowSize = Math.exp(fDelta);
        }

    }

    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double newDelta = windowSize * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else return "";
    }
}
