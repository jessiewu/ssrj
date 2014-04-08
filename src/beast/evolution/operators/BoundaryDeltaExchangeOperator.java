package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.util.Randomizer;

import java.text.DecimalFormat;

/**
 * @author Chieh-Hsi Wu
 */
public class BoundaryDeltaExchangeOperator extends Operator {

    public Input<ParameterList> boundaryListInput = new Input<ParameterList>(
            "boundaryList",
            "A list of start and end positions.",
            Input.Validate.REQUIRED
    );

    public final Input<Boolean> autoOptimizeInput =
            new Input<Boolean>("autoOptimize", "if true, window size will be adjusted during the MCMC run to improve mixing.", true);

    public Input<Double> deltaInput = new Input<Double>("delta", "Magnitude of change for two randomly picked values.", 1.0);
    public Input<DPPointer> pointerInput = new Input<DPPointer>("pointers", "A list of pointers that reference to the correct boundaries.", Input.Validate.REQUIRED);

    private boolean autoOptimize;
    private double delta;
    private ParameterList boundaryList;
    private DPPointer pointers;
    public void initAndValidate() {

        autoOptimize = autoOptimizeInput.get();
        delta = deltaInput.get();
        boundaryList = boundaryListInput.get();
        pointers = pointerInput.get();



    }

    public double proposal(){
        int dim = boundaryList.getDimension();
        if(dim == 1)
            return Double.NEGATIVE_INFINITY;
        int[] boundarySizes = new int[dim];
        for(int i = 0; i < boundarySizes.length;i++){
            boundarySizes[i] = (int)(boundaryList.getValue(i,1) - boundaryList.getValue(i,0))+1;
           // System.out.println(boundarySizes[i]);
        }


        final int dim1 = Randomizer.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = Randomizer.nextInt(dim);

        }

        final int d = Randomizer.nextInt((int) Math.round(delta)) + 1;
        boundarySizes[dim1] = Math.round(boundarySizes[dim1] - d);
        boundarySizes[dim2] = Math.round(boundarySizes[dim2] + d);

        if(boundarySizes[dim1] < 1 || boundarySizes[dim2] < 1)
            return Double.NEGATIVE_INFINITY;

        int startCatIndex = Math.min(dim1,dim2);
        int endCatIndex = Math.max(dim1,dim2);

        int prevStart = 0;
        if(startCatIndex > 0){
            prevStart = (int)boundaryList.getValue(startCatIndex,0);
        }

        for(int i = startCatIndex + 1; i <= endCatIndex; i++ ){
            prevStart += boundarySizes[i - 1];
            boundaryList.setValue(i,0, prevStart);
        }

        int prevEnd = -1;
        if(startCatIndex > 0){
            prevEnd = (int)boundaryList.getValue(startCatIndex - 1,1);
        }
        for(int i = startCatIndex; i < endCatIndex; i++ ){
            prevEnd += boundarySizes[i];
            boundaryList.setValue(i,1, prevEnd);
            //System.out.println(endCatIndex);

        }

        int upperBound;
        int lowerBound;
        int[][] sites = new int[endCatIndex - startCatIndex + 1][];
        //System.out.println(startCatIndex +" "+endCatndex);
        QuietRealParameter[] bounds = new QuietRealParameter[sites.length];
        int l = 0;
        for(int i = startCatIndex; i <= endCatIndex; i++){
            lowerBound = (int)boundaryList.getValue(i,0);
            upperBound = (int)boundaryList.getValue(i,1);
            if(upperBound < lowerBound){
                System.out.println(upperBound+" "+lowerBound);
                for(int j = 0; j< dim; j++){
                    System.out.println(boundaryList.getValue(j,0)+" "+boundaryList.getValue(j,1));
                }
                for(int j = 0; j < dim; j++){
                    System.out.println(boundarySizes[j]);
                }
            }
            sites[l] = new int[upperBound - lowerBound + 1];
            int k = 0;
            for(int j = lowerBound; j <= upperBound;j++ ){
                sites[l][k++] = j;

            }
            bounds[l] = boundaryList.getParameter(i);
            l++;


        }

        pointers.multiPointerChanges(sites,bounds);
        /*for(int i = 0; i < sites.length;i++){
            for(int j = 0; j < sites[i].length; j++){
                System.out.print(sites[i][j]+" "+pointers.indexInList(sites[i][j],boundaryList)+" ");
            }
            System.out.println();
            System.out.println(bounds[i]);
        }  */
        return 0.0;

    }


     @Override
    public double getCoercableParameterValue() {
        return delta;
    }

    @Override
    public void setCoercableParameterValue(final double fValue) {
        delta = fValue;
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
            fDelta += Math.log(delta);
            delta = Math.exp(fDelta);
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
        final double newDelta = delta * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else return "";
    }





}
