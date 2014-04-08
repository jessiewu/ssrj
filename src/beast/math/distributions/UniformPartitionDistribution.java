package beast.math.distributions;

import beast.core.Input;
import beast.core.Valuable;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import org.apache.commons.math.distribution.Distribution;


/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 27/06/13
 * Time: 7:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class UniformPartitionDistribution extends ParametricDistribution {
    public Input<Integer> minCategorySizeInput = new Input<Integer>(
            "minCategorySize",
            "The minimum allowed size of a category in a partition",
            1
    );

    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "The pointers that refers to the boundary list"
    );


    private int minCategorySize;
    private int maxCateogryCount;
    private int partitionLength;
    private double[] logFactorials;

    public void initAndValidate(){
        minCategorySize = minCategorySizeInput.get();
        partitionLength = pointersInput.get().getDimension();
        initiateLogFactorials(partitionLength+1);
        maxCateogryCount = partitionLength/minCategorySize;

    }

    public void initiateLogFactorials(int dim){
        logFactorials = new double[dim];
        logFactorials[0] = 0;
        for(int i = 1; i < dim; i++){
            logFactorials[i] = Math.log(i) + logFactorials[i - 1];

        }

    }

    public Distribution getDistribution(){
        throw new RuntimeException("Not applicable!");
    }

    public double calcLogP(Valuable parameterList){
        ParameterList boundaryList = (ParameterList)parameterList;
        int categoryCount = boundaryList.getDimension();
        if(categoryCount > maxCateogryCount){

            return Double.NEGATIVE_INFINITY;
        }

        int dim = boundaryList.getDimension();
        for(int i = 0; i < dim; i++){
            if(boundaryList.getValue(i,1) - boundaryList.getValue(i,0) + 1 < minCategorySize){
                //System.out.println("WHAT? "+boundaryList.getParameter(i));
                return Double.NEGATIVE_INFINITY;
            }

        }

        double temp = -logChoose(partitionLength - (minCategorySize - 1)*categoryCount-1, categoryCount - 1);
        //System.out.println(categoryCount+" "+Math.exp(temp));
        return temp;
    }

    private double logChoose(int n, int r){
        return logFactorials[n] - logFactorials[r] - logFactorials[n - r];
    }




}
