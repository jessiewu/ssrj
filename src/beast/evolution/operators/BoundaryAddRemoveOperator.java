package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.QuietRealParameter;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
public class BoundaryAddRemoveOperator extends Operator {

    public Input<ParameterList> boundaryListInput = new Input<ParameterList>(
            "boundaryList",
            "A list of start and end positions.",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> pointerInput = new Input<DPPointer>("pointers", "A list of pointers that reference to the correct boundaries.", Input.Validate.REQUIRED);

    private boolean autoOptimize;
    private double delta;
    private ParameterList boundaryList;
    private DPPointer pointers;
    public void initAndValidate() {

        boundaryList = boundaryListInput.get();
        pointers = pointerInput.get();



    }

    public double proposal(){
        try{
            double logq = 0.0;
            int dim = boundaryList.getDimension();
            System.out.println(dim);
            if(dim < 3)
                return Double.NEGATIVE_INFINITY;

            int index1 = Randomizer.nextInt(dim - 1);
            int index2 = index1 + 1;

            ArrayList<Integer> potentialNewCat = new ArrayList<Integer>();
            ArrayList<Integer> potentialNewCatSize = new ArrayList<Integer>();
            for(int i = 0; i < dim; i++){
                if(i != index1 && i != index2){
                    int size = (int)(boundaryList.getValue(i,1) - boundaryList.getValue(i,0) + 1);
                    if(size > 1){
                        potentialNewCat.add(i);
                        System.out.println("weird: "+i);
                        potentialNewCatSize.add(size);
                    }

                }
            }

            int potentialNewCatCount = potentialNewCat.size();
            System.out.println("flag1: "+potentialNewCatCount);
            if(potentialNewCatCount == 0){
                return Double.NEGATIVE_INFINITY;
            }





            //Merging
            boundaryList.setValue(index1,1,boundaryList.getValue(index2,1));
            int removedBoundaryStart = (int)boundaryList.getValue(index2,0);
            int removedBoundaryEnd = (int)boundaryList.getValue(index2,1);
            int mergedBoundaryStart = (int)boundaryList.getValue(index1, 0);
            int[] mergeSites = new int[removedBoundaryEnd - removedBoundaryStart + 1];
            for(int i = 0; i < mergeSites.length;i++){
                mergeSites[i] = removedBoundaryStart+i;
            }
            QuietRealParameter mergeBoundary = boundaryList.getParameter(index1);
            QuietRealParameter removedBoundary = boundaryList.getParameter(index2);
            //pointers.multiPointerChanges(mergeSites,(int)boundaryList.getValue(index1,0));

            //Splitting
            int newCatIndexInList = Randomizer.nextInt(potentialNewCatCount);
            int newCatSize = potentialNewCatSize.get(newCatIndexInList);
            double newCatStart = boundaryList.getValue(potentialNewCat.get(newCatIndexInList),0) + Randomizer.nextInt(newCatSize - 1) + 1;
            double newCatEnd = boundaryList.getValue(potentialNewCat.get(newCatIndexInList),1);
            QuietRealParameter newBoundary = new QuietRealParameter(new Double[]{newCatStart,newCatEnd});
            int[] splitSites = new int[(int)(newCatEnd - newCatStart) + 1];
            int newCatStartInt = (int)newCatStart;
            for(int i = 0; i < splitSites.length; i++){
                splitSites[i]= i + newCatStartInt;
            }

            pointers.multiPointerChanges(
                    new int[][]{mergeSites,splitSites},
                    new QuietRealParameter[]{mergeBoundary,newBoundary}
            );

            boundaryList.setValue(potentialNewCat.get(newCatIndexInList),1,newCatStart-1);
            boundaryList.addParameter(potentialNewCat.get(newCatIndexInList)+1,newBoundary);
            boundaryList.removeParameter(removedBoundary);
            logq -= Math.log(1.0/(dim - 1)*1.0/(potentialNewCatCount)*1.0/(newCatSize - 1));


            int newIndex2 =  boundaryList.indexOf(newBoundary);
            int newIndex1 = newIndex2 - 1;
            potentialNewCat.clear();
            for(int i = 0; i < dim; i++){
                if(i != newIndex1 && i != newIndex2){
                    int size = (int)(boundaryList.getValue(i,1) - boundaryList.getValue(i,0) + 1);
                    if(size > 1){
                        potentialNewCat.add(i);

                    }

                }
            }

            for(int i = 0; i< dim; i++){
                System.out.println(boundaryList.getParameter(i));
            }
            System.out.println("------------------");

            for(int i = 0; i < pointers.getDimension();i++){
                System.out.println(pointers.getParameter(i));
            }




            logq += Math.log(1.0/(dim - 1)*1.0/(potentialNewCat.size())*1.0/(removedBoundaryEnd - mergedBoundaryStart));

            return logq;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }
}
