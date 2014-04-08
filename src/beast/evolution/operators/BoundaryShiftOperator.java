package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
public class BoundaryShiftOperator extends Operator {
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

    public void initAndValidate(){}




    public double proposal(){
        double logq = 0;
        ParameterList boundaryList = boundaryPointsInput.get();
        DPPointer pointers = pointersInput.get();
        if(boundaryList.getDimension() == 1 || pointers.getDimension() == boundaryList.getDimension()){
            return Double.NEGATIVE_INFINITY;
        }
        int lastSplitPointIndex = boundaryList.getDimension() - 1;
        int[] potentialMovablePoints = new int[lastSplitPointIndex];
        int k = 0;
        for(int i = 0; i < lastSplitPointIndex; i++){
            if((boundaryList.getValue(i+1,1) - boundaryList.getValue(i,0))> 1){
                potentialMovablePoints[k++] = i;
            }
        }


        int movePointIndex = potentialMovablePoints[Randomizer.nextInt(k)];
        int moveToLeftCount = (int)(boundaryList.getValue(movePointIndex,1) - boundaryList.getValue(movePointIndex,0));
        int moveToRightCount = (int)(boundaryList.getValue(movePointIndex+1,1) - boundaryList.getValue(movePointIndex+1,0));
        int potentialPosCount =  moveToLeftCount + moveToRightCount;

        int nextIndex = Randomizer.nextInt(potentialPosCount);
        int nextPos;
        pointers = pointersInput.get(this);

        if(nextIndex < moveToLeftCount){
            nextPos = (int)boundaryList.getValue(movePointIndex,0) + nextIndex;
            int temp = (int)boundaryList.getValue(movePointIndex+1,1);
            int[] sitesAffected = new int[temp - nextPos];
            for(int i = 0; i < sitesAffected.length; i++){
                sitesAffected[i] = nextPos + 1 + i;

            }
            pointers.multiPointerChanges(sitesAffected,temp);
        }else{
            nextPos = (int)boundaryList.getValue(movePointIndex+1,0) + (nextIndex - moveToLeftCount);
            int temp = (int)boundaryList.getValue(movePointIndex,1);
            int[] sitesAffected = new int[nextPos-temp];
            for(int i = 0; i < sitesAffected.length; i++){
                sitesAffected[i] = temp+i+1;

            }
            pointers.multiPointerChanges(sitesAffected,temp);


        }
        //System.out.println(boundaryList.getParameter(movePointIndex));
        //System.out.println(boundaryList.getParameter(movePointIndex+1));
        boundaryList = boundaryPointsInput.get(this);
        boundaryList.setValue(movePointIndex,1,nextPos);
        boundaryList.setValue(movePointIndex+1,0,nextPos+1);
        //System.out.println(boundaryList.getParameter(movePointIndex));
        //System.out.println(boundaryList.getParameter(movePointIndex+1));

        int newMoveToLeftCount = (int)(boundaryList.getValue(movePointIndex,1) - boundaryList.getValue(movePointIndex,0));
        int newMoveToRightCount = (int)(boundaryList.getValue(movePointIndex+1,1) - boundaryList.getValue(movePointIndex+1,0));

        //System.out.println("movePointIndex: "+movePointIndex);
        /*int newK = k;
        if(moveToLeftCount == 0 && newMoveToLeftCount > 0){
            newK++;
        }else if(moveToLeftCount > 0 && newMoveToLeftCount == 0  && movePointIndex > 0){
            newK--;

        }

        if(moveToRightCount == 0 && newMoveToRightCount > 0){
            newK++;
        }else if(moveToRightCount > 0 && newMoveToRightCount == 0 && movePointIndex + 1 < boundaryList.getDimension() - 1 ){
            //System.out.println("Hi?");
            newK--;
        }

        //(1/newK*1/potentialPosCount)/(1/k*1/potentialPosCount);
        if(newK == 0){
            throw new RuntimeException(k+" "+newK+" "+newMoveToLeftCount+" "+newMoveToRightCount);
        }  */

        int newK = 0;

        for(int i = 0; i < lastSplitPointIndex; i++){
            if((boundaryList.getValue(i+1,1) - boundaryList.getValue(i,0))> 1){
                newK ++;
            }
        }

        /*if(newK == 0){
            throw new RuntimeException(k+" "+newK+" "+newMoveToLeftCount+" "+newMoveToRightCount);
        }

        System.out.println(k + " "+newK);  */
        //System.out.println(((double)k/(double)newK));
        return ((double)k/(double)newK);
    }
}
