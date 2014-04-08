package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.distributions.DirichletDistribution;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 2/07/13
 * Time: 3:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiLineSplitMergeRetainSampleRJOperator extends Operator {
        public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );



    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "An ParameterList object that contains a list of parameters.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> alphaListInput = new Input<ParameterList>(
                "alphaList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
    );

    public Input<ParameterList> invPListInput = new Input<ParameterList>(
                "invPrList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );

    public Input<ParameterList> modelListInput = new Input<ParameterList>(
                "modelList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );


    public Input<ParameterList> siteModelListInput = new Input<ParameterList>(
                "gammaSiteModelIndicatorList",
                "An ParameterList object that contains a list of parameters.",
                Input.Validate.REQUIRED
        );

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

    public Input<ParametricDistribution> parametricDistributionInput = new Input<ParametricDistribution>(
            "distr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<DirichletDistribution> freqsDistributionInput = new Input<DirichletDistribution>(
            "freqsDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> ratesDistributionInput = new Input<ParametricDistribution>(
            "rateDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> alphaDistributionInput = new Input<ParametricDistribution>(
            "alphaDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> invPrDistributionInput = new Input<ParametricDistribution>(
            "invPrDistr",
            "The distribution to be sampled from.",
            Input.Validate.REQUIRED
    );



    public Input<ParametricDistribution> substModelIndicatorDistrInput = new Input<ParametricDistribution>(
            "substModelIndicatorDistr",
            "The distribution to be sampled from given the current value.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> siteModelIndicatorDistrInput = new Input<ParametricDistribution>(
            "siteModelIndicatorDistr",
            "The distribution to be sampled from given the current value.",
            Input.Validate.REQUIRED
    );




    //public static final double LOG_MERGE_JACOBIAN = Math.log(0.5);

    private ParametricDistribution parametricDistribution;
    private DirichletDistribution freqsDistribution;
    private ParametricDistribution alphaDistribution;
    private ParametricDistribution invPrDistribution;
    private ParametricDistribution ratesDistribution;
    private ParametricDistribution substModelIndicatorDistr;
    private ParametricDistribution siteModelIndicatorDistr;
    private ParameterList boundaryPoints;
    private DPPointer pointers;

    public void initAndValidate(){
        parametricDistribution = parametricDistributionInput.get();
        freqsDistribution = freqsDistributionInput.get();
        ratesDistribution = ratesDistributionInput.get();
        alphaDistribution = alphaDistributionInput.get();
        invPrDistribution = invPrDistributionInput.get();
        substModelIndicatorDistr = substModelIndicatorDistrInput.get();
        siteModelIndicatorDistr = siteModelIndicatorDistrInput.get();
    }

    public double proposal() {
        double logq = 0.0;


        try{
            ParameterList parameterList = parameterListInput.get();
            pointers = pointersInput.get();
            boundaryPoints = boundaryPointsInput.get();
            boolean splitMove;
            splitMove = Randomizer.nextBoolean();
            ///System.out.println(splitMove+" "+parameterList.getDimension());

            if(!splitMove && boundaryPoints.getDimension() == 1){

                 return Double.NEGATIVE_INFINITY;


            }else if (splitMove && boundaryPoints.getDimension() == pointers.getDimension()){

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


        //System.out.println("logq: "+logq);
        return logq;

    }

    private double split(ParameterList parameterList) throws Exception{

        //System.out.println(splitIndex);
        int categoryIndex = -1;
        int currCategoryCount = boundaryPoints.getDimension();
        int l = pointers.getDimension();


        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }

        int[] potentialSplitIndexCounts = new int[suitableCategoryCount];
        int potentialSplitIndexTotal = 0;
        for(int i = 0;i < potentialSplitIndexCounts.length; i++){
            potentialSplitIndexCounts[i] = (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
            potentialSplitIndexTotal+= potentialSplitIndexCounts[i];
        }


        int index = Randomizer.nextInt(potentialSplitIndexTotal)+ 1;

        for(int i = 0; i< potentialSplitIndexCounts.length; i++){
            if(index <= potentialSplitIndexCounts[i]){
                categoryIndex = suitableCategories[i];
                break;
            }
            index -= potentialSplitIndexCounts[i];
        }
        if(categoryIndex < 0){
            throw new RuntimeException("Category index error!");
        }


        int splitIndex = (int)boundaryPoints.getValue(categoryIndex,0)+ index;




        double logq = 0.0;
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        logq += proposeNewParameter(parameterList, categoryIndex, parametricDistribution);
        //System.out.println("a1: "+logq);
        logq += proposeNewParameter(freqsList,categoryIndex, freqsDistribution);
        //System.out.println("a2: "+logq);
        logq += proposeNewParameter(alphaList, categoryIndex, alphaDistribution);
        //System.out.println("a3: "+logq);
        logq += proposeNewParameter(ratesList, categoryIndex, ratesDistribution);
        //System.out.println("a4: "+logq);
        logq += proposeNewParameter(invPList, categoryIndex, invPrDistribution);
        //System.out.println("a5: "+logq);
        logq += proposeNewParameter(modelList, categoryIndex, substModelIndicatorDistr);
        //System.out.println("a6: "+logq);
        logq += proposeNewParameter(siteModelList, categoryIndex, siteModelIndicatorDistr);

        double oldCategoryUpperBound = boundaryPoints.getValue(categoryIndex,0);
        double newCategoryUpperBound = boundaryPoints.getValue(categoryIndex,1);

        double prop = ((double)splitIndex - oldCategoryUpperBound+1.0)/((double)newCategoryUpperBound - oldCategoryUpperBound + 1.0);

        //boundaryPoints.setValue(categoryIndex,1,splitIndex - 1);
        //boundaryPoints.addParameter(categoryIndex + 1, new QuietRealParameter(new Double[]{(double)splitIndex,newCategoryUpperBound}));
        boundaryPoints.splitParameter(categoryIndex,1, (double)(splitIndex - 1),categoryIndex+1, new Double[]{(double)splitIndex,newCategoryUpperBound});
        int[] newCateogrySites = new int[(int)newCategoryUpperBound - splitIndex + 1];
        for(int i = 0; i < newCateogrySites.length;i++ ){
            newCateogrySites[i] = i + splitIndex;
        }
        //System.out.println("sites changed: "+((double)newCateogrySites.length/(newCategoryUpperBound - oldCategoryUpperBound)));
        pointers.multiPointerChanges(newCateogrySites, boundaryPoints.getParameter(categoryIndex + 1));
        //System.out.println(Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount)-parametricDistribution.logDensity(sampleVal) + LOG_SPLIT_JACOBIAN);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+Math.log(pointers.getDimension()-1)-Math.log(currCategoryCount) ;

        //logq += -Math.log(l - currCategoryCount)+Math.log(suitableCategoryCount)+Math.log(potentialSplitIndexCount);
        //logq += Math.log(currCategoryCount)-Math.log(l - currCategoryCount)+ Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);

        logq += Math.log(potentialSplitIndexTotal)-Math.log(currCategoryCount);
        //System.out.println("a7: "+logq);

        //System.out.println("category: "+categoryIndex);
        return logq;
    }






    private double proposeNewParameter(ParameterList parameterList, int categoryIndex, ParametricDistribution distr) throws Exception{
        Double[] oldValues = parameterList.getValues(categoryIndex);
        Double[] newValues =  distr.sample(1)[0];
        parameterList.splitParameter(categoryIndex, categoryIndex+1, newValues);
        return -distr.calcLogP(new RealParameter(newValues));
    }





    private double merge(ParameterList parameterList) throws Exception{
        int l = pointers.getDimension();

        int mergeIndex = Randomizer.nextInt(boundaryPoints.getDimension() - 1);
        int newCategoryCount = boundaryPoints.getDimension() - 1 ;
        double logq = 0;
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();
        ParameterList alphaList = alphaListInput.get();
        ParameterList invPList = invPListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList siteModelList = siteModelListInput.get();
        double temp = boundaryPoints.getValue(mergeIndex+1,1) - boundaryPoints.getValue(mergeIndex,0);
        double temp1 = boundaryPoints.getValue(mergeIndex,0);
        double temp2 = boundaryPoints.getValue(mergeIndex,1);
        double temp3 = boundaryPoints.getValue(mergeIndex+1,0);
        double temp4 = boundaryPoints.getValue(mergeIndex+1,1);
        logq += mergeParameter(parameterList,mergeIndex, parametricDistribution);

        logq += mergeParameter(freqsList,mergeIndex,freqsDistribution);

        logq += mergeParameter(ratesList,mergeIndex,ratesDistribution);

        logq += mergeParameter(alphaList,mergeIndex, alphaDistribution);
        logq += mergeParameter(invPList,mergeIndex, invPrDistribution);
        logq += mergeParameter(modelList,mergeIndex,substModelIndicatorDistr);
        logq += mergeParameter(siteModelList,mergeIndex,siteModelIndicatorDistr);

        int removedCatStart = (int)boundaryPoints.getValue(mergeIndex + 1,0);
        int removedCatEnd = (int)boundaryPoints.getValue(mergeIndex + 1,1);
        int[] sitesToMerge = new int[removedCatEnd - removedCatStart + 1];
        for(int i = 0; i < sitesToMerge.length; i++){
            sitesToMerge[i] = removedCatStart + i;
        }

        //System.out.println("sites changed: "+sitesToMerge.length+" "+((double)sitesToMerge.length/(temp)) + " " +temp1 + " "+temp2+" "+temp3+" " +temp4);
        pointers.multiPointerChanges(sitesToMerge,removedCatStart -1);

        double mergedCategoryUpper =  boundaryPoints.getValue(mergeIndex + 1,1);
        //boundaryPoints.setValue(mergeIndex, 1, mergedCategoryUpper);
        //boundaryPoints.removeParameter(mergeIndex + 1);

        boundaryPoints.mergeParameter(mergeIndex+1,mergeIndex,1,mergedCategoryUpper);

        int[] suitableCategories = new int[boundaryPoints.getDimension()];
        int suitableCategoryCount = 0;
        for(int i = 0; i< suitableCategories.length; i++){
            //System.out.println(boundaryPoints.getParameter(i));
            if(((int)(boundaryPoints.getValue(i,1) - (int)boundaryPoints.getValue(i,0))) >= 1){
                suitableCategories[suitableCategoryCount++] = i;
            }

        }
        //double potentialSplitIndexCount = boundaryPoints.getValue(mergeIndex, 1) - boundaryPoints.getValue(mergeIndex, 0);

        int suitableCategoryCountTotal = 0;
        for(int i = 0; i < suitableCategoryCount;i++){
            suitableCategoryCountTotal += (int)(boundaryPoints.getValue(suitableCategories[i],1) - (int)boundaryPoints.getValue(suitableCategories[i],0));
        }
        //System.out.println(parameterList);
        //logq +=Math.log(l -  newCategoryCount) - Math.log(newCategoryCount) + Math.log(newCategoryCount)-Math.log(pointers.getDimension()-1) ;


        //logq += Math.log(l -  newCategoryCount)-Math.log(suitableCategoryCount) - Math.log(potentialSplitIndexCount);
        //logq += - Math.log(newCategoryCount)+Math.log(l -  newCategoryCount) - Math.log(suitableCategoryCountTotal)+Math.log(newCategoryCount);
        logq += -Math.log(suitableCategoryCountTotal)+Math.log(newCategoryCount);
        return logq;
    }

    public double mergeParameter(ParameterList parameterList, int mergeIndex, ParametricDistribution distr) throws Exception{
        RealParameter oldValues2 = parameterList.getParameter(mergeIndex + 1);

        parameterList.mergeParameter(mergeIndex + 1, mergeIndex);
        return distr.calcLogP(oldValues2);

    }



    public double logit(double p){
        return Math.log(p/(1.0-p));
    }
    public double invLogit(double logit){
        return 1.0/(1.0+Math.exp(-logit));
    }
}
