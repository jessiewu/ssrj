package beast.core;

import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.matrixAlgebra1.Matrix;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
public class ParameterListDiffLogger  extends Plugin implements Loggable{
    public Input<ParameterList> parameterListInput = new Input<ParameterList>("parameterList","A list of parameters and its differences are computed.", Input.Validate.REQUIRED);

    ParameterList parameterList;

    public void initAndValidate(){
        parameterList = parameterListInput.get();
    }

    public void init(PrintStream out){
        String prefix = "diff("+parameterList.getID()+").";
        for(int i = 0; i < parameterList.getDimension(); i++){
            out.print(prefix+i+"\t");
        }

    }


    public void log(int nSample, PrintStream out){
        int dim = parameterList.getDimension();
        for(int i = 0; i < dim; i++){
            out.print((parameterList.getValue(i,1) - parameterList.getValue(i,0) + 1)+"\t");
        }
    }

    public void close(PrintStream out){

    }

}
