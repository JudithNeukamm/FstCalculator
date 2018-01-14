package fst;

import methods.UsefulFunctions;

/**
 * Created by neukamm on 6/22/17.
 */
public class Linearization {

    private final UsefulFunctions functions;

    public Linearization() {
        functions = new UsefulFunctions();

    }


    /**
     * This method calculates Slatkin's linearized Fsts (Slatkin 1995) based on the previously
     * calculated Fst value:
     *
     *              D = Fst / (1-Fst)
     *
     * @param fsts
     * @return linearized Fsts
     */
    public double[][] linearizeWithSlatkin(double[][] fsts){

        double[][] linearizedFsts = new double[fsts.length][fsts[0].length];
        for(int i = 0; i < linearizedFsts.length; i++){
            for(int j = 0; j < linearizedFsts[i].length; j++){
                linearizedFsts[i][j] = fsts[i][j] / (1 - fsts[i][j]);
            }
        }

        return functions.transposeMatrix(linearizedFsts);
    }


    /**
     * This method calculates Reynolds' distance (Reynolds et al. 1983) based on the previously
     * calculated Fst value:
     *
     *              D = -ln(1-Fst)
     *
     * @param fsts
     * @return linearized Fsts
     */
    public double[][] linearizeWithReynolds(double[][] fsts){

        double[][] linearizedFsts = new double[fsts.length][fsts[0].length];
        for(int i = 0; i < linearizedFsts.length; i++){
            for(int j = 0; j < linearizedFsts[i].length; j++){
                linearizedFsts[i][j] = Math.log(1 - fsts[i][j]) * -1;
            }
        }

        return functions.transposeMatrix(linearizedFsts);
    }


}
