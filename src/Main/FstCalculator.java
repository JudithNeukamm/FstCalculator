package Main;

import IO.reader.DistanceTypeParser;
import IO.writer.Writer;
import fst.FstPhiArlequin;
import fst.Linearization;
import methods.Filter;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;

public class FstCalculator {

    private static double gamma = 0.5;
    private static double level_missing_data = 0.05;
    private static double significance = 0.05;
    private static int number_of_permutations = 0;
    private static String missing_data = "N";
    private final HashMap<String, List<String>> data;
    private String distance_method = "Pairwise Difference";
    private Writer writer;
    private String[] groupnames;
    private double[][] fsts_amova;
    private double[][] pvalues;
    private FstPhiArlequin fstPhiArlequin;
    private List<Integer> usableLoci;
    private Linearization linearization;


    public FstCalculator(HashMap<String, List<String>> data,
                         String missing_data,
                         double level_missing_data,
                         String distance_method,
                         int number_of_permutations,
                         double gamma,
                         double significance

    ) {

        this.data = data;
        this.missing_data = missing_data;
        this.distance_method = distance_method;
        this.level_missing_data = level_missing_data;
        this.number_of_permutations = number_of_permutations;
        this.gamma = gamma;
        this.significance = significance;
    }

    public double[][] runCaclulations() throws IOException {
        linearization = new Linearization();

        DistanceTypeParser distanceTypeParser = new DistanceTypeParser();
        Filter filter = new Filter();
        usableLoci = filter.getUsableLoci(data, missing_data, level_missing_data);

        fstPhiArlequin = new FstPhiArlequin(
                usableLoci,
                number_of_permutations,
                distanceTypeParser.parse(distance_method),
                gamma);

        fstPhiArlequin.setData(data);

        fsts_amova = fstPhiArlequin.calculateFst();
        pvalues = null;
        if(number_of_permutations > 0){
            pvalues = fstPhiArlequin.calculatePermutatedFst();
        }

        groupnames = fstPhiArlequin.getGroupnames();

        return fsts_amova;

    }


    public void writeResultToFile(String filepath) throws IOException {


        writer = new Writer(number_of_permutations);

        writer.writeResultsFstToString(
                fsts_amova,
                pvalues,
                groupnames,
                usableLoci,
                level_missing_data,
                significance
        );
        writer.addDistanceMatrixToResult(
                fstPhiArlequin.getDistanceCalculator().getDistancematrix_d()
        );

        writer.addLinerarizedFstMatrix(
                linearization.linearizeWithSlatkin(fsts_amova),
                "# Linearized Fst values (Slatkin)."
        );
        writer.addLinerarizedFstMatrix(
                linearization.linearizeWithReynolds(fsts_amova),
                "# Linearized Fst values (Reynold)."
        );

        writer.writeResultsToFile(filepath + "/resultsFstStatisticsAMOVA.tsv");


    }

    public String getResultString(){
        return writer.getResult_as_string();
    }

    public String[] getGroupnames() {
        return groupnames;
    }
}
