import IO.reader.DistanceTypeParser;
import IO.writer.Writer;
import fst.FstAMOVA;
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


    public FstCalculator(HashMap<String, List<String>> data,
                         String missing_data,
                         double level_missing_data,
                         String distance_method,
                         int number_of_permutations,
                         double gamma

                         ) {

        this.data = data;
        this.missing_data = missing_data;
        this.distance_method = distance_method;
        this.level_missing_data = level_missing_data;
        this.number_of_permutations = number_of_permutations;
        this.gamma = gamma;
    }

    public double[][] runCaclulations() throws IOException {
        return runAmova();
    }

    private double[][] runAmova() throws IOException {

        Linearization linearization = new Linearization();

        DistanceTypeParser distanceTypeParser = new DistanceTypeParser();
        Filter filter = new Filter();
        List<Integer> usableLoci = filter.getUsableLoci(data, missing_data, level_missing_data);

        FstAMOVA standardAMOVA = new FstAMOVA(
                usableLoci,
                number_of_permutations,
                distanceTypeParser.parse(distance_method),
                gamma);

        standardAMOVA.setData(data);

        double[][] fsts_amova = standardAMOVA.calculateFst();
        double[][] pvalues = standardAMOVA.calculatePermutatedFst();


        Writer writer = new Writer(number_of_permutations);

        writer.writeResultsFstToString(
                fsts_amova,
                pvalues,
                standardAMOVA.getGroupnames(),
                usableLoci,
                level_missing_data,
                significance
        );
        writer.addDistanceMatrixToResult(
                standardAMOVA.getDistanceCalculator().getDistancematrix_d()
        );

        writer.addLinerarizedFstMatrix(
                linearization.linearizeWithSlatkin(fsts_amova),
                "# Linearized Fst values (Slatkin)."
        );
        writer.addLinerarizedFstMatrix(
                linearization.linearizeWithReynolds(fsts_amova),
                "# Linearized Fst values (Reynold)."
        );

        writer.writeResultsToFile("resultsFstStatisticsAMOVA.tsv");

        return fsts_amova;

    }



}
