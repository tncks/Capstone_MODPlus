package modi;

import msutil.ProtCutter;

public class Mutables {
    public static ProtCutter protease = ProtCutter.getCutter("Trypsin");
    public static double fragmentTolerance = 0.6; // never change after set param invocation
    public static PTMDB variableModifications;
    public static PTMDB fixedModifications;
    public static double[] reporterMassOfIsobaricTag = null;
    public static double massToleranceForDenovo = 0.3;

    public static int    maxNoOfC13 = 0;
    public static double precursorAccuracy = 0.5;
    public static double precursorTolerance = 0.5;
    public static double gapTolerance = 0.6;
    public static double gapAccuracy = 1.6;
    public static double nonModifiedDelta = massToleranceForDenovo;






    public static boolean fEqual(double v1, double v2) {
        return Math.abs(v1 - v2) <= Mutables.fragmentTolerance;
    }


    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public static boolean isInModifiedRange(double v) {
        if (Constants.minModifiedMass - Mutables.gapTolerance < v && v < Constants.maxModifiedMass + Mutables.gapTolerance)
            return true;
        else return Math.abs(v) <= Mutables.gapTolerance;
    }

    public static boolean isWithinTolerance(double calc, double obsv, double tol) {

        if (Constants.minNoOfC13 == 0 && Mutables.maxNoOfC13 == 0) {
            return !(Math.abs(calc - obsv) > tol);
        } else {
            double tempError = obsv - calc;
            int isoerr = Constants.round(tempError / Constants.IsotopeSpace);
            if (isoerr < Constants.minNoOfC13 || Mutables.maxNoOfC13 < isoerr) return false;
            return !(Math.abs(tempError - isoerr * Constants.IsotopeSpace) > Mutables.precursorAccuracy);
        }
    }

    public static boolean isWithinAccuracy(double err) {
        if (Mutables.gapAccuracy > 0.5) return true;
        int isoerr = Constants.round(err / Constants.IsotopeSpace);
        return !(Math.abs(err - isoerr * Constants.IsotopeSpace) > Mutables.gapAccuracy);
    }

    public static void printGlobalMutables() {
        System.out.println("*****************************");
        System.out.println("*****************************");
        System.out.println("maxNoOfC13: " + maxNoOfC13); //analyze
        System.out.println("precursorAccuracy: " + precursorAccuracy);
        System.out.println("precursorTolerance: " + precursorTolerance);
        System.out.println("gapTolerance: " + gapTolerance);
        System.out.println("gapAccuracy: " + gapAccuracy);
        System.out.println("nonModifiedDelta: " + nonModifiedDelta);
        System.out.println("*****************************");
        System.out.println("*****************************");
    }


}
