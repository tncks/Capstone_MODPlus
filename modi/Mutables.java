package modi;

import msutil.ProtCutter;

public class Mutables {
    public static int			targetDecoy=0;

    public static String 		firstSearchProgram = "";

    public static double		NTERM_FIX_MOD = 0;
    public static double		CTERM_FIX_MOD = 0;

    public static ProtCutter protease = ProtCutter.getCutter("Trypsin");
    public static int			numberOfEnzymaticTermini = 2;
    public static int			missCleavages = 2;

    public static int			minNoOfC13 = 0;
    public static int			maxNoOfC13 = 0;
    public static int			rangeForIsotopeIncrement = 0;

    public static double		alkylatedToCys = 0;
    public static String		alkylationMethod;

    public static double		precursorAccuracy = 0.5;
    public static double		precursorTolerance = 0.5;
    public static double		PPMTolerance = 0;
    public static double		fragmentTolerance = 0.6;
    public static double		gapTolerance = 0.6;
    public static double		gapAccuracy = 1.6;


    public static PTMDB 		variableModifications;
    public static PTMDB 		fixedModifications;
    public static double		minModifiedMass = -precursorTolerance;
    public static double		maxModifiedMass = precursorTolerance;
    public static boolean		canBeModifiedOnFixedAA = false;



    //for De novo sequencing
    public static double		massToleranceForDenovo = 0.3;


    public static int			minNumOfPeaksInWindow = 4;


    public static boolean		Leu_indistinguishable_Ile = true;
    public static boolean		Lys_indistinguishable_Qln = true;



    public static int			maxPTMPerGap		= 2;
    public static int 			maxPTMPerPeptide	= 4;

    public static String		PTM_FILE_NAME = "PTMDB.xml";

    public static double		nonModifiedDelta = massToleranceForDenovo;

    public static String		isobaricTag = "";
    public static double[]		reporterMassOfIsobaricTag = null;

    public static String		enrichedModification = "";


    public static int getMaxPTMOccurrence( int seqLength ){
        if( seqLength > 10 ) return 1;
        return Mutables.maxPTMPerGap;
    }

    public static boolean	fEqual(double v1, double v2){
        return Math.abs(v1 - v2) <= Mutables.fragmentTolerance;
    }


    @SuppressWarnings("BooleanMethodIsAlwaysInverted")
    public static boolean		isInModifiedRange(double v ){
        if( Mutables.minModifiedMass-Mutables.gapTolerance < v && v < Mutables.maxModifiedMass+Mutables.gapTolerance ) return true;
        else return Math.abs(v) <= Mutables.gapTolerance;
    }

    public static boolean isWithinTolerance(double calc, double obsv, double tol){

        if( Mutables.minNoOfC13 ==0 && Mutables.maxNoOfC13 == 0 ) {
            return !(Math.abs(calc - obsv) > tol);
        }
        else {
            double tempError = obsv - calc;
            int isoerr = Constants.round( tempError / Constants.IsotopeSpace );
            if( isoerr < Mutables.minNoOfC13 || Mutables.maxNoOfC13 < isoerr ) return false;
            return !(Math.abs(tempError - isoerr * Constants.IsotopeSpace) > Mutables.precursorAccuracy);
        }
    }

    public static boolean isWithinAccuracy(double err){
        if( Mutables.gapAccuracy > 0.5 ) return true;
        int isoerr = Constants.round( err / Constants.IsotopeSpace );
        return !(Math.abs(err - isoerr * Constants.IsotopeSpace) > Mutables.gapAccuracy);
    }


}
