package scaniter;


import modi.Constants;
import modi.PTMDB;


public class ScanContext__ {

    private final double[] reporterMassOfIsobaricTag;  // marked as ing. for now. (setSpectrum())
    private final double ppmTolerance;  // marked as ing. for now. (setSpectrum())
    private final int rangeForIsotopeIncrement;  // marked as ing. for now. (setSpectrum())
    private double massToleranceForDenovo; //important marked
    private double gapAccuracy; //important marked
    private double nonModifiedDelta; //important marked
    private double precursorTolerance; //important marked
    private double precursorAccuracy; //important marked
    private int maxNoOfC13; //important marked
    private double gapTolerance; //important marked
    private double minModifiedMass;
    private double maxModifiedMass;
    private PTMDB variableModifications;
    private PTMDB fixedModifications;
    private ScanContext__(
            double precursorTolerance,
            double precursorAccuracy,
            double gapTolerance,
            double gapAccuracy,
            double nonModifiedDelta,
            int maxNoOfC13,
            double massToleranceForDenovo,
            double[] reporterMassOfIsobaricTag,
            double ppmTolerance,
            int rangeForIsotopeIncrement,
            double minModifiedMass, double maxModifiedMass, PTMDB variableModifications, PTMDB fixedModifications) {
        this.precursorTolerance = precursorTolerance;
        this.precursorAccuracy = precursorAccuracy;
        this.gapTolerance = gapTolerance;
        this.gapAccuracy = gapAccuracy;
        this.nonModifiedDelta = nonModifiedDelta;
        this.maxNoOfC13 = maxNoOfC13;
        this.massToleranceForDenovo = massToleranceForDenovo;
        this.reporterMassOfIsobaricTag = reporterMassOfIsobaricTag;
        this.ppmTolerance = ppmTolerance;
        this.rangeForIsotopeIncrement = rangeForIsotopeIncrement;
        this.minModifiedMass = minModifiedMass;
        this.maxModifiedMass = maxModifiedMass;
        this.variableModifications = variableModifications;
        this.fixedModifications = fixedModifications;
    }

    public static ScanContext__ fromConstants() {
        return new ScanContext__(
                Constants.precursorTolerance,
                Constants.precursorAccuracy,
                Constants.gapTolerance,
                Constants.gapAccuracy,
                Constants.nonModifiedDelta,
                Constants.maxNoOfC13,
                Constants.massToleranceForDenovo,
                Constants.reporterMassOfIsobaricTag,
                Constants.PPMTolerance,
                Constants.rangeForIsotopeIncrement,
                Constants.minModifiedMass, Constants.maxModifiedMass,
                Constants.variableModifications, Constants.fixedModifications
        );
    }

    public PTMDB getVariableModifications() {
        return variableModifications;
    }

    public PTMDB getFixedModifications() {
        return fixedModifications;
    }

    public double getMinModifiedMass() {
        return minModifiedMass;
    }

    public void setMinModifiedMass(double minModifiedMass) {
        this.minModifiedMass = minModifiedMass;
    }

    public double getMaxModifiedMass() {
        return maxModifiedMass;
    }

    public void setMaxModifiedMass(double maxModifiedMass) {
        this.maxModifiedMass = maxModifiedMass;
    }

    // getter 메서드들
    public double getPrecursorTolerance() {
        return precursorTolerance;
    }

    public void setPrecursorTolerance(double precursorTolerance) {
        this.precursorTolerance = precursorTolerance;
    }

    public double getPrecursorAccuracy() {
        return precursorAccuracy;
    }

    public void setPrecursorAccuracy(double precursorAccuracy) {
        this.precursorAccuracy = precursorAccuracy;
    }

    public double getGapTolerance() {
        return gapTolerance;
    }

    public void setGapTolerance(double _gapTolerance) {
        this.gapTolerance = _gapTolerance;
    }

    public double getGapAccuracy() {
        return gapAccuracy;
    }

    public void setGapAccuracy(double gapAccuracy) {
        this.gapAccuracy = gapAccuracy;
    }

    public double getNonModifiedDelta() {
        return nonModifiedDelta;
    }

    public void setNonModifiedDelta(double nonModifiedDelta) {
        this.nonModifiedDelta = nonModifiedDelta;
    }

    public int getMaxNoOfC13() {
        return maxNoOfC13;
    }

    public void setMaxNoOfC13(int maxNoOfC13) {
        this.maxNoOfC13 = maxNoOfC13;
    }

    public double getMassToleranceForDenovo() {
        return massToleranceForDenovo;
    }

    public void setMassToleranceForDenovo(double massToleranceForDenovo) {
        this.massToleranceForDenovo = massToleranceForDenovo;
    }

    public double[] getReporterMassOfIsobaricTag() {
        return reporterMassOfIsobaricTag;
    }

    public double getPpmTolerance() {
        return ppmTolerance;
    }

    public int getRangeForIsotopeIncrement() {
        return rangeForIsotopeIncrement;
    }

    // 새로운 값으로 업데이트된 ScanContext 반환하는 with 메서드
    public ScanContext__ withPrecursorTolerance(double precursorTolerance) {
        return new ScanContext__(
                precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.minModifiedMass, this.maxModifiedMass, this.variableModifications, this.fixedModifications
        );
    }

    public ScanContext__ withPrecursorAccuracy(double precursorAccuracy) {
        return new ScanContext__(
                this.precursorTolerance, precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.minModifiedMass, this.maxModifiedMass, this.variableModifications, this.fixedModifications
        );
    }

    //


    public ScanContext__ withMaxNoOfC13(int maxNoOfC13) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, maxNoOfC13, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.minModifiedMass, this.maxModifiedMass, this.variableModifications, this.fixedModifications
        );
    }

    public ScanContext__ withGapTolerance(double gapTolerance) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.minModifiedMass, this.maxModifiedMass, this.variableModifications, this.fixedModifications
        );
    }

    public ScanContext__ withNonModifiedDelta(double nonModifiedDelta) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                nonModifiedDelta, this.maxNoOfC13, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.minModifiedMass, this.maxModifiedMass, this.variableModifications, this.fixedModifications
        );
    }

}
