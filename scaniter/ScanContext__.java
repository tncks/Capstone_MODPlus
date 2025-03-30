package scaniter;


import modi.Constants;


public class ScanContext__ {
    private final double gapAccuracy;
    private final double fragmentTolerance;
    private final double massToleranceForDenovo;
    private final double[] reporterMassOfIsobaricTag;
    private final double proton;
    private final double isotopeSpace;
    private final double ppmTolerance;
    private final int rangeForIsotopeIncrement;
    private final int maxPTMPerPeptide;
    private double nonModifiedDelta;
    private double precursorTolerance;
    private double precursorAccuracy;
    private int    maxNoOfC13;
    private double gapTolerance;
    private double minModifiedMass;
    private double maxModifiedMass;
    private ScanContext__(
            double precursorTolerance,
            double precursorAccuracy,
            double gapTolerance,
            double gapAccuracy,
            double nonModifiedDelta,
            int maxNoOfC13,
            double fragmentTolerance,
            double massToleranceForDenovo,
            double[] reporterMassOfIsobaricTag,
            double proton,
            double isotopeSpace,
            double ppmTolerance,
            int rangeForIsotopeIncrement,
            int maxPTMPerPeptide, double minModifiedMass, double maxModifiedMass) {
        this.precursorTolerance = precursorTolerance;
        this.precursorAccuracy = precursorAccuracy;
        this.gapTolerance = gapTolerance;
        this.gapAccuracy = gapAccuracy;
        this.nonModifiedDelta = nonModifiedDelta;
        this.maxNoOfC13 = maxNoOfC13;
        this.fragmentTolerance = fragmentTolerance;
        this.massToleranceForDenovo = massToleranceForDenovo;
        this.reporterMassOfIsobaricTag = reporterMassOfIsobaricTag;
        this.proton = proton;
        this.isotopeSpace = isotopeSpace;
        this.ppmTolerance = ppmTolerance;
        this.rangeForIsotopeIncrement = rangeForIsotopeIncrement;
        this.maxPTMPerPeptide = maxPTMPerPeptide;
        this.minModifiedMass = minModifiedMass;
        this.maxModifiedMass = maxModifiedMass;
    }

    public static ScanContext__ fromConstants() {
        return new ScanContext__(
                Constants.precursorTolerance,
                Constants.precursorAccuracy,
                Constants.gapTolerance,
                Constants.gapAccuracy,
                Constants.nonModifiedDelta,
                Constants.maxNoOfC13,
                Constants.fragmentTolerance,
                Constants.massToleranceForDenovo,
                Constants.reporterMassOfIsobaricTag,
                Constants.Proton,
                Constants.IsotopeSpace,
                Constants.PPMTolerance,
                Constants.rangeForIsotopeIncrement,
                Constants.maxPTMPerPeptide, Constants.minModifiedMass, Constants.maxModifiedMass
        );
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

    public double getFragmentTolerance() {
        return fragmentTolerance;
    }

    public double getMassToleranceForDenovo() {
        return massToleranceForDenovo;
    }

    public double[] getReporterMassOfIsobaricTag() {
        return reporterMassOfIsobaricTag;
    }

    public double getProton() {
        return proton;
    }

    public double getIsotopeSpace() {
        return isotopeSpace;
    }

    public double getPpmTolerance() {
        return ppmTolerance;
    }

    public int getRangeForIsotopeIncrement() {
        return rangeForIsotopeIncrement;
    }

    public int getMaxPTMPerPeptide() {
        return maxPTMPerPeptide;
    }

    // 새로운 값으로 업데이트된 ScanContext 반환하는 with 메서드
    public ScanContext__ withPrecursorTolerance(double precursorTolerance) {
        return new ScanContext__(
                precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.maxPTMPerPeptide, this.minModifiedMass, this.maxModifiedMass
        );
    }

    public ScanContext__ withPrecursorAccuracy(double precursorAccuracy) {
        return new ScanContext__(
                this.precursorTolerance, precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.maxPTMPerPeptide, this.minModifiedMass, this.maxModifiedMass
        );
    }

    //


    public ScanContext__ withMaxNoOfC13(int maxNoOfC13) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.maxPTMPerPeptide, this.minModifiedMass, this.maxModifiedMass
        );
    }

    public ScanContext__ withGapTolerance(double gapTolerance) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, gapTolerance, this.gapAccuracy,
                this.nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.maxPTMPerPeptide, this.minModifiedMass, this.maxModifiedMass
        );
    }

    public ScanContext__ withNonModifiedDelta(double nonModifiedDelta) {
        return new ScanContext__(
                this.precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
                nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
                this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
                this.rangeForIsotopeIncrement, this.maxPTMPerPeptide, this.minModifiedMass, this.maxModifiedMass
        );
    }

}
