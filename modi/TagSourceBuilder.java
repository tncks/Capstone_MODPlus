package modi;


import scaniter.ScanContext__;

public class TagSourceBuilder {

    public	TagSourceBuilder() {}
    public TagPool buildTagPool( Spectrum sourceSpec, ScanContext__ context ) {
        if( sourceSpec == null ) return null;

        for(int i=0; i<sourceSpec.size(); i++){
            if( Math.abs(sourceSpec.get(i).getMass() - sourceSpec.getPrecursor() ) < 2 ) {
                sourceSpec.remove(i);
                i--;
                continue;
            }
            sourceSpec.get(i).setIndex(i);
        }//temporal

        sourceSpec.normalizeIntensityLocally();

        int extra = ( sourceSpec.getCharge() > 2 && Constants.INSTRUMENT_TYPE != Constants.msms_type.QTOF )? 2 : 0;
        sourceSpec.peakSelection(Constants.selectionWindowSize, Constants.minNumOfPeaksInWindow+extra, context.getMassToleranceForDenovo() );

        return sourceSpec.generateTags(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain, context.getMassToleranceForDenovo());
    }


}

