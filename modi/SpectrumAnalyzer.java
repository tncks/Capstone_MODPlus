package modi;

import msutil.PGraph;
import msutil.ProtCutter;
import processedDB.RetrivedPeptideMap;
import processedDB.TagTrie;
import scaniter.ScanContext__;

import java.util.LinkedList;


public class SpectrumAnalyzer {

    public SpectrumAnalyzer() {
    }


    public TagChainPool buildTagChain(MatchedTagPool matchedTags, ScanContext__ context) {
        TagChainPool tagChainPool = new TagChainPool();
        tagChainPool.buildTagChainPool(matchedTags, context);
        return tagChainPool;
    }

    public boolean interpretTagChain(PTMDB ptmDB, TagChainPool tcPool, PGraph graph) {
        Spectrum sourceSpectrum = null;
        boolean specAnnotated = false;

        for (LinkedList<TagChain> tagChainList : tcPool.values()) {
            for (int k = 0; k < tagChainList.size(); k++) {
                TagChain tc = tagChainList.get(k);

                boolean allGapAnnotated = true;
                if (sourceSpectrum == null) {
                    sourceSpectrum = tc.getSourceSpectrum();
                }
                Peptide pep = tc.getMatchedPeptide();
                for (SpecInterpretation si : tc) {
                    if (!(si instanceof Gap gap)) continue;
                    PTMSearchResult interpretation = ptmDB.searchPTM(pep.subSequence(gap.getStart(), gap.getEnd() + 1),
                            gap.getOffset(), gap.getPosition());

                    if (!interpretation.isInterpreted()) {
                        gap.setInterpreted(false);
                        allGapAnnotated = false;
                        tc.setAllGapAnnotated(false);
                        break;
                    } else gap.setInterpreted(true);

                    gap.setInterpretation(interpretation, graph);
                }

                if (allGapAnnotated) {
                    tc.setAllGapAnnotated(true);
                    specAnnotated = true;
                } else {
                    tagChainList.remove(k);
                    k--;
                }
            }
        }
        return specAnnotated;
    }

    public MatchedTagPool extendedBuildMatchedTagPool(TagPool primitiveTags, double motherMass,
                                               TagTrie ixPDB, ProtCutter enzyme, int NTT, ScanContext__ context) {
        if (primitiveTags == null || ixPDB == null)
            return null;

        double minDelta = (context.getMinModifiedMass() < 0)? context.getMinModifiedMass() - context.getGapTolerance() : -context.getGapTolerance();
        double maxDelta = (context.getMaxModifiedMass() > 0)? context.getMaxModifiedMass() + context.getGapTolerance() : +context.getGapTolerance();

        TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

        int realTag = 0;
        double orbMass = motherMass - Constants.H2O;
        RetrivedPeptideMap searchResults = new RetrivedPeptideMap();
        for (Tag tag : longTags) {

            RetrivedPeptideMap bRes = ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, tag.getBIonNtermOffset() - Constants.NTERM_FIX_MOD, tag,
                    tag.getBIonCtermOffset() - Constants.CTERM_FIX_MOD, IonDirection.B_DIRECTION, minDelta, maxDelta, Constants.gapTolerance);
            searchResults.combine(bRes);

            Tag reverseTag = tag.reverseTag();
            RetrivedPeptideMap yRes = ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, reverseTag.getYIonNtermOffset() - Constants.NTERM_FIX_MOD, reverseTag,
                    reverseTag.getYIonCtermOffset() - Constants.CTERM_FIX_MOD, IonDirection.Y_DIRECTION, minDelta, maxDelta, Constants.gapTolerance);
            searchResults.combine(yRes);
            realTag++;

            if (realTag > Constants.MAX_TAG_SIZE * 2) break;
        }
        return searchResults.convertToMatchedTagPool(primitiveTags.extract(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain), context);
    }

}