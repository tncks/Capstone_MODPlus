package modi;

import msutil.PGraph;
import processedDB.TagTrie;
import scaniter.ScanContext__;

import java.util.ArrayList;

public class DynamicDBHolder {

    public DynamicDBHolder() {}

    public ArrayList<AnsPeptide> dynamicMODeye(TagTrie dynamicDB, PGraph graph, TagPool tPool, ScanContext__ context) {
        SpectrumAnalyzer szer = new SpectrumAnalyzer();
        MatchedTagPool matchedList = szer.extendedBuildMatchedTagPool(tPool, graph.getCorrectedMW(),
                dynamicDB, Constants.protease, Constants.numberOfEnzymaticTermini, context);

        TagChainPool tcPool = new TagChainPool();
        tcPool.putAll(szer.buildTagChain(matchedList, context));
        tcPool.discardPoorTagChain();

        boolean specAnnotated = false;
        if (tcPool.size() != 0) {
            specAnnotated = szer.interpretTagChain(context.getVariableModifications(), tcPool, graph);
        }

        ArrayList<AnsPeptide> cands = new ArrayList<>();
        if (tcPool.size() != 0 && specAnnotated) {
            cands = tcPool.getAnswerPeptides(graph, context);
        }
        return cands;
    }
}
