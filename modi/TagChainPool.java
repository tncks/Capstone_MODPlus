package modi;

import java.util.*;

import msutil.*;
import scaniter.ScanContext__;

public class TagChainPool extends TreeMap<Peptide, LinkedList<TagChain>> {

	public LinkedList<TagChain> buildTagChainList(Map.Entry<Peptide, LinkedList<MatchedTag>> tagPool){
		initTagMerge(tagPool);
		return enumTagChain(tagPool);
	}

	public MatchedTag getMatchedB2Tag(Tag tag, Peptide pept)
	{
		MatchedTag b2m = new MatchedTag(tag);
		b2m.matchedPeptide = pept;
		b2m.start = 0;
		b2m.end = 1;
		b2m.tagSequence = pept.subSequence(0, 2);
		b2m.direction = IonDirection.B_DIRECTION;
		return b2m;
	}

	public void initTagMerge( Map.Entry<Peptide, LinkedList<MatchedTag>> entry )
	{
		LinkedList<MatchedTag> tagList = entry.getValue();
		tagList.sort(new SpecInterpretationComparator());
		Spectrum ccspec = tagList.getFirst().getSourceSpectrum();

		if( tagList.getFirst().start != 0 ){
			Peptide pept= entry.getKey();
			Tag gto= ccspec.getB2Tag(pept.subSequence(0, 2));
			if( gto != null ){
				MatchedTag b2mTag= getMatchedB2Tag(gto, pept);
				tagList.addFirst( b2mTag );
			}
		}//*/

		int initSize  = tagList.size();
		for(int i=0; i<initSize-1; i++)
		{
			MatchedTag seed = tagList.get(i);
			for(int j=i+1; j<initSize; j++)
			{
				MatchedTag tag = tagList.get(j);
				if( seed.getRelativePosition(tag) == RelativePosition.SEPERATED )
					break;

				if( seed.getDirection() == tag.getDirection() &&
						Constants.fEqual( seed.getOffset(), tag.getOffset() ) ){
					MatchedTag temp= new MatchedTag(seed);
					temp.extend(tag);
					if( temp.size() != temp.sequence().size()+1 )
						continue;
					seed.extend(tag);	// extend seed from extendable TagList
					tagList.remove(j);
					j--;
					initSize--;
					continue;
				}

				if( seed.isComplementarySame(tag) ){
					tagList.remove(j);		// remove seed from matchedTagList
					j--;
					initSize--;
				}
			}
		}

		int countOfLongTag= 0;
		MatchedTag sLongTag= null;
		for(int k=0; k<tagList.size(); k++){
			tagList.get(k).setScore();
			if( tagList.get(k).sequence().size() > 2 ){
				countOfLongTag++;
				sLongTag= tagList.get(k);
			}
		}

//		topept++;
		if( tagList.size() > Constants.maxTagPerPept ){
			//	maxHitperPept++;
			tagList.sort(Collections.reverseOrder(new TagComparator()));
			for(int k=Constants.maxTagPerPept; k<tagList.size(); k++){
				tagList.remove(k--);
			}
			tagList.sort(new SpecInterpretationComparator());
		}//*/

		if( countOfLongTag == 1 && tagList.size() < 4 ){
			int rev= 0;
			MatchedTag newTag= new MatchedTag(sLongTag);
			if( sLongTag.getDirection() == IonDirection.B_DIRECTION ){
				//	if( !Constants.fEqual(sLongTag.getNTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(s), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getStart(), s);
				}
				//	if( !Constants.fEqual(sLongTag.getCTermOffset(), 0) )

				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(sLongTag.size()-s-1), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getEnd(), s);
				}
			}
			else{
				//	if( !Constants.fEqual(sLongTag.getCTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(s), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getEnd(), s);
				}
				//	if( !Constants.fEqual(sLongTag.getNTermOffset(), 0) )
				{
					rev++;
					int s;
					for(s=1; sLongTag.sequence().size()-rev > 2 ;s++){
						if( ccspec.isConfidentPeak(sLongTag.get(sLongTag.size()-s-1), 1) ) break;
						rev++;
					}
					newTag.trim(newTag.getStart(), s);
				}
			}
			if( rev > 0 ){
				if( newTag.size() > 2 ) tagList.add(newTag);
			}
		}
	}

	public LinkedList<TagChain> enumTagChain(Map.Entry<Peptide, LinkedList<MatchedTag>> entry)
	{
		LinkedList<MatchedTag> tagList = entry.getValue();
		tagList.sort(new SpecInterpretationComparator());	// sort matched tag list by start position

		LinkedList<TagChain> tagChainList = new LinkedList<>();

		//	System.out.println("TAGCHAIN.JAVA : " + entry.getKey()+ " " + tagList.size() );
		for( MatchedTag tag : tagList ) {
			TagChain t = new TagChain(entry.getKey(), tag.getSourceSpectrum());
			t.add(tag);
			tagChainList.add(t);
			//		System.out.println(tag);
		}//generate base tag_chain
		//	System.out.println(entry.getKey() + " " + tagList.size());

		// tag chain enumeration
		ArrayList<TagChain> addedTagChain = new ArrayList<>();
		int start = 0, end = tagChainList.size(), addedCount;

		while(true)
		{
			ListIterator<TagChain> listIt = tagChainList.listIterator(start);
			while(listIt.hasNext())
			{
				TagChain curTC = listIt.next();
				if( !(curTC.last() instanceof MatchedTag seed) )
					continue;

				for( MatchedTag tag : tagList ) {
					if( seed.compareTo(tag) > 0 || curTC.useSamePeak(tag) ) continue;

					combineTagChains(addedTagChain, curTC, tag);
				}
			}
			addedCount = addedTagChain.size();
			tagChainList.addAll(addedTagChain);
			addedTagChain.clear();
			start = end;
			end += addedCount;
			if( addedCount == 0 ) break;
		}

		// make gap for all tag chain list
		ListIterator<TagChain> listIt = tagChainList.listIterator();
		double topSCore = 0;
		while( listIt.hasNext() ){
			TagChain curTC = listIt.next();
			if( curTC.makeGap() ){
				curTC.setTagChainScore();
				if( curTC.score < 0 ) listIt.remove();
				if( curTC.score > topSCore )
					topSCore = curTC.score;
			}
			else{ listIt.remove(); }
		}

		listIt = tagChainList.listIterator();
		while( listIt.hasNext() ){
			TagChain curTC = listIt.next();
			if( curTC.score < topSCore * Constants.tagChainPruningRate ){
				listIt.remove();
			}
		}

		if( tagChainList.size() > Constants.maxTagChainPerPept ){
			Collections.sort( tagChainList );
			for(int i=Constants.maxTagChainPerPept; i<tagChainList.size(); i++)
				tagChainList.remove(i--);
		}

		return tagChainList;
	}


	public void  combineTagChains (ArrayList<TagChain> newTCList, TagChain baseTC, MatchedTag tag) {
		if( !(baseTC.last() instanceof MatchedTag seed) )
			return;

		TagChain newTC = (TagChain)baseTC.clone();

		RelativePosition ir= seed.getRelativePosition(tag);
		if( ir == RelativePosition.SEPERATED ){
			newTC.add(tag);
			newTCList.add(newTC);
		}

		else if( ir == RelativePosition.ADJACENT ){

			if(	seed.size() < 3 )
				return;

			if(	seed.sequence().size() < Constants.minTagLengthPeptideShouldContain &&
					tag.sequence().size()< Constants.minTagLengthPeptideShouldContain )
				return;

			newTC.remove(newTC.last());

			if( LinkableTags(seed, tag) ){
				if( seed.size() < tag.size() )
					newTCList.add( extendTagChain(newTC, seed, tag.getTrimedTag( tag.getStart(), 1 )) );
				else
					newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), 1 ), tag) );
			}
			else{
				if( tag.size() > 2 ) newTCList.add( extendTagChain(newTC, seed, tag.getTrimedTag( tag.getStart(), 1 )) );
				newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), 1 ), tag) );
			}
		}

		else if( ir == RelativePosition.OVERLAP ){

			if( seed.sequence().size()< Constants.minTagLengthPeptideShouldContain &&
					tag.sequence().size()< Constants.minTagLengthPeptideShouldContain )
				return;

			newTC.remove(newTC.last());

			int overLabPart= seed.getEnd() - tag.getStart() + 1;

			newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), overLabPart ),
					tag.getTrimedTag( tag.getStart(), overLabPart )) );
			if( LinkableTags(seed, tag) ){ return; }

			if( tag.sequence().size() > overLabPart+1 )
				newTCList.add( extendTagChain(newTC, seed,
						tag.getTrimedTag( tag.getStart(), overLabPart+1 )) );
			if( seed.sequence().size() > overLabPart+1 )
				newTCList.add( extendTagChain(newTC, seed.getTrimedTag( seed.getEnd(), overLabPart+1 ),
						tag) );
		}

		else if( ir == RelativePosition.INCLUDING ){

			if(	seed.sequence().size() < Constants.minTagLengthPeptideShouldContain &&
					tag.sequence().size() < Constants.minTagLengthPeptideShouldContain )
				return;
			if( Constants.fEqual( seed.getOffset(), tag.getOffset()) ) return;
			if( !seed.isLikelyChild(tag) ) return;

			int rev= 0;
			MatchedTag newTag;
			if( seed.score > tag.score )
				newTag= new MatchedTag(seed);
			else
				newTag= new MatchedTag(tag);

			if( seed.getDirection() == IonDirection.B_DIRECTION ){
				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.getFirst(), 1) ) break;
					newTag.trim(newTag.getStart(), 1);
					rev++;
				}

				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.getLast(), 1) ) break;
					newTag.trim(newTag.getEnd(), 1);
					rev++;
				}
			}
			else {
				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.getFirst(), 1) ) break;
					newTag.trim(newTag.getEnd(), 1);
					rev++;
				}

				while( newTag.sequence().size() > 2 ) {
					if( newTag.sourceSpectrum.isConfidentPeak(newTag.getLast(), 1) ) break;
					newTag.trim(newTag.getStart(), 1);
					rev++;
				}
			}

			if( rev > 0 ){
				newTC.remove(newTC.last());
				newTC.add(newTag);
				newTCList.add(newTC);
			}
		}//INCLUDING*/
	}

	public TagChain extendTagChain (TagChain baseTC, MatchedTag one, MatchedTag two){
		TagChain newTC = (TagChain)baseTC.clone();
		newTC.add(one);
		newTC.add(two);
		return newTC;
	}

	public boolean LinkableTags(MatchedTag one, MatchedTag two){
		//ADJACENT tags citation
		if( one.getDirection() == two.getDirection() ) return false;
		return Constants.fEqual(one.getNTermOffset(), two.getNTermOffset());
	}


	public void buildTagChainPool(MatchedTagPool matchedTags)
	{
		if( matchedTags == null || matchedTags.isEmpty()){
			return;
		}
		
		LinkedList<TagChain> tcList;
		
		Iterator<Map.Entry<Peptide, LinkedList<MatchedTag>>> it = matchedTags.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<MatchedTag>> entry = it.next();
			tcList = buildTagChainList(entry);
			if(tcList.isEmpty()){ continue; }
			Peptide matchedPeptide = new Peptide(entry.getKey());
			put(matchedPeptide, tcList);
		}
	}

	public void discardPoorTagChain()	// should be modified
	{
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();

		// calculating best score
		double bestScore = 0;
		// TagChain bet = null;
		while( it.hasNext() ) {
			for(TagChain tc : it.next().getValue())
			{				
				double curScore = tc.getScore();
				if( bestScore < curScore ){
					bestScore = curScore;
					// bet = tc;
				}
			}
		}
	
		// discarding tag chain		
		it = this.entrySet().iterator();
		ArrayList<Peptide> deletePepList = new ArrayList<>();
		
		while(it.hasNext())
		{
			LinkedList<TagChain> tcList = it.next().getValue();
			Peptide pep = tcList.getFirst().getMatchedPeptide();
			ListIterator<TagChain> listIt = tcList.listIterator();
			while(listIt.hasNext()) {
				TagChain tc = listIt.next();
				if( tc.getScore() <= bestScore * Constants.tagChainPruningRate ){
					listIt.remove();
				}	
			}
			
			if(tcList.isEmpty())
				deletePepList.add(pep);
		}
		
		for( Peptide pep : deletePepList )
			this.remove(pep);	
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		output.append("Matched Tag Pool");
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> it = this.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<Peptide, LinkedList<TagChain>> entry = it.next();
			output.append("\n").append(entry.getKey()).append("\n");
			for(TagChain tc : entry.getValue())
			{
				output.append(tc).append("\n");
			}
		}
		return output.toString();
	}

	public int round(double a){
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}

	public boolean isWithinTolerance(double calc, double obsv, double tol, ScanContext__ context){

		if( Constants.minNoOfC13 == 0 && context.getMaxNoOfC13() == 0 ) {

			return !(Math.abs(calc - obsv) > tol);
		}
		else {
			double tempError = obsv - calc;
			int isoerr = round( tempError / context.getIsotopeSpace() );

			if( isoerr < Constants.minNoOfC13 || context.getMaxNoOfC13() < isoerr ) return false;
			return !(Math.abs(tempError - isoerr * context.getIsotopeSpace()) > context.getPrecursorAccuracy());
		}
	}


	public int getModEyeRankScore( String peptide, double[] ptms, PGraph graph, ScanContext__ context ){	//for modeye preliminary ranking
		IonGraph iGraph;
		if( Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF ) iGraph = new TOFGraph(peptide, ptms, graph);
		else iGraph= new TRAPGraph(peptide, ptms, graph);

		if( !isWithinTolerance(iGraph.getCalculatedMW(), graph.getObservedMW(), context.getPrecursorTolerance(), context) )
			return -1;

		iGraph.setScore(graph);
		return iGraph.getRankScore();
	}

	public ArrayList<AnsPeptide> getAnswerPeptides(PGraph graph, ScanContext__ context){
		AnsHeap answerPepts = new AnsHeap();
		
		Iterator<Map.Entry<Peptide, LinkedList<TagChain>>> entries = this.entrySet().iterator();
		while(entries.hasNext()){
			Map.Entry<Peptide, LinkedList<TagChain>> entry = entries.next();
			String pept = entry.getKey().toString();
			LinkedList<TagChain> alignedTagChainList = entry.getValue();
			if(alignedTagChainList.isEmpty()) continue;
			
			HashSet<PTMCombination> ptmComb = new HashSet<>();
			for( TagChain tc : alignedTagChainList ){
				if( tc.allGapAnnotated ) ptmComb.addAll( tc.getPTMCombination() );
			}
			
			Iterator<PTMCombination> iter = ptmComb.iterator();
			while( iter.hasNext() ){
				PTMCombination p = iter.next();
	
				int s = getModEyeRankScore(pept, p.ptms, graph, context);
				if( s < 0 ) continue;
				AnsPeptide candidate = new AnsPeptide(entry.getKey(), p.ptmComb, p.ptms, p.ptmList, s);
				answerPepts.add( candidate ); 
			}
		}
		
		return answerPepts.getFinalList(graph, context);
	}
	
}
	
class TagChainListComparator implements Comparator<Map.Entry<Peptide, LinkedList<TagChain>>>
{
	public int compare(Map.Entry<Peptide, LinkedList<TagChain>> tc1, Map.Entry<Peptide, LinkedList<TagChain>> tc2)
	{
		LinkedList<TagChain> tcList1 = tc1.getValue();
		LinkedList<TagChain> tcList2 = tc2.getValue();
		assert(tcList1 != null && tcList2 != null); 
		if(tcList1.isEmpty() && tcList2.isEmpty())
			return 0;
		else if(tcList1.isEmpty())
			return 1;
		else if(tcList2.isEmpty())
			return -1;
		
		if(tcList1 == null || tcList2 == null || tcList1.isEmpty() || tcList2.isEmpty())
			return 0;
		double score1 = tcList1.getFirst().getScore();
		double score2 = tcList2.getFirst().getScore();
		if(score1 > score2)
			return -1;
		else if(score1 == score2)
			return 0;
		else
			return 1;
	}
	
	public boolean equal(LinkedList<TagChain> tcList1, LinkedList<TagChain> tcList2)
	{
		assert(tcList1 != null && tcList2 != null && !tcList1.isEmpty() && !tcList2.isEmpty());
		return tcList1.getFirst().getScore() == tcList2.getFirst().getScore();
	}
	
}
