package modi;

import scaniter.ScanContext__;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeSet;

public class TagChain extends TreeSet<SpecInterpretation> implements Comparable<TagChain> {
	final Peptide		matchedPeptide;
	final Spectrum	sourceSpectrum;
	double 		score = 0;
	boolean		allGapAnnotated = true;
	double 		mostAbundantBYPeakIntensity;
	int 		longestTagLen= 0; 
	int 		tagHit= 0;
	int 		tagCoverage= 0;

	public Peptide	getMatchedPeptide()			{ return matchedPeptide; }
	
	public	TagChain(Peptide matchedPeptide, Spectrum sourceSpectrum)
	{
		this.matchedPeptide = matchedPeptide; 
		this.sourceSpectrum = sourceSpectrum;
	}

	public Spectrum getSourceSpectrum() {
		return this.sourceSpectrum;
	}
	
	public	void	setAllGapAnnotated(boolean annotated)	{ allGapAnnotated = annotated; }

	public	int	 compareTo(TagChain tc)
	{ 
		double score1 = this.getScore();
		double score2 = tc.getScore();
		
		if(score1 < score2)
			return 1;
		else if(score1 == score2)
			return 0;
		else
			return -1;
	}

	public void setMostAbundantBYPeakIntensity()
	{
		ArrayList<Peak> theoPeaks = new ArrayList<>();
		for(SpecInterpretation t : this)
			theoPeaks.addAll(t.getTheoreticalPeaks());
		PeakListComparator comparator = new PeakListComparator(sourceSpectrum, theoPeaks);
		for(PeakPair pp : comparator.getSharedPeaks())
		{
			if(pp.getSecond().getPeakProperty() == PeakProperty.B_ION ||
					pp.getSecond().getPeakProperty() == PeakProperty.Y_ION)
			{
				if(this.mostAbundantBYPeakIntensity <  pp.getFirst().getIntensity())
					this.mostAbundantBYPeakIntensity = pp.getFirst().getIntensity();
			}
		}
	}
	
	public boolean makeGap(ScanContext__ context) {
		int peptLastIndex = matchedPeptide.size()-1;

		setMostAbundantBYPeakIntensity();
		
		// check if matchedPeptide is Protein N-Term/C-term : should move to Peptide
		boolean protNTerm = false, protCTerm = false;
		
		if( matchedPeptide.getSrcProteinInfo().getFirst().getSrcProteinID() == 0 ) protNTerm = true;
		else if( matchedPeptide.getSrcProteinInfo().getFirst().getSrcProteinID() == 2 ) protCTerm = true;

		int start = 0, end;
		ArrayList<Gap> gapList = new ArrayList<>();
		PTMPosition position;
		double prevOffset = Constants.NTERM_FIX_MOD;
        double bStart = Constants.B_ION_OFFSET;
		double yStart;
		boolean hasQTag = false;
		for(SpecInterpretation element : this)
		{
			assert(element instanceof MatchedTag);
			MatchedTag curTag = (MatchedTag)element;
			if( curTag.size() > Constants.minTagLengthPeptideShouldContain ) hasQTag = true;
			double motherMass = curTag.getSourceSpectrum().getCorrectedMW();
			end = element.getStart()-1;
			
			if(curTag.getDirection() == IonDirection.B_DIRECTION)
				yStart = curTag.getFirst().getComplementMass(motherMass);
			else
				yStart = curTag.getLast().getMass();
			
			if(start <= end)
			{
				// determine gap position
				if(start == 0 && protNTerm )
					position = PTMPosition.PROTEIN_N_TERM;
				else if(start == 0 && !protNTerm )
					position = PTMPosition.ANY_N_TERM;
				else
					position = PTMPosition.ANYWHERE;
				
				double offset = curTag.getNTermOffset() - prevOffset;
				
				Gap tpGap = new Gap(matchedPeptide, start, end, bStart, yStart, offset, 
						position, sourceSpectrum, this.mostAbundantBYPeakIntensity);
				gapList.add(tpGap);
				
				prevOffset = curTag.getNTermOffsetByLastPeak();
            }
			if(curTag.getDirection() == IonDirection.B_DIRECTION)
				bStart = curTag.getLast().getMass();
			else
				bStart = curTag.getFirst().getComplementMass(curTag.getSourceSpectrum().getCorrectedMW());
			start = element.getEnd() + 1;
		}
		if( !hasQTag ) return false;
		
		// make C-term gap 
		if(start <= peptLastIndex)
		{
			end = peptLastIndex;
			if(protCTerm)
				position = PTMPosition.PROTEIN_C_TERM;
			else
				position = PTMPosition.ANY_C_TERM;
			
			Gap tpGap = new Gap(matchedPeptide, start, end, bStart, Constants.Y_ION_OFFSET, 
					((MatchedTag)this.last()).getCTermOffset()-Constants.CTERM_FIX_MOD, position, sourceSpectrum, this.mostAbundantBYPeakIntensity);
			gapList.add( tpGap );
		}
		
		for(Gap gap : gapList){	
			if( !gap.isReasonable(context) )
				return false;
			this.add(gap);
		}
		return true;
	}

	public int round(double a){
		if( a > 0 ) return (int)(a + 0.5);
		else return (int)(a - 0.5);
	}

	public void setTagChainScore(ScanContext__ context) {
		
		HashSet<Integer> modtype = new HashSet<>();

        this.score= 0;
		for( SpecInterpretation t : this ){
			if(t instanceof Gap gap) {
                if( Math.abs( gap.getOffset() ) > context.getGapTolerance() ) {
					modtype.add( round(gap.getOffset()) );
				}
			}
			if(t instanceof MatchedTag tag){
				this.score -= 1;

                tagCoverage += tag.sequence().size();
							
				for( Peak p: tag ) this.score += p.getProbability();			
			}
		}
		
		this.score -= modtype.size();
		
		if( this.tagHit == 1 && this.longestTagLen < Constants.minTagLengthPeptideShouldContain ){
			if(modtype.isEmpty()) {
				this.score = 0; //because this is unmodified
			}
		}//*/
	}
	
	public double getScore(){		// should be optimized later
		return score;
	}
	

	
	public String getTagChainSequence()
	{
		StringBuffer output = new StringBuffer();
		for(SpecInterpretation t : this)
			output.append(t.getSequenceStr());

		return output.toString();
	}
	




	public boolean useSamePeak(MatchedTag t) {
		for(SpecInterpretation element : this)
		{
			if(!(element instanceof MatchedTag)) continue;
			if( ((MatchedTag)element).useSamePeak(t) ) {
				return true;
			}
		}
		return false;
	}
	
	public boolean add(SpecInterpretation si){
	//	System.out.println(si);
		if( si instanceof MatchedTag ){
			tagHit++;
			if( longestTagLen < ((MatchedTag)si).sequence().size() )
				longestTagLen= ((MatchedTag)si).sequence().size();
		}			
		return super.add(si);
	}
	
	public boolean remove(SpecInterpretation si){
		if( si instanceof MatchedTag ){	tagHit--; }
		return super.remove(si);
	}
	

	
	
	public ArrayList<PTMCombination> getPTMCombination(){
		
		ArrayList<PTMCombination> answers = new ArrayList<>();
		ArrayList<Gap> gapList = new ArrayList<>();
		
		int madeSize = 1;
		for( SpecInterpretation t : this ){		
			if(t instanceof Gap gap){
                if(gap.getGapInterpretation().isEmpty()) continue;
				gapList.add(gap);
				madeSize *= gap.getGapInterpretation().size();
			}
		}
		
		if(gapList.isEmpty()){
			answers.add( new PTMCombination( matchedPeptide.size() ) );
			return answers;
		}
		
		int poolSize = 1;
		int[] indexTrace = new int[gapList.size()];
		while( true ){

			PTMCombination ptmComb = new PTMCombination( matchedPeptide.size() );
			StringBuffer tempComb = new StringBuffer();
			int ix = 0, modCount = 0;
			for( Gap gap : gapList ){
				PTMRun prun = gap.getGapInterpretation().get(indexTrace[ix]);
				modCount += prun.size();
				tempComb.append( prun.toString(matchedPeptide, gap.getStart()) );
				prun.setPTMMass(ptmComb.ptms, gap.getStart());
				prun.setPTMs(ptmComb.ptmList, gap.getStart());
				ix++;
			}

			// 안 바꾸어도 됨
			if( modCount <= Constants.maxPTMPerPeptide ) {
				ptmComb.ptmComb = tempComb.toString();
				answers.add(ptmComb);
			}
			if( poolSize++ ==  madeSize ) break;	
	
			indexTrace[--ix]++;
			while( ix > -1 && indexTrace[ix] == gapList.get(ix).getGapInterpretation().size() ){
				indexTrace[ix--] = 0;
				indexTrace[ix]++;
			}
		}
	
		return answers;
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer(this.getTagChainSequence());
	
		output.append("(");
		for(SpecInterpretation t : this)
		{
			if(t instanceof Gap)
				output.append(new DecimalFormat("#.###").format(((Gap)t).getOffset())+ ",");
		}
		output.deleteCharAt(output.length()-1);
		output.append(") " + getScore());
		
		for (SpecInterpretation si : this)
		{
			if (!(si instanceof Gap gap)) continue;
            if (gap.getGapInterpretation()==null || gap.getGapInterpretation().isEmpty()) continue;
			output.append("\nGap interpretation for "+gap.getStart()+"~"+gap.getEnd()+":");
			output.append(gap.getGapInterpretation());
		}
		
		return output.toString();
	}
}
