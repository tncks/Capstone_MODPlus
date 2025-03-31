package modi;

import java.util.ArrayList;
import java.util.Collections;

public class ScanCap implements Comparable<ScanCap> {
	private final String 	title;
	private final double 	pmz;
	private final double 	neutralMW;
	private final int 	charge;
	private int 	scanNo;
	private long 	offset;
	
	private static final double tolerance= Constants.massToleranceForDenovo;
	
	public ScanCap(String title, double pmz, int charge){
	
		this.title 		= title;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public ScanCap(String title, int sn, double pmz, int charge){		
		this.title 		= title;	
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public void setOffset(long offset){ this.offset = offset; }
	public String getTitle(){ return title; }
	public int getScanNumber(){ return scanNo; }
	public double getObservedMW(){ return neutralMW; }
	public double getPMZ(){ return pmz; }
	public int getCharge(){ return charge; }
	public long getOffset(){ return offset; }
	

	
	private static class RawP implements Comparable<RawP> {
		final double mz;
		final double it;
		public RawP(double m, double i){
			mz=m;
			it=i;
		}	
		public int compareTo(RawP p) {
			if( mz > p.mz ) return 1;
			else if( mz < p.mz ) return -1;
			else return 0;
		}
	}
	
	public int compareTo(ScanCap s) 
	{
		if( this.neutralMW > s.neutralMW ) return 1;
		else if( this.neutralMW < s.neutralMW ) return -1;
		
		if( this.charge > s.charge ) return 1;
		else if( this.charge < s.charge ) return -1;
		else return 0;
	}
	
	private void xprocessingiTRAQ(ArrayList<RawP> rawPL){
		
		ArrayList<RawP> reporters = new ArrayList<>();
		reporters.add( new RawP(114.1105, Constants.fragmentTolerance) );//reporterIon114
		reporters.add( new RawP(115.1074, Constants.fragmentTolerance) );//reporterIon115
		reporters.add( new RawP(116.1107, Constants.fragmentTolerance) );//reporterIon116
		reporters.add( new RawP(117.1141, Constants.fragmentTolerance) );//reporterIon117
		
		reporters.add( new RawP(Constants.NTERM_FIX_MOD + Constants.Proton, Constants.fragmentTolerance) );//iTRAQ TAG
		
		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - Constants.NTERM_FIX_MOD + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawP(compItraqTag, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(secondIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(thirdIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			reporters.add( new RawP(forthIso, Constants.fragmentTolerance) );//precursor without iTRAQ TAG//*/

//			reporters.add( new RawP(compItraqTag, (3*Constants.IsotopeSpace/fragCS)+Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			for(int i=1; i<=Constants.maxNoOfC13; i++){
				reporters.add( new RawP(compItraqTag-i*Constants.IsotopeSpace/fragCS, Constants.fragmentTolerance) );//precursor without iTRAQ TAG
			}//*/

			if( ++fragCS >= this.charge ) break;
		}

		Collections.sort(reporters);
		
		int start = 0;
		for( RawP rp : reporters ){
			for (int i=start; i<rawPL.size(); i++){
				if( rawPL.get(i).mz < rp.mz-rp.it ) continue;
				else if( rawPL.get(i).mz > rp.mz+rp.it ) {
					start = i;
					break;
				}				
				rawPL.remove(i);
				i--;
			}
		}
	
	}
}
