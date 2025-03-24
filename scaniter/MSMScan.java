package scaniter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;

import modi.Constants;
import modi.Peak;
import modi.Spectrum;
import msutil.MSMass;

public class MSMScan {
	
	static final int minPeaksCount = 4;
	static final double minMW = 8 * MSMass.getMinAAMass() + 18;
	
	private final String 		title;
	private final int 		specIndex;
	private final int 		scanNo;
	private final double 		pmz;
	private final double 		neutralMW;
	private final int 		charge;
	private Spectrum 	peaklist;
	private static final double tolerance= Constants.massToleranceForDenovo;
	
	//private double 	fragmentTol = 0;
	private double 	precursorTolerance = 0;
	private double 	precursorAccuracy= 0;
	private double 	gapTolerance = 0;
	private double 	nonModifiedDelta = 0;
	private int		maxNoOfC13 = 0;
	
	public MSMScan(int index, double pmz, int charge){
		this.title 		= "";
		this.specIndex	= index;
		this.scanNo		= 0;	
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}
	
	public MSMScan(String title, int index, int sn, double pmz, int charge){		
		this.title 		= title;
		this.specIndex  = index;
		this.scanNo		= sn;
		this.pmz 		= pmz;
		this.charge 	= charge;
		this.neutralMW 	= (pmz - Constants.Proton)*charge;
	}

	public double 	getObservedMW(){ return neutralMW; }

	public String 	getHeader(){ return String.format("%d\t%.4f\t%d\t%d\t%s",
													specIndex, neutralMW, charge, scanNo, title); }
	
	public Spectrum getSpectrum() { 
		Constants.precursorTolerance= precursorTolerance;
		Constants.precursorAccuracy	= precursorAccuracy;
		Constants.gapTolerance 		= gapTolerance;
		Constants.gapAccuracy 		= precursorAccuracy + 2*Constants.fragmentTolerance;
		Constants.nonModifiedDelta 	= nonModifiedDelta;
		Constants.maxNoOfC13 		= maxNoOfC13;
		return peaklist; 
	}

	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {
		
		if( neutralMW < minMW || Constants.maxPeptideMass < neutralMW ) return false; 
				
		if( Constants.reporterMassOfIsobaricTag != null ) removeReporterIons(rawPL, Constants.reporterMassOfIsobaricTag);
		
		if( Constants.rangeForIsotopeIncrement != 0 ) maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );
		else maxNoOfC13 = Constants.maxNoOfC13;
		
		if( Constants.PPMTolerance != 0 ) precursorAccuracy = Constants.PPMtoDalton( neutralMW, Constants.PPMTolerance );
		else precursorAccuracy = Constants.precursorAccuracy;
		
		precursorTolerance = precursorAccuracy + maxNoOfC13*Constants.IsotopeSpace;
		
		int index = 0;
		Spectrum spectrum = new Spectrum( this.pmz, this.charge, this.title );
		
		double basePeakIntensity=0, TIC=0;
		double tarMass=0, tarInten=0;
		for( RawPeak rp : rawPL ) {		
			double mass = rp.mz;
			double intensity = rp.it;
			if( intensity <= 0 || mass <= 0 ) continue;
			if( mass > neutralMW ) continue;
			
			if( ( mass - tarMass ) < tolerance ){
				double sum = tarInten + intensity;
				tarMass = tarMass*(tarInten/sum)+ mass*(intensity/sum);
				tarInten += intensity;
				spectrum.get(index-1).set(tarMass, tarInten);
			}
			else{
				spectrum.add( new Peak(index++, mass, intensity) );
				tarMass = mass;
				tarInten = intensity; 
			}
			TIC += intensity;
			if( tarInten > basePeakIntensity )
				basePeakIntensity= tarInten;
		}	
		spectrum.setExtraInformation( basePeakIntensity, TIC );
		
		gapTolerance = Constants.fragmentTolerance*2;
		nonModifiedDelta = (precursorTolerance < Constants.massToleranceForDenovo)? precursorTolerance : Constants.massToleranceForDenovo;
				
		if( precursorTolerance > gapTolerance ) gapTolerance += precursorTolerance;
		
		if( spectrum.size() < minPeaksCount ) peaklist = null; 
		else peaklist = spectrum;
		
		return (peaklist!=null);
	}
	
	private void removeReporterIons( ArrayList<RawPeak> rawPL, double[] removedMasses ){
	
		ArrayList<RawPeak> reporters = new ArrayList<>();
		for(int i=1; i<removedMasses.length; i++)
			reporters.add( new RawPeak(removedMasses[i], Constants.fragmentTolerance) );
		
		reporters.add( new RawPeak(removedMasses[0] + Constants.Proton, Constants.fragmentTolerance) );

		int fragCS = 1;
		while( true ){
			double compItraqTag = (this.neutralMW - removedMasses[0] + Constants.Proton*fragCS)/fragCS;
	
			double secondIso = compItraqTag + Constants.IsotopeSpace/fragCS;
			double thirdIso  = secondIso + Constants.IsotopeSpace/fragCS;
			double forthIso  = thirdIso + Constants.IsotopeSpace/fragCS;
			
			reporters.add( new RawPeak(compItraqTag, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(secondIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(thirdIso, Constants.fragmentTolerance) );
			reporters.add( new RawPeak(forthIso, Constants.fragmentTolerance) );

			for(int i=1; i<=maxNoOfC13; i++){
				reporters.add( new RawPeak(compItraqTag-i*Constants.IsotopeSpace/fragCS, Constants.fragmentTolerance) );
			}
			if( ++fragCS >= this.charge ) break;
		}

		Collections.sort(reporters);
		
		int start = 0;
		for( RawPeak rp : reporters ){
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

	private static ArrayList<Integer> getCharge(String csline){
		ArrayList<Integer> cslist = new ArrayList<>();
		
		StringTokenizer csTok = new StringTokenizer(csline, "|");
		while( csTok.hasMoreTokens() ){
			String csStr = csTok.nextToken();
			int t_cs;
			int st = 0, ed = 1;
			for(int i=st; i<csStr.length(); i++){
				if( Character.isDigit( csStr.charAt(i) ) ){
					st = i;
					ed = st+1;
					break;
				}
			}
			for(int i=ed; i<csStr.length(); i++){
				if( !Character.isDigit( csStr.charAt(i) ) ){
					ed = i;
					break;
				}
			}
			try {
				t_cs= Integer.parseInt( csStr.substring(st,ed) );
			} catch (NumberFormatException e) {
				t_cs = 0;
			}
			if( t_cs!= 0 ) cslist.add(t_cs);
		}
		return cslist;
	}
	
}

























