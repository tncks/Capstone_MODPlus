package scaniter;

import java.util.ArrayList;
import java.util.Collections;

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


	/// /* * * * */ // problematic
	public Spectrum getSpectrum() {
		//Constants.PPMTolerance is also important factor (it has code dependency)
		Constants.precursorTolerance= precursorTolerance; /*dependency check temporally finished*/ //문제 상황 파악 어느정도 완료되었다는 뜻
		Constants.precursorAccuracy	= precursorAccuracy; /*dependency check temporally finished*/
		Constants.gapTolerance 		= gapTolerance;     /*dependency check temporally finished*/
		Constants.gapAccuracy 		= precursorAccuracy + 2*Constants.fragmentTolerance; /*dependency check temporally finished*/
		Constants.nonModifiedDelta 	= nonModifiedDelta; /*dependency check temporally finished*/
		Constants.maxNoOfC13 		= maxNoOfC13; /*dependency check temporally finished*/
		return peaklist; 
	}/*****/



	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {
		
		if( neutralMW < minMW || Constants.maxPeptideMass < neutralMW ) return false; 
				
		if( Constants.reporterMassOfIsobaricTag != null ) removeReporterIons(rawPL, Constants.reporterMassOfIsobaricTag);
		
		if( Constants.rangeForIsotopeIncrement != 0 ) maxNoOfC13 = (int)Math.ceil( neutralMW / Constants.rangeForIsotopeIncrement );
		else maxNoOfC13 = Constants.maxNoOfC13;
		
		if( Constants.PPMTolerance != 0 ) precursorAccuracy = PPMtoDalton( neutralMW, Constants.PPMTolerance );
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
		nonModifiedDelta = (precursorTolerance < Constants.massToleranceForDenovo)? precursorTolerance : Constants.massToleranceForDenovo; ///////
				
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

	public double PPMtoDalton(double mass, double ppm)
	{
		return mass/1000000*ppm;
	}

	/*private static ArrayList<Integer> getCharge(String csline){
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
	}*/
	
}



// 기존 Constants.java 코드는 어느정도 수정을 할 것인지 생각하기
// set_parameter() 코드의 의존성에도 수정 시 주의하기
// Temporal for now
// 스캔 컨텍스트 클래스 - Constants 대신 사용할 설정 정보를 담는 불변 객체
class ScanContext__ {
	private final double precursorTolerance;
	private final double precursorAccuracy;
	private final double gapTolerance;
	private final double gapAccuracy;
	private final double nonModifiedDelta;
	private final int    maxNoOfC13;  // 여기까진 필수인데  (precursorTolerance ~ maxNoOfC13) 까지

	private final double fragmentTolerance;
	private final double massToleranceForDenovo;
	private final double[] reporterMassOfIsobaricTag;
	private final double proton;
	private final double isotopeSpace;
	private final int    ppmTolerance;
	private final double rangeForIsotopeIncrement;
	private final double maxPeptideMass;

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
				(int)Constants.PPMTolerance, ///////////
				Constants.rangeForIsotopeIncrement,
				Constants.maxPeptideMass
		);
	}

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
			int ppmTolerance,  ///////////
			double rangeForIsotopeIncrement,
			double maxPeptideMass) {
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
		this.maxPeptideMass = maxPeptideMass;
	}

	// getter 메서드들
	public double getPrecursorTolerance() { return precursorTolerance; }
	public double getPrecursorAccuracy() { return precursorAccuracy; }
	public double getGapTolerance() { return gapTolerance; }
	public double getGapAccuracy() { return gapAccuracy; }
	public double getNonModifiedDelta() { return nonModifiedDelta; }
	public int    getMaxNoOfC13() { return maxNoOfC13; }
	public double getFragmentTolerance() { return fragmentTolerance; }
	public double getMassToleranceForDenovo() { return massToleranceForDenovo; }
	public double[] getReporterMassOfIsobaricTag() { return reporterMassOfIsobaricTag; }
	public double getProton() { return proton; }
	public double getIsotopeSpace() { return isotopeSpace; }
	public int    getPpmTolerance() { return ppmTolerance; }
	public double getRangeForIsotopeIncrement() { return rangeForIsotopeIncrement; }
	public double getMaxPeptideMass() { return maxPeptideMass; }

	// 새로운 값으로 업데이트된 ScanContext 반환하는 with 메서드
	public ScanContext__ withPrecursorTolerance(double precursorTolerance) {
		return new ScanContext__(
				precursorTolerance, this.precursorAccuracy, this.gapTolerance, this.gapAccuracy,
				this.nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
				this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
				this.rangeForIsotopeIncrement, this.maxPeptideMass
		);
	}

	public ScanContext__ withPrecursorAccuracy(double precursorAccuracy) {
		return new ScanContext__(
				this.precursorTolerance, precursorAccuracy, this.gapTolerance, this.gapAccuracy,
				this.nonModifiedDelta, this.maxNoOfC13, this.fragmentTolerance, this.massToleranceForDenovo,
				this.reporterMassOfIsobaricTag, this.proton, this.isotopeSpace, this.ppmTolerance,
				this.rangeForIsotopeIncrement, this.maxPeptideMass
		);
	}

}


//리팩토링된 MSMScan 클래스 (예시)
class MSMScan_new {

	static final int minPeaksCount = 4;
	static final double minMW = 8 * MSMass.getMinAAMass() + 18;

	private final String title;
	private final int specIndex;
	private final int scanNo;
	private final double pmz;
	private final double neutralMW;
	private final int charge;
	private Spectrum peaklist;  // 여기까진 동일한데
	private static final double tolerance = Constants.massToleranceForDenovo; // 이것도 컨텍스트로 이동 가능

	// 인스턴스별 컨텍스트
	private ScanContext__ context;

	public MSMScan_new(int index, double pmz, int charge) {
		this.title = "";
		this.specIndex = index;
		this.scanNo = 0;
		this.pmz = pmz;
		this.charge = charge;
		this.neutralMW = (pmz - Constants.Proton) * charge;
		this.context = ScanContext__.fromConstants(); // 초기 컨텍스트 설정
	}

	public MSMScan_new(String title, int index, int sn, double pmz, int charge) {
		this.title = title;
		this.specIndex = index;
		this.scanNo = sn;
		this.pmz = pmz;
		this.charge = charge;
		this.neutralMW = (pmz - Constants.Proton) * charge;
		this.context = ScanContext__.fromConstants(); // 초기 컨텍스트 설정
	}

	public double getObservedMW() { return neutralMW; }

	public String getHeader() {
		return String.format("%d\t%.4f\t%d\t%d\t%s", specIndex, neutralMW, charge, scanNo, title);
	}

	// Important (code changed)
	public Spectrum getSpectrum() {
		return peaklist;
	}


	// Important (code changed)
	public boolean setSpectrum(ArrayList<RawPeak> rawPL) {

		if (neutralMW < minMW || Constants.maxPeptideMass < neutralMW) return false;


		double[] reporterMassOfIsobaricTag = context.getReporterMassOfIsobaricTag();
		if (reporterMassOfIsobaricTag != null) {
			removeReporterIons(rawPL, reporterMassOfIsobaricTag);
		}


		int maxNoOfC13;
		double rangeForIsotopeIncrement = context.getRangeForIsotopeIncrement();
		if (rangeForIsotopeIncrement != 0) {
			maxNoOfC13 = (int)Math.ceil(neutralMW / rangeForIsotopeIncrement);
		} else {
			maxNoOfC13 = context.getMaxNoOfC13();
		}

		double precursorAccuracy;
		int ppmTolerance = context.getPpmTolerance();
		if (ppmTolerance != 0) {
			precursorAccuracy = PPMtoDalton(neutralMW, ppmTolerance);
		} else {
			precursorAccuracy = context.getPrecursorAccuracy();
		}

		double precursorTolerance = precursorAccuracy + maxNoOfC13 * context.getIsotopeSpace();

		int index = 0;
		Spectrum spectrum = new Spectrum(this.pmz, this.charge, this.title);

		double basePeakIntensity = 0, TIC = 0;
		double tarMass = 0, tarInten = 0;

		// 스펙트럼 데이터 처리 (변경 없음)
		for (RawPeak rp : rawPL) {
			double mass = rp.mz;
			double intensity = rp.it;
			if (intensity <= 0 || mass <= 0) continue;
			if (mass > neutralMW) continue;

			if ((mass - tarMass) < tolerance) {
				double sum = tarInten + intensity;
				tarMass = tarMass * (tarInten / sum) + mass * (intensity / sum);
				tarInten += intensity;
				spectrum.get(index - 1).set(tarMass, tarInten);
			} else {
				spectrum.add(new Peak(index++, mass, intensity));
				tarMass = mass;
				tarInten = intensity;
			}
			TIC += intensity;
			if (tarInten > basePeakIntensity) {
				basePeakIntensity = tarInten;
			}
		}
		spectrum.setExtraInformation(basePeakIntensity, TIC);

		double gapTolerance = context.getFragmentTolerance() * 2;
		double nonModifiedDelta = (precursorTolerance < context.getMassToleranceForDenovo()) ?
				precursorTolerance : context.getMassToleranceForDenovo();

		if (precursorTolerance > gapTolerance) {
			gapTolerance += precursorTolerance;
		}

		// 결과 저장 및 새 컨텍스트 생성
		if (spectrum.size() < minPeaksCount) {
			peaklist = null;
		} else {
			peaklist = spectrum;
		}

		// 계산된 값으로 새 컨텍스트 생성
		this.context = ScanContext__.fromConstants()
				.withPrecursorTolerance(precursorTolerance)
				.withPrecursorAccuracy(precursorAccuracy)
		// 다른 필드도 필요에 따라 업데이트
		;

		return (peaklist != null);
	}

	private void removeReporterIons(ArrayList<RawPeak> rawPL, double[] removedMasses) {
		ArrayList<RawPeak> reporters = new ArrayList<>();
		for (int i = 1; i < removedMasses.length; i++) {
			reporters.add(new RawPeak(removedMasses[i], context.getFragmentTolerance()));
		}

		reporters.add(new RawPeak(removedMasses[0] + context.getProton(), context.getFragmentTolerance()));

		int fragCS = 1;
		while (true) {
			double compItraqTag = (this.neutralMW - removedMasses[0] + context.getProton() * fragCS) / fragCS;

			double secondIso = compItraqTag + context.getIsotopeSpace() / fragCS;
			double thirdIso = secondIso + context.getIsotopeSpace() / fragCS;
			double forthIso = thirdIso + context.getIsotopeSpace() / fragCS;

			reporters.add(new RawPeak(compItraqTag, context.getFragmentTolerance()));
			reporters.add(new RawPeak(secondIso, context.getFragmentTolerance()));
			reporters.add(new RawPeak(thirdIso, context.getFragmentTolerance()));
			reporters.add(new RawPeak(forthIso, context.getFragmentTolerance()));

			for (int i = 1; i <= context.getMaxNoOfC13(); i++) {
				reporters.add(new RawPeak(compItraqTag - i * context.getIsotopeSpace() / fragCS, context.getFragmentTolerance()));
			}
			if (++fragCS >= this.charge) break;
		}

		Collections.sort(reporters);

		int start = 0;
		for (RawPeak rp : reporters) {
			for (int i = start; i < rawPL.size(); i++) {
				if (rawPL.get(i).mz < rp.mz - rp.it) continue;
				else if (rawPL.get(i).mz > rp.mz + rp.it) {
					start = i;
					break;
				}
				rawPL.remove(i);
				i--;
			}
		}
	}

	public double PPMtoDalton(double mass, double ppm) {
		return mass / 1000000 * ppm;
	}



	//

	// 스펙트럼과 함께 컨텍스트 정보도 반환하는 새 메서드
	public SpectrumWithContext getSpectrumWithContext() {
		return new SpectrumWithContext(peaklist, context);
	}

	//
	//
	// 컨텍스트와 스펙트럼을 함께 포함하는 클래스
	public static class SpectrumWithContext {
		private final Spectrum spectrum;
		private final ScanContext__ context;

		public SpectrumWithContext(Spectrum spectrum, ScanContext__ context) {
			this.spectrum = spectrum;
			this.context = context;
		}

		public Spectrum getSpectrum() { return spectrum; }
		public ScanContext__ getContext() { return context; }
	}


}



// ScanContext 클래스 -> basic 한 형태. record 스타일.
// No usage -> future usage, on board stage remaining
// Temporal (inspecting in progress)
final class ScanContext {
	private  double precursorTolerance;
	private  double precursorAccuracy;
	private  double gapTolerance;
	private  double gapAccuracy;
	private  double nonModifiedDelta;
	private  int    maxNoOfC13;
	private  double maxPeptideMass;
	private  double fragmentTolerance;






	public ScanContext(double precursorTolerance, double precursorAccuracy, double gapTolerance, double gapAccuracy, double nonModifiedDelta, int maxNoOfC13,
					   double maxPeptideMass, double fragmentTolerance
	) {
		this.precursorTolerance = precursorTolerance;
		this.precursorAccuracy = precursorAccuracy;
		this.gapTolerance = gapTolerance;
		this.gapAccuracy = gapAccuracy;
		this.nonModifiedDelta = nonModifiedDelta;
		this.maxNoOfC13 = maxNoOfC13;
		this.maxPeptideMass = maxPeptideMass;
		this.fragmentTolerance = fragmentTolerance;
	}



	public double getPrecursorTolerance() {
		return precursorTolerance;
	}

	public void setPrecursorTolerance(double precursorTolerance) {
		this.precursorTolerance = precursorTolerance;
	}

	public double getMaxNoOfC13() {
		return maxNoOfC13;
	}

	public void setMaxNoOfC13(int maxNoOfC13) {
		this.maxNoOfC13 = maxNoOfC13;
	}

	public double getMaxPeptideMass() {
		return maxPeptideMass;
	}

	public void setMaxPeptideMass(double maxPeptideMass) {
		this.maxPeptideMass = maxPeptideMass;
	}

	public double getFragmentTolerance() {
		return fragmentTolerance;
	}

	public void setFragmentTolerance(double fragmentTolerance) {
		this.fragmentTolerance = fragmentTolerance;
	}



	public double getPrecursorAccuracy() {
		return precursorAccuracy;
	}



	public double getGapTolerance() {
		return gapTolerance;
	}



	public double getGapAccuracy() {
		return gapAccuracy;
	}



	public double getNonModifiedDelta() {
		return nonModifiedDelta;
	}

	public void setPrecursorAccuracy(double precursorAccuracy) {
		this.precursorAccuracy = precursorAccuracy;
	}

	public void setGapTolerance(double gapTolerance) {
		this.gapTolerance = gapTolerance;
	}

	public void setGapAccuracy(double gapAccuracy) {
		this.gapAccuracy = gapAccuracy;
	}

	public void setNonModifiedDelta(double nonModifiedDelta) {
		this.nonModifiedDelta = nonModifiedDelta;
	}
}
//