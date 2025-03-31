package scaniter;

import java.util.ArrayList;
import java.util.Collections;

import modi.Constants;
import modi.Peak;
import modi.Spectrum;
import msutil.MSMass;

//리팩토링된 MSMScan 클래스 (실 코드로 적용됨. but not verified yet)
public class MSMScan {

    static final int minPeaksCount = 4;
    static final double minMW = 8 * MSMass.getMinAAMass() + 18;

    private final String title;
    private final int specIndex;
    private final int scanNo;
    private final double pmz;
    private final double neutralMW;
    private final int charge;
    private Spectrum peaklist;  // 여기까진 동일한데

    // 인스턴스별 컨텍스트
    private ScanContext__ context;

    public MSMScan(int index, double pmz, int charge) {
        this.title = "";
        this.specIndex = index;
        this.scanNo = 0;
        this.pmz = pmz;
        this.charge = charge;
        this.neutralMW = (pmz - Constants.Proton) * charge;
        this.context = ScanContext__.fromConstants(); // 초기 컨텍스트 설정
    }

    public MSMScan(String title, int index, int sn, double pmz, int charge) {
        this.title = title;
        this.specIndex = index;
        this.scanNo = sn;
        this.pmz = pmz;
        this.charge = charge;
        this.neutralMW = (pmz - Constants.Proton) * charge;
        this.context = ScanContext__.fromConstants(); // 초기 컨텍스트 설정
    }

    public double getObservedMW() {
        return neutralMW;
    }

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

        int rangeForIsotopeIncrement = context.getRangeForIsotopeIncrement();
        if (rangeForIsotopeIncrement != 0) {
            context.setMaxNoOfC13((int) Math.ceil(neutralMW / rangeForIsotopeIncrement));
        }

        double ppmTolerance = context.getPpmTolerance();
        if (ppmTolerance != 0) context.setPrecursorAccuracy(PPMtoDalton(neutralMW, ppmTolerance));


        context.setPrecursorTolerance(context.getPrecursorAccuracy() + context.getMaxNoOfC13() * Constants.IsotopeSpace);

        int index = 0;
        Spectrum spectrum = new Spectrum(this.pmz, this.charge, this.title);

        double basePeakIntensity = 0, TIC = 0;
        double tarMass = 0, tarInten = 0;

        // 스펙트럼 데이터 처리 (변경 없음)
        final double tol = context.getMassToleranceForDenovo();
        for (RawPeak rp : rawPL) {
            double mass = rp.mz;
            double intensity = rp.it;
            if (intensity <= 0 || mass <= 0) continue;
            if (mass > neutralMW) continue;

            if ((mass - tarMass) < tol) {
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

        context.setGapTolerance(Constants.fragmentTolerance * 2);
        if(context.getPrecursorTolerance() < context.getMassToleranceForDenovo()) context.setNonModifiedDelta(context.getPrecursorTolerance());
        else context.setNonModifiedDelta(context.getMassToleranceForDenovo());

        if (context.getPrecursorTolerance() > context.getGapTolerance()) {
            double tempo = context.getGapTolerance();
            context.setGapTolerance(tempo + context.getPrecursorTolerance());
        }

        if (spectrum.size() < minPeaksCount) {
            peaklist = null;
        } else {
            peaklist = spectrum;
        }

        // 계산된 값으로 새 컨텍스트 생성
        this.context = ScanContext__.fromConstants()
                .withPrecursorTolerance(context.getPrecursorTolerance())
                .withPrecursorAccuracy(context.getPrecursorAccuracy())
                .withMaxNoOfC13(context.getMaxNoOfC13())
                .withGapTolerance(context.getGapTolerance())
                .withNonModifiedDelta(context.getNonModifiedDelta())
        // 다른 필드도 필요에 따라 업데이트
        ;

        return (peaklist != null);
    }

    private void removeReporterIons(ArrayList<RawPeak> rawPL, double[] removedMasses) {
        ArrayList<RawPeak> reporters = new ArrayList<>();
        for (int i = 1; i < removedMasses.length; i++) {
            reporters.add(new RawPeak(removedMasses[i], Constants.fragmentTolerance));
        }

        reporters.add(new RawPeak(removedMasses[0] + Constants.Proton, Constants.fragmentTolerance));

        int fragCS = 1;
        while (true) {
            double compItraqTag = (this.neutralMW - removedMasses[0] + Constants.Proton * fragCS) / fragCS;

            double secondIso = compItraqTag + Constants.IsotopeSpace / fragCS;
            double thirdIso = secondIso + Constants.IsotopeSpace / fragCS;
            double forthIso = thirdIso + Constants.IsotopeSpace / fragCS;

            reporters.add(new RawPeak(compItraqTag, Constants.fragmentTolerance));
            reporters.add(new RawPeak(secondIso, Constants.fragmentTolerance));
            reporters.add(new RawPeak(thirdIso, Constants.fragmentTolerance));
            reporters.add(new RawPeak(forthIso, Constants.fragmentTolerance));

            for (int i = 1; i <= context.getMaxNoOfC13(); i++) {
                reporters.add(new RawPeak(compItraqTag - i * Constants.IsotopeSpace / fragCS, Constants.fragmentTolerance));
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

        public Spectrum getSpectrum() {
            return spectrum;
        }

        public ScanContext__ getContext() {
            return context;
        }
    }


}
