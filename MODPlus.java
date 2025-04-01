import java.io.*;
import java.text.DateFormat;
import java.util.*;

import moda.DPPeptide;
import modi.*;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import moda.DPHeap;
import moda.MultiMOD;
import moda.OneMOD;
import msutil.IsobaricTag;
import msutil.MSMass;
import msutil.PGraph;
import msutil.ProtCutter;
import processedDB.*;
import scaniter.MSMScan;
import scaniter.ScanIterator;

public class MODPlus {
    static boolean dynamicPMCorrection = false, multiBlind = true;
    static int numHeatedPeptides = 50;
    private static final String[] message = {
            "[Error] Cannot read any MS/MS scan from input dataset.\r\n" +
                    "[Error] Check consistency between input file and its format.",

            "[Error] Cannot read any protein from input database.\r\n" +
                    "[Error] Check input fasta format.",

            "[Error] One fixed modification per amino acid can be allowed.\r\n" +
                    "[Error] Check specfied fixed modifications.",

            "[Error] Unsupported character set in your search parameter",

            "[Error] Required field is empty.\r\n" +
                    "[Error] Required fields : MS/MS Data, Database",

            "[Error] Wrong usage.\r\n" +
                    "[Error] Re-confirm it.",

            "[Error] Not defined"
    };


    public static void main(String[] args) throws Exception {
        Constants.engine = "modplus";
        Constants.engineVersion = "hyu";


        int availableCores = Runtime.getRuntime().availableProcessors();


        System.out.println("************************************************************************************");
        System.out.println("Modplus (version " + Constants.engineVersion + ") - Identification of post-translational modifications");
        System.out.println("Release Date: 2025");
        System.out.println("Available CPU Cores: " + availableCores);
        System.out.println("************************************************************************************");
        System.out.println();


        run(args[0]);
    }


    protected static int set_parameter(String Prixparam) throws Exception {

        System.out.println("Reading parameter.....");
        System.out.println("DEBUG: cpu-AvailableProcessors: " + Runtime.getRuntime().availableProcessors());
        System.out.println();

        Document doc;
        try {
            doc = new SAXBuilder().build(Prixparam);
        } catch (JDOMException e) {
            System.out.println(message[3]);
            return 1;
        } catch (IOException e) {
            System.out.println(message[5]);
            return 5;
        }

        Element search = doc.getRootElement();
        Constants.runDate = DateFormat.getDateInstance().format(new Date());
        if (search.getAttributeValue("user") != null) {
            Constants.runUser = search.getAttributeValue("user");
        }
        if (search.getAttributeValue("title") != null) {
            Constants.runTitle = search.getAttributeValue("title");
        } else Constants.runTitle = String.valueOf(System.currentTimeMillis());

        Element dataset = search.getChild("dataset");
        if (dataset != null) {
            Constants.SPECTRUM_LOCAL_PATH = dataset.getAttributeValue("local_path");
            if (Constants.SPECTRUM_LOCAL_PATH == "") {
                System.out.println(message[4]);
                return 4;
            }

            String type = dataset.getAttributeValue("format");
            if (type.compareToIgnoreCase("mgf") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MGF;
            else if (type.compareToIgnoreCase("pkl") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.PKL;
            else if (type.compareToIgnoreCase("ms2") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MS2;
            else if (type.compareToIgnoreCase("dta") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.DTA;
            else if (type.compareToIgnoreCase("mzxml") == 0)
                Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MZXML;
            else if (type.compareToIgnoreCase("zip") == 0)
                Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.ZIPDTA;

            Constants.INSTRUMENT_NAME = dataset.getAttributeValue("instrument");
            if (Constants.INSTRUMENT_NAME.equals("QTOF")) Constants.INSTRUMENT_TYPE = Constants.msms_type.QTOF;
            else Constants.INSTRUMENT_TYPE = Constants.msms_type.TRAP;
        }
        System.out.print("Input datasest : " + Constants.SPECTRUM_LOCAL_PATH);
        System.out.println(" (" + Constants.SPECTRA_FILE_TYPE + " type)");

        Element database = search.getChild("database");
        if (database != null) {
            Constants.PROTEIN_DB_LOCAL_PATH = database.getAttributeValue("local_path");
            if (Constants.PROTEIN_DB_LOCAL_PATH == "") {
                System.out.println(message[4]);
                return 4;
            }
        }
        System.out.println("Input database : " + Constants.PROTEIN_DB_LOCAL_PATH);

        Element enzyme = search.getChild("enzyme");//DEPRECATED
        if (enzyme != null) {
            String enzymeName = enzyme.getAttributeValue("name");
            String cut = enzyme.getAttributeValue("cut");
            String sence = enzyme.getAttributeValue("sence");
            Mutables.protease = new ProtCutter(enzymeName, cut, sence);
        }//*/

        Element com_enzyme = search.getChild("combined_enzyme");
        if (com_enzyme != null) {
            String enzymeName = com_enzyme.getAttributeValue("name");
            String nn = com_enzyme.getAttributeValue("nterm_cleave");
            String cc = com_enzyme.getAttributeValue("cterm_cleave");
            Mutables.protease = new ProtCutter(enzymeName, nn, cc, true);
        }

        Mutables.alkylatedToCys = 0; //DEPRECATED
        Element cys_alkylated = search.getChild("cys_alkylated");
        if (cys_alkylated != null) {
            Mutables.alkylationMethod = cys_alkylated.getAttributeValue("name");
            Mutables.alkylatedToCys = Double.parseDouble(cys_alkylated.getAttributeValue("massdiff"));
            AminoAcid.modifiedAminoAcidMass('C', Mutables.alkylatedToCys);
            MSMass.modifiedAminoAcidMass('C', Mutables.alkylatedToCys);
        }

        Element instrument_resolution = search.getChild("instrument_resolution");
        if (instrument_resolution != null) {
            Constants.MSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("ms")) == 0) ? 1 : 0;
            if (Constants.MSResolution == 1) System.out.println("High resolution MS!!");
            Constants.MSMSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("msms")) == 0) ? 1 : 0;
            if (Constants.MSMSResolution == 1) System.out.println("High resolution MS/MS!!");
        }

        Element parameters = search.getChild("parameters");
        if (parameters != null) {
            Element param;
            param = parameters.getChild("enzyme_constraint");
            if (param != null) {
                Mutables.missCleavages = Integer.parseInt(param.getAttributeValue("max_miss_cleavages"));
                Mutables.numberOfEnzymaticTermini = Integer.parseInt(param.getAttributeValue("min_number_termini"));
                if (Mutables.numberOfEnzymaticTermini > 2) Mutables.numberOfEnzymaticTermini = 2;
            }

            param = parameters.getChild("isotope_error");
            if (param != null) {
                if (param.getAttributeValue("min_C13_number") != null)
                    Mutables.minNoOfC13 = Integer.parseInt(param.getAttributeValue("min_C13_number"));

                if (param.getAttributeValue("max_C13_number") != null)
                    Mutables.maxNoOfC13 = Integer.parseInt(param.getAttributeValue("max_C13_number"));

                if (Mutables.maxNoOfC13 == 0 && param.getAttributeValue("increment_per_dalton") != null)
                    Mutables.rangeForIsotopeIncrement = Integer.parseInt(param.getAttributeValue("increment_per_dalton"));
            }

            param = parameters.getChild("peptide_mass_tol");
            if (param != null) {
                if (param.getAttributeValue("unit").compareToIgnoreCase("ppm") == 0) {
                    Mutables.PPMTolerance = Double.parseDouble(param.getAttributeValue("value"));
                } else {
                    Mutables.precursorTolerance = Mutables.precursorAccuracy = Double.parseDouble(param.getAttributeValue("value"));
                }
            }

            param = parameters.getChild("fragment_ion_tol");
            if (param != null) {
                Mutables.fragmentTolerance = Double.parseDouble(param.getAttributeValue("value"));
            }
            param = parameters.getChild("modified_mass_range");
            if (param != null) {
                Mutables.minModifiedMass = Double.parseDouble(param.getAttributeValue("min_value"));
                Mutables.maxModifiedMass = Double.parseDouble(param.getAttributeValue("max_value"));
            }
        }

        Element protocol = search.getChild("protocol");
        if (protocol != null) {
            System.out.print("Protocol Description: ");
            Element isobaric = protocol.getChild("isobaric_labeling");
            if (isobaric != null) {
                if (isobaric.getAttributeValue("reagent") != null) {
                    Mutables.isobaricTag = isobaric.getAttributeValue("reagent");
                    Mutables.reporterMassOfIsobaricTag = IsobaricTag.getReporterMasses(isobaric.getAttributeValue("reagent"));
                    if ("".compareTo(Mutables.isobaricTag) != 0)
                        System.out.print(Mutables.isobaricTag + " Labelled" + ((Mutables.reporterMassOfIsobaricTag == null) ? " (NOT Supported)" : " (Supported)"));
                }
            }
            Element modEnrich = protocol.getChild("modification_enrichment");
            if (modEnrich != null) {
                if (modEnrich.getAttributeValue("mod") != null) {
                    Mutables.enrichedModification = modEnrich.getAttributeValue("mod");
                    if ("".compareTo(Mutables.enrichedModification) != 0)
                        System.out.print(" & " + Mutables.enrichedModification + " Enriched" +
                                (("Acetyl".compareToIgnoreCase(Mutables.enrichedModification) == 0 || "Phospho".compareToIgnoreCase(Mutables.enrichedModification) == 0) ? " (Supported)" : " (NOT Supported)"));
                }
            }
            System.out.println();
        }

        Element modifications = search.getChild("modifications");
        Mutables.variableModifications = new PTMDB();
        Mutables.fixedModifications = new PTMDB();

        if (modifications != null) {
            double[] fixedAA = new double[26];

            Element fixed = modifications.getChild("fixed");
            if (fixed != null) {
                if (Mutables.fixedModifications.setFixedModificatinos(fixed, fixedAA) == 0) {
                    System.out.println(message[2]);
                    return 2;
                }
            }
            if (!Mutables.fixedModifications.isEmpty())
                System.out.println("Fixed modifications : " + Mutables.fixedModifications.size() + " selected");

            Element variable = modifications.getChild("variable");
            if (variable != null) {
                Mutables.PTM_FILE_NAME = variable.getAttributeValue("local_path");
                boolean canBeModifiedOnFixedAA = variable.getAttributeValue("canBeModifiedOnFixedAA").equals("1");
                Mutables.canBeModifiedOnFixedAA = canBeModifiedOnFixedAA;
                if (Mutables.PTM_FILE_NAME != null) {
                    Mutables.variableModifications.setVariableModificatinos(Mutables.PTM_FILE_NAME, fixedAA, canBeModifiedOnFixedAA);
                }
                Mutables.variableModifications.setVariableModificatinos(variable, fixedAA, canBeModifiedOnFixedAA);

                if (canBeModifiedOnFixedAA) {
                    for (PTM p : Mutables.fixedModifications) {
                        Mutables.variableModifications.add(
                                new PTM(Mutables.variableModifications.size(), "De-" + p.getName(), "",
                                        -p.getMassDifference(), 0, p.getResidue(), p.getPTMPosition(), (p.getAbbAA() == 'C') ? 1 : 0));
                    }
                }
                if (variable.getAttributeValue("multi_mods") != null && variable.getAttributeValue("multi_mods").equals("0")) {
                    Mutables.maxPTMPerGap = Mutables.maxPTMPerPeptide = 1;
                }
            }
            if (!Mutables.variableModifications.isEmpty()) {
                System.out.print("Variable modifications : " + Mutables.variableModifications.size() + " selected (");
                Mutables.variableModifications.setPTMDiagnosticIon();
                if (Mutables.maxPTMPerPeptide == 1) System.out.println("one modification per peptide)");
                else System.out.println("multiple modifications per peptide)");
            }
        }
        Mutables.variableModifications.constructPTMLookupTable();

        Element decoy_search = search.getChild("decoy_search");
        if (decoy_search != null) {
            if (decoy_search.getAttributeValue("checked") != null) {
                if ("1".compareTo(decoy_search.getAttributeValue("checked")) == 0) {
                    Mutables.targetDecoy = 1;
                    System.out.println("Decoy search checked");
                }
            }
        }

        Element multistages_search = search.getChild("multistages_search");
        if (multistages_search != null) {
            if (multistages_search.getAttributeValue("checked") != null) {
                if ("1".compareTo(multistages_search.getAttributeValue("checked")) == 0) {
                    Mutables.firstSearchProgram = multistages_search.getAttributeValue("program");
                    System.out.println("MultiStages Search checked " + Mutables.firstSearchProgram);
                }
            }
        }

        Element mod_map = search.getChild("mod_map");
        if (mod_map != null) {
            if (mod_map.getAttributeValue("checked") != null) {
                if ("1".compareTo(mod_map.getAttributeValue("checked")) == 0) {
                    System.out.println("MODMap checked");
                }
            }
        }

        Constants.adjustParameters();

        System.out.println();
        return 0;
    }


    public static void run(String arg) throws Exception {


        try {
            if (set_parameter(arg) != 0) return;
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        try {
            Constants.MAX_TAG_SIZE = 100;
            Constants.minTagLength = 2;
            Constants.minTagLengthPeptideShouldContain = 3;
            Constants.tagChainPruningRate = 0.4;
            File analPath = new File(Constants.SPECTRUM_LOCAL_PATH);
            if (analPath.isDirectory()) {
                String type = Constants.SPECTRA_FILE_TYPE.toString().toLowerCase();
                for (File file : analPath.listFiles()) {
                    if (file.getName().endsWith(type)) {
                        Constants.SPECTRUM_LOCAL_PATH = file.getPath();
                        System.out.println("Input dataset: " + Constants.SPECTRUM_LOCAL_PATH);
                        modplus_mod_search();
                    }
                }
                System.out.println("End of process");
            } else {
                modplus_mod_search();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    static int modplus_mod_search() throws Exception{
        System.out.println("Starting MODPlus for modification search!");
        long startTime= System.currentTimeMillis();

        ScanIterator scaniter = ScanIterator.get( Constants.SPECTRUM_LOCAL_PATH, Constants.SPECTRA_FILE_TYPE );
        StemTagTrie  ixPDB = new StemTagTrie( Constants.PROTEIN_DB_LOCAL_PATH );
        final boolean considerIsotopeErr = ( Mutables.maxNoOfC13 !=0 || Mutables.precursorTolerance > 0.50001 )? true : false;
        if( scaniter == null || scaniter.size() == 0 || ixPDB.getSizeOfEntries() == 0 ){
            return 1;
        }
        System.out.println();

        int index = 1;
        int iterSize= scaniter.size();
        String identifier = Constants.SPECTRUM_LOCAL_PATH;
        identifier = identifier.substring(0, identifier.lastIndexOf('.'));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(identifier+".modplus.txt")));
        while( scaniter.hasNext() ){

            ArrayList<MSMScan> chargedSpectra = scaniter.getNext();
            System.out.println("MODPlus | " + (index++) + "/" + iterSize);

            int selected = -1;
            ArrayList<AnsPeptide> candidates = null;

            for(int i=0; i<chargedSpectra.size(); i++){
                Spectrum spectrum = chargedSpectra.get(i).getSpectrum();
                PGraph graph = spectrum.getPeakGraph();
                spectrum.setCorrectedParentMW( graph.correctMW( dynamicPMCorrection ) );
                TagPool tPool = null;

                if(spectrum != null) {

                    for(int t=0; t<spectrum.size(); t++){
                        if( Math.abs(spectrum.get(t).getMass() - spectrum.getPrecursor() ) < 2 ) {
                            spectrum.remove(t); ///////////////////////////////////////////////////////////
                            t--;
                            continue;
                        }
                        spectrum.get(t).setIndex(t);
                    }

                    spectrum.normalizeIntensityLocally();

                    int extra = ( spectrum.getCharge() > 2 && Constants.INSTRUMENT_TYPE != Constants.msms_type.QTOF )? 2 : 0;
                    spectrum.peakSelection(Constants.selectionWindowSize, Mutables.minNumOfPeaksInWindow+extra );

                    tPool = spectrum.generateTags(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain, Mutables.massToleranceForDenovo);

                }


                DPHeap heatedPepts = OneMOD.getHeatedPeptides( ixPDB, graph, tPool, considerIsotopeErr );
                DPHeap tepidPepts  = null;
                if( Mutables.maxPTMPerPeptide > 1 )
                    if( heatedPepts == null || !heatedPepts.isConfident() ) {
                        tepidPepts = heatedPepts;
                        heatedPepts = MultiMOD.getHeatedPeptides( ixPDB, graph, tPool, dynamicPMCorrection );
                    }

                if( heatedPepts == null ) continue;

                HeatedDB bitDB =  getHeatedDB( ixPDB, heatedPepts, tepidPepts );
                TagTrie bitTrie = bitDB.getPartialDB(ixPDB);

                ArrayList<AnsPeptide> tp= dynamicMODeye( bitTrie, graph, tPool );
                if( tp.size() > 0 ) {
                    if( candidates == null || candidates.get(0).compareTo( tp.get(0) ) == 1 ) {
                        candidates = tp;
                        selected = i;
                    }
                }
            }

            if( selected != -1 ) {
                MSMScan scan = chargedSpectra.get(selected);

                HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap = new HashMap<>();
                out.println( ">>"+scaniter.getFileName()+"\t"+scan.getHeader() );

                for( int k=0; k<candidates.size(); k++ ){
                    String tpSeq = candidates.get(k).getPeptideSequence();
                    ArrayList<PeptideMatchToProtein> matchedProteins = seqToProtMap.get(tpSeq);

                    if( matchedProteins == null ){
                        matchedProteins = ixPDB.getMatchProteins(tpSeq);
                        seqToProtMap.put(tpSeq, matchedProteins);
                    }
                    out.println( candidates.get(k).toMODPlus(scan.getObservedMW(), matchedProteins) );
                }
                out.println();
            }
        }
        out.close();

        System.out.println("[MOD-Plus] Elapsed Time : " + (System.currentTimeMillis()-startTime)/1000 + " Sec" );
        return 0;
    }
    private static class ResultEntry {
        final MSMScan scan;
        final ArrayList<AnsPeptide> candidates;
        final HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap;

        ResultEntry(MSMScan scan, ArrayList<AnsPeptide> candidates,
                    HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap) {
            this.scan = scan;
            this.candidates = candidates;
            this.seqToProtMap = seqToProtMap;
        }
    }

    private static class SpectrumAnalyzer {

        private SpectrumAnalyzer() {
        }


        TagChainPool buildTagChain(MatchedTagPool matchedTags) {
            TagChainPool tagChainPool = new TagChainPool();
            tagChainPool.buildTagChainPool(matchedTags);
            return tagChainPool;
        }

        boolean interpretTagChain(PTMDB ptmDB, TagChainPool tcPool, PGraph graph) {
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

        MatchedTagPool extendedBuildMatchedTagPool(TagPool primitiveTags, double motherMass,
                                                   TagTrie ixPDB, ProtCutter enzyme, int NTT) {
            if (primitiveTags == null || ixPDB == null)
                return null;

            double minDelta = (Mutables.minModifiedMass < 0) ? Mutables.minModifiedMass - Mutables.gapTolerance : -Mutables.gapTolerance;
            double maxDelta = (Mutables.maxModifiedMass > 0) ? Mutables.maxModifiedMass + Mutables.gapTolerance : +Mutables.gapTolerance;
            TagPool longTags = primitiveTags.extractAbove(Constants.minTagLengthPeptideShouldContain);

            int realTag = 0;
            double orbMass = motherMass - Constants.H2O;
            RetrivedPeptideMap searchResults = new RetrivedPeptideMap();
            for (Tag tag : longTags) {

                RetrivedPeptideMap bRes = ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, tag.getBIonNtermOffset() - Mutables.NTERM_FIX_MOD, tag,
                        tag.getBIonCtermOffset() - Mutables.CTERM_FIX_MOD, IonDirection.B_DIRECTION, minDelta, maxDelta, Mutables.gapTolerance);
                searchResults.combine(bRes);

                Tag reverseTag = tag.reverseTag();
                RetrivedPeptideMap yRes = ixPDB.getRetrivedPeptides(orbMass, enzyme, NTT, reverseTag.getYIonNtermOffset() - Mutables.NTERM_FIX_MOD, reverseTag,
                        reverseTag.getYIonCtermOffset() - Mutables.CTERM_FIX_MOD, IonDirection.Y_DIRECTION, minDelta, maxDelta, Mutables.gapTolerance);
                searchResults.combine(yRes);
                realTag++;

                if (realTag > Constants.MAX_TAG_SIZE * 2) break;
            }
            return searchResults.convertToMatchedTagPool(primitiveTags.extract(Constants.minTagLength, Constants.minTagLengthPeptideShouldContain));
        }

    }


    private static HeatedDB getHeatedDB(StemTagTrie stemDB, DPHeap candidates, DPHeap tepids) {
        HeatedDB matchedBits = new HeatedDB();
        int count = 0;
        for (DPPeptide dp : candidates) {
            if (dp.getScore() < 1) break;
            String modapept = dp.getPeptide();
            int pro_start = dp.getProtein();
            ProtDatabase proDB = stemDB.get(dp.getStem());
            matchedBits.add(proDB.getProteinIdentity(pro_start), dp.getStem(), pro_start, pro_start + modapept.length());
            if (++count == numHeatedPeptides) break;
        }

        count = 0;
        if (tepids != null) {
            for (DPPeptide dp : tepids) {
                if (dp.getScore() < 1) break;
                String modapept = dp.getPeptide();
                int pro_start = dp.getProtein();
                ProtDatabase proDB = stemDB.get(dp.getStem());
                matchedBits.add(proDB.getProteinIdentity(pro_start), dp.getStem(), pro_start, pro_start + modapept.length());
                if (++count == 10) break;
            }
        }
        return matchedBits;
    }


    private static ArrayList<AnsPeptide> dynamicMODeye(TagTrie dynamicDB, PGraph graph, TagPool tPool) throws Exception {
        SpectrumAnalyzer szer = new SpectrumAnalyzer();
        MatchedTagPool matchedList = szer.extendedBuildMatchedTagPool(tPool, graph.getCorrectedMW(),
                dynamicDB, Mutables.protease, Mutables.numberOfEnzymaticTermini);

        TagChainPool tcPool = new TagChainPool();
        tcPool.putAll(szer.buildTagChain(matchedList));
        tcPool.discardPoorTagChain();

        boolean specAnnotated = false;
        if (tcPool.size() != 0) {
            specAnnotated = szer.interpretTagChain(Mutables.variableModifications, tcPool, graph);
        }

        ArrayList<AnsPeptide> cands = new ArrayList<>();
        if (tcPool.size() != 0 && specAnnotated) {
            cands = tcPool.getAnswerPeptides(graph);
        }
        return cands;
    }

}





