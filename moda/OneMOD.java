package moda;

import processedDB.*;

import msutil.MSMass;
import msutil.PGraph;
import msutil.PRM;
import msutil.PRMforHighCharge;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

import modi.Constants;
import modi.TagPool;
import scaniter.ScanContext__;

public class OneMOD {

    private static AtomicInteger bestOnlineScore = new AtomicInteger(2);

    public OneMOD(){}

    public DPHeap getHeatedPeptides(StemTagTrie stemDB, PGraph graph, TagPool tPool, boolean dynamicPMCorrection, ScanContext__ context) {
        bestOnlineScore.set(2);
        DPHeap annotation = null;
        final double MW = graph.getCorrectedMW();

        for (TagTrie stem : stemDB) {
            CandidateContainer cpool = (new DBSearch()).construct_onemod_cpool(tPool, MW, stem, context);

            DPHeap sanno = null;
            if (dynamicPMCorrection)
                sanno = run_dynamic_mass_mode(cpool, graph, stem,context);
            else sanno = run_static_mass_mode(cpool, graph, stem,context);

            if (annotation == null) annotation = sanno;
            else annotation.insertAll(sanno);
        }

        if (annotation.evaluate(graph, context) < 1) return null;

        return annotation;
    }









    public DPHeap run_static_mass_mode(CandidateContainer cpool, PGraph graph, TagTrie ixPDB, ScanContext__ context) {
        int poolsize = cpool.size();
        if (poolsize == 0) {
            DPHeap emptyList = new DPHeap();
            emptyList.setStemNo(ixPDB.getStemNo());
            return emptyList;
        }

        PRM prmTable = (graph.getCharge() > 2)
                ? new PRMforHighCharge(graph)
                : new PRM(graph);


        DPHeap topList = new DPHeap();
        MODPeptide[] processArray = Arrays.copyOf(cpool.getList(), poolsize);


        final int len = processArray.length;
        for (int i = 0; i < len; i++) {
            static_single_align(processArray[i], prmTable, ixPDB, topList,context);
        }

        topList.setStemNo(ixPDB.getStemNo());
        return topList;
    }

    public DPHeap run_dynamic_mass_mode(CandidateContainer cpool, PGraph graph, TagTrie ixPDB, ScanContext__ context) {
        int poolsize = cpool.size();
        if (poolsize == 0) {
            DPHeap emptyList = new DPHeap();
            emptyList.setStemNo(ixPDB.getStemNo());
            return emptyList;
        }

        PRM prmTable = (graph.getCharge() > 2)
                ? new PRMforHighCharge(graph)
                : new PRM(graph);


        DPHeap topList = new DPHeap();
        MODPeptide[] processArray = Arrays.copyOf(cpool.getList(), poolsize);


        final int len = processArray.length;
        for (int i = 0; i < len; i++) {
            dynamic_single_align(processArray[i], prmTable, ixPDB, topList,context);
        }

        topList.setStemNo(ixPDB.getStemNo());
        return topList;
    }

    public void static_single_align(MODPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList, ScanContext__ context) {
        double observedMass = prmTable.getPeptMass();
        String peptide = entry.getPeptide(ixPDB);


        int rowMax = 2, colMax = peptide.length() + 1;
        MatCell[][] specMatrix = new MatCell[rowMax][colMax];

        for (int n = 0; n < colMax; n++) {
            specMatrix[0][n] = new MatCell(2);
            specMatrix[1][n] = new MatCell(0);
        }

        DPPeptide temp = null;
        double delta = observedMass - MSMass.getPepMass(peptide) - MODaConst.TERMINALMOD;
        double nTermDeletion = 0.0, cTermDeletion = 0.0;
        int npi = 0, cpi = colMax;
        for (int i = entry.getStart(); i <= entry.getLeft(); i++) {
            cTermDeletion = 0.0;
            cpi = colMax;
            for (int j = entry.getEnd(); j >= entry.getRight(); j--) {

                int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease);
                if (noOfET >= Constants.numberOfEnzymaticTermini) {
                    double massRange = delta + nTermDeletion + cTermDeletion;

                    // context.getMaxModifiedMass() context.getMinModifiedMass() context.getPrecursorTolerance()
                    if ((massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) ||
                            Math.abs(massRange) < Constants.precursorTolerance) {
                        double ptm = MODaConst.ptmUnit.getPtmMass(massRange);
                        int intptm = Constants.round(massRange);
                        double pmzErr = massRange - ptm;

                        specMatrix[0][npi].setMass(0, 0., 0);
                        specMatrix[0][npi].refresh();
                        specMatrix[1][npi].setMass(0, ptm, intptm);
                        specMatrix[1][npi].refresh();

                        double cellMass = Constants.NTERM_FIX_MOD;
                        for (int n = npi + 1; n < cpi; n++) {
                            cellMass += MSMass.getAAMass(peptide.charAt(n - 1));
                            specMatrix[0][n].setMass(cellMass, 0., 0);
                            specMatrix[0][n].refresh();
                            specMatrix[1][n].setMass(cellMass, ptm, intptm);
                            specMatrix[1][n].refresh();
                        }

                        specMatrix[0][cpi - 1].mass += Constants.CTERM_FIX_MOD;
                        specMatrix[1][cpi - 1].mass += Constants.CTERM_FIX_MOD;

                        temp = dynamic_programming(peptide.substring(npi, cpi - 1), observedMass - pmzErr, rowMax, npi, cpi,
                                specMatrix, prmTable, pmzErr,context);
                        if (temp != null) {
                            temp.setProteinAndNTT(entry.getStart(), noOfET);
                            topList.insert(temp);
                        }
                    }
                }
                cpi--;
                if (cpi - npi < 4) break;
                cTermDeletion += MSMass.getAAMass(peptide.charAt(cpi - 1));
            }
            nTermDeletion += MSMass.getAAMass(peptide.charAt(npi++));
        }
    }

    public void dynamic_single_align(MODPeptide entry, PRM prmTable, TagTrie ixPDB, DPHeap topList, ScanContext__ context) {
        double observedMass = prmTable.getPeptMass();
        String peptide = entry.getPeptide(ixPDB);


        int rowMax = 2, colMax = peptide.length() + 1;
        MatCell[][] specMatrix = new MatCell[rowMax][colMax];

        for (int n = 0; n < colMax; n++) {
            specMatrix[0][n] = new MatCell(2);
            specMatrix[1][n] = new MatCell(0);
        }

        DPPeptide temp = null;
        double delta = observedMass - MSMass.getPepMass(peptide) - MODaConst.TERMINALMOD;
        double nTermDeletion = 0., cTermDeletion = 0.;
        int npi = 0, cpi = colMax;

        for (int i = entry.getStart(); i <= entry.getLeft(); i++) {
            cTermDeletion = 0.;
            cpi = colMax;
            for (int j = entry.getEnd(); j >= entry.getRight(); j--) {

                int noOfET = ixPDB.getNTTOfPeptide(i, j, Constants.protease);
                if (noOfET >= Constants.numberOfEnzymaticTermini) {
                    double massRange = delta + nTermDeletion + cTermDeletion;
                    // context.getMaxModifiedMass() context.getMinModifiedMass() context.getPrecursorTolerance()
                    if ((massRange < Constants.maxModifiedMass && massRange > Constants.minModifiedMass) ||
                            Math.abs(massRange) < Constants.precursorTolerance) {
                        double ptm = MODaConst.ptmUnit.getPtmMass(massRange);
                        int intptm = Constants.round(massRange);  // round ()
                        double pmzErr = massRange - ptm;
                        ptm += MODaConst.maxIsotopeError;
                        intptm += MODaConst.maxIntIsotopeError;

                        specMatrix[0][npi].setMass(0, 0., 0); //init start-columns
                        specMatrix[0][npi].refresh();
                        specMatrix[1][npi].setMass(0, ptm, intptm);
                        specMatrix[1][npi].refresh();

                        double cellMass = Constants.NTERM_FIX_MOD;
                        for (int n = npi + 1; n < cpi; n++) {
                            cellMass += MSMass.getAAMass(peptide.charAt(n - 1));
                            specMatrix[0][n].setMass(cellMass, 0., 0);
                            specMatrix[0][n].refresh();
                            specMatrix[1][n].setMass(cellMass, ptm, intptm);
                            specMatrix[1][n].refresh();
                        }
                        specMatrix[0][cpi - 1].mass += Constants.CTERM_FIX_MOD;
                        specMatrix[1][cpi - 1].mass += Constants.CTERM_FIX_MOD;

                        temp = DPwithMassCorrection(peptide.substring(npi, cpi - 1), observedMass - pmzErr, rowMax, npi, cpi,
                                specMatrix, prmTable, pmzErr,context);
                        temp.setProteinAndNTT(entry.getStart(), noOfET);
                        topList.insert(temp);
                    }
                }
                cpi--;
                if (cpi - npi < 4) break;
                cTermDeletion += MSMass.getAAMass(peptide.charAt(cpi - 1));
            }
            nTermDeletion += MSMass.getAAMass(peptide.charAt(npi++));
        }
    }

    public DPPeptide DPwithMassCorrection(String peptide, double obsMass, int rowMax, int smStart, int smEnd,
                                                  MatCell[][] specMatrix, PRM prmTable, double pmzErr, ScanContext__ context) {

        DPPeptide best = new DPPeptide(), temp = null;

        double massCorrection = MODaConst.maxIsotopeError;
        for (int i = 0; i < MODaConst.isotopePointsToBeCorrected; i++) {
            temp = dynamic_programming(peptide, obsMass + massCorrection, rowMax, smStart, smEnd,
                    specMatrix, prmTable, pmzErr - massCorrection,context);
            if (temp != null && temp.score > best.score) {
                best = temp;
            }

            massCorrection += MODaConst.isotopeUnit;
            for (int n = smStart; n < smEnd; n++) {
                specMatrix[0][n].refresh();
                specMatrix[1][n].refresh();
                specMatrix[1][n].correctMass(MODaConst.isotopeUnit);
            }
        }
        return best;
    }

    public DPPeptide dynamic_programming(String peptide, double obsMass, int rowMax, int smStart, int smEnd,
                                                 MatCell[][] specMatrix, PRM prmTable, double pmzErr, ScanContext__ context) {

        int colMax = smEnd - smStart, endingTag = 1;
        // context.getMaxModifiedMass()
        if (specMatrix[endingTag][smStart].nominalDelta > Constants.maxModifiedMass) return null;

        double upperLimit = obsMass + Constants.fragmentTolerance;
        specMatrix[0][smStart].isAAJump = 1;

        MatCell currNode, prevNode;
        for (int n = smStart + 1; n < smEnd; n++) {
            for (int m = 0; m < rowMax; m += endingTag) {

                currNode = specMatrix[m][n];
                if (currNode.mass > upperLimit) continue;

                double max = MODaConst.baseScore;
                for (int d = 0; d < rowMax; d += endingTag) {

                    prevNode = specMatrix[d][n - 1];
                    if (prevNode.isAAJump == -1) continue;

                    if (m == d) {
                        if (max <= prevNode.score) {
                            max = prevNode.score;
                            currNode.setParent(prevNode);
                            currNode.isAAJump = 1;
                        }
                    } else if (m != 0) {

                        if (currNode.mass - prevNode.mass < MODaConst.minimumDistance) continue;

                        if (max < prevNode.score) {
                            max = prevNode.score;
                            currNode.setParent(prevNode);
                            currNode.isAAJump = 0;
                        }
                    }
                }
                if (max < 0) continue;
                currNode.score = max + prmTable.getScore(currNode.mass, pmzErr);
            }
        }

        MatCell initNode = specMatrix[0][smStart], tarNode = specMatrix[endingTag][smEnd - 1];
        double idScore = (tarNode.nominalDelta == 0) ? tarNode.score : tarNode.score - Constants.rNorm[0];
        if (idScore < bestOnlineScore.get() / 2) return null;

        double[] ptms = new double[colMax - 1];
        double[] matchedList = new double[colMax];
        int ixx = colMax - 1;
        while (tarNode != initNode) {
            matchedList[ixx--] = tarNode.mass;
            if (tarNode.nominalDelta != tarNode.parent.nominalDelta) {
                ptms[ixx] = tarNode.delta - tarNode.parent.delta;
            }
            tarNode = tarNode.parent;
        }

        //Solve Symmetric path problem
        int forward = 1, backward = colMax - 2;
        double PMCorr = obsMass + Constants.H2O;
        while (forward < backward) {
            double symmetric = matchedList[forward] + matchedList[backward];
            if (Constants.fEqual(symmetric, PMCorr)) {
                idScore -= prmTable.getScore(matchedList[forward], pmzErr);
                forward++;
                backward--;
            } else if (symmetric > PMCorr)
                backward--;
            else forward++;
        }
        if (idScore > bestOnlineScore.get()) bestOnlineScore.set((int) idScore);
        return new DPPeptide(peptide, (int) idScore, ptms, smStart);
    }

}


