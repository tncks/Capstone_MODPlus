package moda;

import java.util.LinkedList;

import modi.Constants;
import msutil.*;

public class DPHeap extends LinkedList<DPPeptide> {


    private static final int Capacity = 100; // max heap capacity //MODplus

    public DPHeap() {
        for (int i = 0; i < Capacity; i++)
            this.add(new DPPeptide());
    }

    public void setStemNo(int stem) {
        for (DPPeptide dp : this)
            dp.stem = stem;
    }

    public boolean isConfident() {
        int cnt = 0;
        for (DPPeptide dp : this) {
            if (dp.isConfident()) return true;
            if (++cnt >= 10) break;
        }
        return false;
    }

    public boolean insert(DPPeptide dp) {

        if (dp.score < 1 || dp.compareTo(this.getLast()) == 1) return false;

        int i = 0;
        for (DPPeptide x : this) {
            if (dp.isSame(x)) break;
            if (x.compareTo(dp) == 1) {
                this.add(i, dp);
                this.removeLast();
                break;
            }
            i++;
        }
        return true;
    }


    public boolean insertInternal(DPPeptide dp) {
        if (dp.score < 1 || dp.compareTo(this.getLast()) == 1) return false;

        int i = 0;
        for (DPPeptide x : this) {
            if (dp.isSame(x)) break;
            if (x.compareTo(dp) == 1) {
                this.add(i, dp);
                this.removeLast();
                break;
            }
            i++;
        }
        return true;
    }

    public void insertAll(DPHeap heap) {

        heap.forEach(this::insertInternal);
    }


    public int evaluate(PGraph graph) { //for moda
        int i, maxScore = 0;
        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score < 1) {
                this.remove(i);
                i--;
                continue;
            }
            this.get(i).evaluatePSM(graph);
            if (maxScore < this.get(i).score) maxScore = this.get(i).score;
        }


        final int _maxScore = maxScore;
        this.forEach(dp -> {
            if (dp.score == _maxScore) dp.score -= 1;
        });

        this.sort(new DPPeptideRefinedComparator());
        if (this.size() == 0 || this.get(0).score < 1) return 0;
        return i;
    }

    public int reArr(PGraph graph) throws NullPointerException { //for modplus
        int i;
        int maxScore = 0;
        DPPeptide ele;
        for (i = 0; i < this.size(); i++) {
            ele = this.get(i);
            if (ele.score < 1) {
                this.remove(i);
                i--;
                continue;
            }
            this.get(i).score = getModARankScore(ele.peptide, ele.ptms, graph);
            ele = this.get(i);
            if (maxScore < ele.score) maxScore = ele.score;
        }

        this.sort(new DPPeptideRefinedComparator());
        int cut = maxScore / 2;
        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score < cut) this.removeRange(i, this.size());
        }

        if (this.size() == 0 || this.get(0).score < 1) return 0;
        return i;
    }

    //MODa Scoring
    public int getModARankScore(String peptide, double[] ptms, PGraph graph) { //for modA preliminary ranking
        IonGraph iGraph;
        if (Constants.INSTRUMENT_TYPE == Constants.msms_type.QTOF) iGraph = new TOFGraph(peptide, ptms, graph);
        else iGraph = new TRAPGraph(peptide, ptms, graph);

        iGraph.setScore(graph);
        return iGraph.getRankScore();
    }

}













