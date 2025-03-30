package modi;

import moda.DPHeap;
import moda.DPPeptide;
import processedDB.HeatedDB;
import processedDB.ProtDatabase;
import processedDB.StemTagTrie;

public class DBHolder {
    public DBHolder() {}

    public HeatedDB getHeatedDB(StemTagTrie stemDB, DPHeap candidates, DPHeap tepids, int numHeatedPeptides_f) {
        HeatedDB matchedBits = new HeatedDB();
        int count = 0;
        for (DPPeptide dp : candidates) {
            if (dp.getScore() < 1) break;
            String modapept = dp.getPeptide();
            int pro_start = dp.getProtein();
            ProtDatabase proDB = stemDB.get(dp.getStem());
            matchedBits.add(proDB.getProteinIdentity(pro_start), dp.getStem(), pro_start, pro_start + modapept.length());
            if (++count == numHeatedPeptides_f) break;
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
}
