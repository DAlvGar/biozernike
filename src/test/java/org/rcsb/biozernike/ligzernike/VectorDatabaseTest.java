package org.rcsb.biozernike.ligzernike;

import java.util.Random;
import org.junit.Test;

public class VectorDatabaseTest {

    public static double[] randomArray(int length) {
        double[] randomlist = new double[length];
        Random random = new Random();
        for (int i = 0; i < length; i++) {
            randomlist[i] = random.nextDouble();
        }
        return randomlist;
    }

    @Test
    public void testVectorDB() {
        double[] moments = randomArray(444);
        double[] geom = randomArray(13);
        double[] fp = randomArray(72);
        System.err.println(geom[0]+","+geom[1]+","+geom[2]+","+geom[3]+"...");
        VectorDatabaseImpl VectorDatabase = new VectorDatabaseImpl();
        VectorDatabase.createDatabase();
        VectorDatabase.insertTarget("Example Target");
        int molid = VectorDatabase.insertMolecule("Example Target", "Example Molecule",true, "smiles");
        int confid = VectorDatabase.insertConformer(molid, 2, geom, "ctab", "Originaltitle");
        int mapid = VectorDatabase.insertVolumeMap(confid, "Example Map");
        VectorDatabase.insertMoments(mapid, 1, moments);
        VectorDatabase.insertFPVector(mapid, fp);
        VectorDatabase.queryData();
        VectorDatabase.closeConnection();
    }
}
