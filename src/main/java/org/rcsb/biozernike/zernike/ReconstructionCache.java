package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;

import java.io.*;
import java.util.Map;

public class ReconstructionCache {
    private static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_arr = null;

    public static String getCachePath(int DIM) {
        String currentUsersHomeDir = System.getProperty("user.home");
        return currentUsersHomeDir + File.separator + "LigZernike-CacheReconstruction-" + DIM + ".coefs";

    }

    public static boolean cacheExists(int DIM) {
        return new File(getCachePath(DIM)).exists();
    }

    public static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> readZpCache(int DIM)
            throws IOException, ClassNotFoundException {
        if (zp_arr == null) {
            System.out.println("Reading cache...");

            ObjectInputStream ois = new ObjectInputStream(new FileInputStream(getCachePath(DIM)));

            Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals = (Map<Integer, Map<Integer, Map<Integer, Complex[]>>>) ois
                    .readObject();
            ois.close();

            zp_arr = zp_vals;
            System.out.println("Finished reading cache");
        }
        return zp_arr;
    }

    public static void writeZpCache(Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals, int DIM)
            throws Exception {

        // ObjectOutputStream oos = new ObjectOutputStream(new
        // FileOutputStream("test.coefs"));
        System.out.println("Writing cache...");
        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(getCachePath(DIM)));
        oos.writeObject(zp_vals);
        oos.close();
        System.out.println("Done writing cache " + getCachePath(DIM));
    }

}
