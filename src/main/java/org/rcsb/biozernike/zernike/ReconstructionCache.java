package org.rcsb.biozernike.zernike;
import org.rcsb.biozernike.complex.Complex;

import java.io.*;
import java.util.Map;

public class ReconstructionCache {
    private static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_arr = null;

    public static boolean cacheExists(int DIM){
        return new File("Reconstruction"+DIM+".coefs").exists();
    }

    public  static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> readZpCache(int DIM) throws IOException, ClassNotFoundException {
        if (zp_arr==null) {
            System.out.println("Reading cache...");

            ObjectInputStream ois = new ObjectInputStream(new FileInputStream("Reconstruction"+DIM+".coefs"));

            Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals =
                    (Map<Integer, Map<Integer, Map<Integer, Complex[]>>>) ois.readObject();
            ois.close();

            zp_arr = zp_vals;
            System.out.println("Finished reading cache");
        }
        return zp_arr;
    }

    public static void writeZpCache(Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals, int DIM) throws Exception {

//        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream("test.coefs"));
        System.out.println("Writing cache...");
        ObjectOutputStream oos = new  ObjectOutputStream(new FileOutputStream("Reconstruction"+DIM+".coefs"));
        oos.writeObject(zp_vals);
        oos.close();
        String op = new File("Reconstruction"+DIM+".coefs").getAbsolutePath();
        System.out.println("Done writing cache "+ op);
    }

}
