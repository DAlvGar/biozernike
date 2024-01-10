package org.rcsb.biozernike.ligzernike;

import org.rcsb.biozernike.zernike.ZernikeMomentsIO;
import org.rocksdb.Options;
import org.rocksdb.RocksDB;
import org.rocksdb.RocksDBException;

import java.util.List;
import java.util.stream.Collectors;

import org.rcsb.biozernike.complex.Complex;

public class RocksDBInterface {
    private RocksDB db;

    public RocksDBInterface(String dbPath) {
        RocksDB.loadLibrary(); // Ensure RocksDB native library is loaded

        Options options = new Options()
                .setCreateIfMissing(true)
                .setKeepLogFileNum(1) // Set to keep only one log file
                .setMaxLogFileSize(100 * 1024 * 1024); // Set maximum log file size to 100MB

        try {
            db = RocksDB.open(options, dbPath);
        } catch (RocksDBException e) {
            e.printStackTrace();
        }
    }

    public void put(String key, String value) {
        try {
            byte[] byteKey = key.getBytes();
            byte[] byteValue = value.getBytes();
            db.put(byteKey, byteValue);
        } catch (RocksDBException e) {
            e.printStackTrace();
        }
    }

    public void putMoments(String key, List<Complex> complexList) {
        put(key, ZernikeMomentsIO.complexListToString(complexList));
    }

    public void putFP(String key, List<Double> fpList) {
        String doubleConcat = fpList.stream().map(String::valueOf).collect(Collectors.joining(","));
        put(key, doubleConcat);
    }

    public String get(String key) {
        try {
            byte[] byteKey = key.getBytes();
            byte[] byteValue = db.get(byteKey);
            return byteValue != null ? new String(byteValue) : null;
        } catch (RocksDBException e) {
            e.printStackTrace();
            return null;
        }
    }

    public List<Complex> getMoments(String key) {
        return ZernikeMomentsIO.readComplex(get(key));
    }

    public List<Double> getFP(String key) {
        return ZernikeMomentsIO.readDoubleFromString(get(key));
    }

    public void close() {
        if (db != null) {
            db.close();
        }
    }

}