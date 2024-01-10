package org.rcsb.biozernike;

import org.junit.Test;
import org.rcsb.biozernike.ligzernike.RocksDBInterface;


public class RocksDBInterfaceTest {
    @Test
    public void writeRocksDB() {
        String dbPath = "rocksdb_test.db"; // Change this to your desired database path
        RocksDBInterface dbInterface = new RocksDBInterface(dbPath);

        // Adding content to the database
        dbInterface.put("key1", "value1");
        dbInterface.put("key2", "value2");

        // Retrieving content from the database
        String value1 = dbInterface.get("key1");
        String value2 = dbInterface.get("key2");

        System.out.println("Retrieved values:");
        System.out.println("key1: " + value1);
        System.out.println("key2: " + value2);

        // Don't forget to close the database when done
        dbInterface.close();
    }

    @Test
    public void writeMomentsRocksDB() {
        String dbPath = "rocksdb_test.db"; // Change this to your desired database path
        RocksDBInterface dbInterface = new RocksDBInterface(dbPath);

        // Adding content to the database
        dbInterface.put("key1", "value1");
        dbInterface.put("key2", "value2");

        // Retrieving content from the database
        String value1 = dbInterface.get("key1");
        String value2 = dbInterface.get("key2");

        System.out.println("Retrieved values:");
        System.out.println("key1: " + value1);
        System.out.println("key2: " + value2);

        // Don't forget to close the database when done
        dbInterface.close();
    }

    @Test
    public void readRocksDB() {
        String dbPath = "rocksdb_test.db"; // Change this to your desired database path
        RocksDBInterface dbInterface = new RocksDBInterface(dbPath);

        // Retrieving content from the database
        String value1 = dbInterface.get("key1");
        String value2 = dbInterface.get("key2");

        System.out.println("Retrieved values:");
        System.out.println("key1: " + value1);
        System.out.println("key2: " + value2);

        // Don't forget to close the database when done
        dbInterface.close();
    }

}
