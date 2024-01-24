package org.rcsb.biozernike.ligzernike;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.time.DurationFormatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.jena.sparql.function.library.leviathan.radiansToDegrees;
import org.rcsb.biozernike.AlignmentResult;
import org.rcsb.biozernike.InvariantNorm;
import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.molecules.SDFReader;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.zernike.ZernikeMoments;
import org.rcsb.biozernike.zernike.ZernikeMomentsIO;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.forester.development.neTest;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

/**
 * Hello world!
 *
 */
@Command(name = "LigZernike", mixinStandardHelpOptions = true, description = "Tool to work with volumetric files and 3D zernikes for ligands", version = "1.0")
public class LigZernike {
    private static final double[] PERCENTILES_FOR_GEOM = { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0 };
    private static final int[] INVARIANT_ORDERS = { 2, 4, 6 };

    public Volume loadMoleculeVolume(String filename, double minCap, double maxCap, boolean flip, boolean normalize,
            double multiplier)
            throws Exception {
        Volume v1 = OpenDXIO.read(filename);
        prepareVolume(v1, minCap, maxCap, flip, normalize, multiplier);
        return v1;
    }

    public String getMatrixString(Matrix4d m, String name) {
        double[] row = new double[4];
        String matrixString = name + " = [";
        for (int i = 0; i < 4; i++) {
            m.getRow(i, row);
            for (double v : row) {
                matrixString += v + ",";
            }
        }
        matrixString += "];\n";
        return matrixString;
    }

    public void writeStructurePDB(Structure s, String fileName) throws FileNotFoundException {
        try (PrintWriter out = new PrintWriter(fileName)) {
            out.println(s.toPDB());
        }
    }

    public Structure loadPDB(String filename) throws Exception {
        PDBFileReader pdbreader = new PDBFileReader();
        Structure structure = pdbreader.getStructure(filename);
        return structure;
    }

    public static void prepareVolume(Volume v1, double minCap, double maxCap, boolean flip, boolean normalize,
            double multiplier) {
        if (flip)
            v1.flipValues(); // Turn negative to positive and viceversa
        v1.capMax(maxCap); // throw away positive values, cap to zero and keep negatives
        v1.capMin(minCap);
        if (normalize)
            v1.applyContourAndNormalize(minCap, multiplier); // Eliminate values below 0.2 and normalize
        v1.updateCenter();
    }

    public static void calMomentsAndStore_rocksdb(String key, Volume volume, int order, RocksDBInterface momentsDB,
            RocksDBInterface FPDB) {
        InvariantNorm n = new InvariantNorm(volume, order);
        List<Double> fp = n.getFingerprint();
        List<Complex> moments = ZernikeMoments
                .flattenMomentsComplex(n.getMoments().getOriginalMoments());
        momentsDB.putMoments(key, moments);
        FPDB.putFP(key, fp);
    }

    public static void calMomentsAndStore_psql(Integer conformerID, String mapName, Volume volume, Integer order,
            VectorDatabaseImpl db) {
        InvariantNorm n = new InvariantNorm(volume, order);
        List<Double> fp = n.getFingerprint();
        Double i = fp.get(1);

        if (Double.isNaN(i) || Double.isInfinite(i) || i == null) {
            return; // Moments don't exist, skip volumes insertion
        }

        // Store map in conformer
        int volumeID = db.getVolumeId(conformerID, mapName);
        if (volumeID == -1) {
            volumeID = db.insertVolumeMap(conformerID, mapName);
            db.insertFPVector(volumeID, listToArray(fp));
        }

        for (int invOrder : INVARIANT_ORDERS) {
            List<Double> imoments = n.getInvariants(invOrder);
            i = imoments.get(1);
            if (Double.isNaN(i) || Double.isInfinite(i) || i == null) {
                continue;
            }
            // System.out.println("Moments: " + imoments.get(0) + "," + imoments.get(1));
            db.insertMoments(volumeID, invOrder, listToArray(imoments));
        }
    }

    public static void calMomentsAndStore_txt(String key, Volume volume, int order, String outPrefix) {
        String dest = outPrefix + key;
        InvariantNorm n = new InvariantNorm(volume, order);
        List<Double> fp = n.getFingerprint();
        try {
            n.getMoments().write(dest + "_moments.txt", false);
            ZernikeMomentsIO.writeDouble(dest + "_fp.txt", fp);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private static double[] listToArray(List<Double> i) {
        return i.stream()
                .mapToDouble(d -> d)
                .toArray();
    }

    /**
     * Convert molecules in sdf file to zernike moments
     * 
     * @param maps      - Maps to generate {int column: string mapname}. Map
     *                  'hydroele' will be split in positive and negative
     * @param sdfPath   - SDF file with atomic logp values assigned in the V block
     *                  section. each column will identify a map. Generated by
     *                  pharmscreen.
     * @param momentsDB - rocksdb interface to store the full moments set (complex
     *                  numbered), indexed by a key generated like:
     *                  moleculename_mapname
     * @param FPDB      - rocksdb interface to store invariant fingerprint
     *                  (doubles), same indexing as moments.
     * @param order     - Zernike order to use
     * @param nThreads  - number of threads to use. this still does not escale
     *                  linearly though, but some gain in speed can be achieved with
     *                  multiple threads
     */
    public static void process_sdf_(HashMap<Integer, String> maps, String sdfPath, RocksDBInterface momentsDB,
            RocksDBInterface FPDB, int order, int nThreads, boolean writeFiles, boolean printDX, String outPrefix) {

        try (IteratingSDFReader reader = SDFReader.read(sdfPath)) {
            ExecutorService executor = Executors.newFixedThreadPool(nThreads);
            FieldCalculator calculator = new FieldCalculator();
            while (reader.hasNext()) {
                IAtomContainer molecule = reader.next();
                int[] fieldIDs = maps.keySet().stream().mapToInt(i -> i).toArray();
                int NTypes = fieldIDs.length;

                Runnable task = () -> {

                    List<Volume> volumeList = calculator.projectMultiField(molecule, fieldIDs);
                    for (int i = 0; i < NTypes; i++) {
                        String m = maps.get(fieldIDs[i]);
                        String key = molecule.getTitle() + "_" + m;
                        Volume volume = volumeList.get(i);
                        if (printDX) {
                            try {
                                OpenDXIO.write(outPrefix + key + ".dx", volume);
                            } catch (IOException e) {
                                // TODO Auto-generated catch block
                                e.printStackTrace();
                            }
                        }
                        if (m.contains("hydroele")) {
                            Volume pos = new Volume(volume); // Copy
                            // NEGATIVE PART and POSITIVE part apart
                            prepareVolume(volume, 0, 20, true, true, 1.0);
                            prepareVolume(pos, 0, 20, false, true, 1.0);
                            if (writeFiles) {
                                calMomentsAndStore_txt(key, volume, order, outPrefix);
                                calMomentsAndStore_txt(key + "_pos", pos, order, outPrefix);
                            } else {
                                calMomentsAndStore_rocksdb(key, volume, order, momentsDB, FPDB);
                                calMomentsAndStore_rocksdb(key + "_pos", pos, order, momentsDB, FPDB);
                            }
                        } else if (m.contains("hydro")) {
                            prepareVolume(volume, 0, 20, true, true, 1.0);
                            if (writeFiles) {
                                calMomentsAndStore_txt(key, volume, order, outPrefix);
                            } else {
                                calMomentsAndStore_rocksdb(key, volume, order, momentsDB, FPDB);
                            }
                        } else {
                            // hbond maps
                            prepareVolume(volume, 0, 100, false, false, 10);
                            if (writeFiles) {
                                calMomentsAndStore_txt(key, volume, order, outPrefix);
                            } else {
                                calMomentsAndStore_rocksdb(key, volume, order, momentsDB, FPDB);
                            }
                        }
                    }
                };

                executor.submit(task);
            }
            // Close the SDF reader
            reader.close();

            executor.shutdown();
            try {
                executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double getRadius(Point3d[] reprPoints) {

        // Find the minimum and maximum values along all points
        double minX = Double.MAX_VALUE;
        double minY = Double.MAX_VALUE;
        double minZ = Double.MAX_VALUE;
        double maxX = Double.MIN_VALUE;
        double maxY = Double.MIN_VALUE;
        double maxZ = Double.MIN_VALUE;

        for (Point3d point : reprPoints) {
            minX = Math.min(minX, point.x);
            minY = Math.min(minY, point.y);
            minZ = Math.min(minZ, point.z);
            maxX = Math.max(maxX, point.x);
            maxY = Math.max(maxY, point.y);
            maxZ = Math.max(maxZ, point.z);
        }

        return sqrt(pow((maxX - minX), 2) + pow((maxY - minY), 2) + pow((maxZ - minZ), 2)) / 2.;
    }

    private static double[] calcGeometricDescriptor(IAtomContainer molecule) {
        Point3d[] reprPoints = SDFReader.getAtomPositions(molecule);

        List<Double> geomDescriptorList = new ArrayList<>();
        geomDescriptorList.add(getRadius(reprPoints));

        // geomDescriptorList.add(volume.getResiduesNominalWeight());

        Point3d centerPoint = new Point3d(0, 0, 0);

        for (Point3d selPoint : reprPoints) {
            centerPoint.add(selPoint);
        }
        centerPoint.scale(1 / (double) reprPoints.length);

        double[] distances = new double[reprPoints.length];
        for (int iPoint = 0; iPoint < reprPoints.length; iPoint++) {
            distances[iPoint] = centerPoint.distance(reprPoints[iPoint]);
        }

        StandardDeviation standardDeviation = new StandardDeviation();
        Skewness skewness = new Skewness();
        Kurtosis kurtosis = new Kurtosis();
        Percentile percentile = new Percentile();
        percentile.setData(distances);

        for (double p : PERCENTILES_FOR_GEOM) {
            geomDescriptorList.add(percentile.evaluate(p));
        }

        geomDescriptorList.add(standardDeviation.evaluate(distances));
        geomDescriptorList.add(skewness.evaluate(distances));
        geomDescriptorList.add(kurtosis.evaluate(distances));

        return listToArray(geomDescriptorList);
    }

    private static class RegisterTask extends ArrayList<IAtomContainer> implements Runnable {
        private final ConcurrentHashMap<IAtomContainer, String> map;
        private final HashMap<Integer, String> volmaps;
        private final VectorDatabaseImpl db;
        private final Integer targetId;
        private final FieldCalculator calculator = new FieldCalculator();
        private final Integer order;
        private final int[] fieldIDs;
        private final int NTypes;

        public RegisterTask(ConcurrentHashMap<IAtomContainer, String> map, 
                            HashMap<Integer, String> volmaps,
                            VectorDatabaseImpl db, Integer targetId, Integer order) {
            this.map = map;
            this.db = db;
            this.targetId = targetId;
            this.volmaps = volmaps;
            this.order=order;
            this.fieldIDs = this.volmaps.keySet().stream().mapToInt(i -> i).toArray();
            this.NTypes = this.fieldIDs.length;
        }

        @Override
        public void run() {
            try {
                for (IAtomContainer m : this) {
                    String molTitle = m.getTitle();
                    Map<String, Object> molInfo = SDFReader.parseString(molTitle);
                    int molid = this.db.getMoleculeId(this.targetId, (String) molInfo.get("molid"));
                    if (molid == -1) {
                        String smiles = SDFReader.getCanonicalSmiles(m);
                        molid = this.db.insertMolecule(this.targetId,  (String) molInfo.get("molid"),
                                (boolean) molInfo.get("active"), smiles);
                    }
                    int confid = this.db.getConformerId(molid, (Integer) molInfo.get("conformer"));
                    if (confid == -1){
                        String ctab = SDFReader.getCTABString(m);
                        confid = this.db.insertConformer(molid, (Integer) molInfo.get("conformer"),
                                calcGeometricDescriptor(m), ctab, molTitle);
                    }
                    List<Volume> volumeList = this.calculator.projectMultiField(m, fieldIDs);
                    for (int i = 0; i < NTypes; i++) {
                        String mapname = this.volmaps.get(fieldIDs[i]);
                        // String key = molecule.getTitle() + "_" + m;
                        Volume volume = volumeList.get(i);
                        if (mapname.contains("hydroele")) {
                            Volume pos = new Volume(volume); // Copy
                            // NEGATIVE PART and POSITIVE part apart
                            prepareVolume(volume, 0, 20, true, true, 1.0);
                            prepareVolume(pos, 0, 20, false, true, 1.0);
                            calMomentsAndStore_psql(confid, mapname, volume, order, db);
                            calMomentsAndStore_psql(confid, mapname + "_pos", pos, order, db);
                        } else if (mapname.contains("hydro")) {
                            prepareVolume(volume, 0, 20, true, true, 1.0);
                            calMomentsAndStore_psql(confid, mapname, volume, order, db);
                        } else {
                            // hbond maps
                            prepareVolume(volume, 0, 100, false, false, 10);
                            calMomentsAndStore_psql(confid, mapname, volume, order, db);
                        }
                    }
                    // been processed
                    map.put(m, "Done");
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Process sdf file and insert data directly to a PSQL database
     * 
     * @param maps     maps to generate
     * @param sdfPath  path to sdf file to process
     * @param order    zernike order
     * @param nThreads
     * @throws CDKException
     */
    public static void process_sdf_PSQL(HashMap<Integer, String> maps, String target, String sdfPath,
            int order, int nThreads, VectorDatabaseImpl db) throws CDKException {

        // Create target in DB if not exists
        int targetID = db.insertTarget(target);

        try (IteratingSDFReader reader = SDFReader.read(sdfPath)) {
            ExecutorService pool = Executors.newFixedThreadPool(nThreads);
            ConcurrentHashMap<IAtomContainer, String> map = new ConcurrentHashMap<>();
            RegisterTask task = new RegisterTask(map, maps, db, targetID, order);
            while (reader.hasNext()) {
                IAtomContainer molecule = reader.next();
                task.add(molecule);
                if (task.size() >= 10 || !reader.hasNext()){
                    pool.submit(task);
                    task = new RegisterTask(map, maps, db, targetID, order);
                }
            }
            // Close the SDF reader
            reader.close();

            pool.shutdown();
            try {
                if (!pool.awaitTermination(800, TimeUnit.MILLISECONDS)) {
                    pool.shutdownNow();
                } 
            } catch (InterruptedException e) {
                pool.shutdownNow();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Command(name = "process_sdf", description = "Generates projections of molecules in sdf file(s) with hydrophobic values precomputed to moments and fingerprints inserted into rocksdbs.")
    public void process_sdf_psql(@Parameters(description = "SDF input") String sdfFile,
            @Option(names = {
                    "--target" }, required = true, description = "Target name") String target,
            @Option(names = {
                    "--folder" }, description = "Flag if input is a folder, will read all SDF files inside if that's the case. Default: FALSE (treat single file)") boolean isFolder,
            @Option(names = {
                    "--threads" }, defaultValue = "1", description = "Number of threads to use") int nThreads,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order (0 to 20)") int order)
            throws IOException {

        HashMap<Integer, String> maps = new HashMap<Integer, String>();
        maps.put(2, "hydroele");
        maps.put(3, "hydrocav");
        maps.put(5, "hbond_Donors");
        maps.put(6, "hbond_Acceptors");
        // maps.put("hydrototal", 1);
        // maps.put("hydrovdw", 4);

        VectorDatabaseImpl db = new VectorDatabaseImpl();

        // Process hydroele and cav
        if (isFolder) {
            List<String> result;
            try (Stream<Path> walk = Files.walk(Paths.get(sdfFile))) {
                result = walk
                        .filter(p -> !Files.isDirectory(p)) // not a directory
                        .map(p -> p.toString()) // convert path to string
                        .filter(f -> f.endsWith("sdf")) // check end with
                        .collect(Collectors.toList()); // collect all matched to a List
            }
            if (result != null) {
                int counter = 0;
                int N = result.size();
                long time0 = System.currentTimeMillis();
                for (String sdf : result) {

                    long time = System.currentTimeMillis();
                    counter++;
                    try {
                        process_sdf_PSQL(maps, target, sdf, order, nThreads, db);
                    } catch (CDKException e) {
                        // TODO Auto-generated catch block
                        e.printStackTrace();
                    }
                    long completedIn = System.currentTimeMillis() - time;
                    System.out.println("DONE WITH FILE (" + counter + "/" + N + ")" + sdf);
                    System.out.println(DurationFormatUtils.formatDuration(completedIn, "HH:mm:ss:SS"));
                }
                long completedAll = System.currentTimeMillis() - time0;
                System.out.println("DONE PROCESSING SDF FOLDER " + sdfFile);
                System.out.println(DurationFormatUtils.formatDuration(completedAll, "HH:mm:ss:SS"));
            } else {
                System.out.println("NOTHING DONE ON " + sdfFile);
            }
        } else {
            try {
                long time0 = System.currentTimeMillis();
                process_sdf_PSQL(maps, target, sdfFile, order, nThreads, db);
                long completedAll = System.currentTimeMillis() - time0;
                System.out.println(DurationFormatUtils.formatDuration(completedAll, "HH:mm:ss:SS"));
            } catch (CDKException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            System.out.println("DONE PROCESSING SDF FILE " + sdfFile);
        }
        db.closeConnection();
    }

    @Command(name = "process_sdf_ext", description = "Generates projections of molecules in sdf file(s) with hydrophobic values precomputed to moments and fingerprints inserted into rocksdbs.")
    public void process_sdf(@Parameters(description = "SDF input") String sdfFile,
            @Option(names = {
                    "--folder" }, description = "Flag if input is a folder, will read all SDF files inside if that's the case. Default: FALSE (treat single file)") boolean isFolder,
            @Option(names = {
                    "--momentsdb" }, defaultValue = "moments.db", description = "RocksDB to store moments") String momentsDBPath,
            @Option(names = {
                    "--fpdb" }, defaultValue = "fp.db", description = "RocksDB to store invariant fingerprint") String fpDBPath,
            @Option(names = {
                    "--threads" }, defaultValue = "1", description = "Number of threads to use") int nThreads,
            @Option(names = {
                    "--out" }, defaultValue = "", description = "Prefix for text files containing the moments. Will print also the DX with same prefix grids.") String outPrefix,
            @Option(names = {
                    "--dx" }, description = "Print also the DX projection volumes? This is before the processing.") boolean printDX,
            @Option(names = {
                    "--writefiles" }, description = "Write to output files instead of rocksdb. Will use --out prefix") boolean writeFiles,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order (0 to 20)") int order)
            throws IOException {

        HashMap<Integer, String> maps = new HashMap<Integer, String>();
        maps.put(2, "hydroele");
        maps.put(3, "hydrocav");
        maps.put(5, "hbond_Donors");
        maps.put(6, "hbond_Acceptors");
        // maps.put("hydrototal", 1);
        // maps.put("hydrovdw", 4);

        RocksDBInterface momentsDB = null;
        RocksDBInterface FPDB = null;

        if (!writeFiles) {
            momentsDB = new RocksDBInterface(momentsDBPath);
            FPDB = new RocksDBInterface(fpDBPath);
        }

        // Process hydroele and cav
        if (isFolder) {
            List<String> result;
            try (Stream<Path> walk = Files.walk(Paths.get(sdfFile))) {
                result = walk
                        .filter(p -> !Files.isDirectory(p)) // not a directory
                        .map(p -> p.toString()) // convert path to string
                        .filter(f -> f.endsWith("sdf")) // check end with
                        .collect(Collectors.toList()); // collect all matched to a List
            }
            ;
            if (result != null) {
                int counter = 0;
                int N = result.size();
                for (String sdf : result) {
                    counter++;
                    process_sdf_(maps, sdf, momentsDB, FPDB, order, nThreads, writeFiles, printDX, outPrefix);
                    System.out.println("DONE WITH FILE (" + counter + "/" + N + ")" + sdf);
                }
                System.out.println("DONE PROCESSING SDF FOLDER " + sdfFile);
            } else {
                System.out.println("NOTHING DONE ON " + sdfFile);
            }
        } else {
            process_sdf_(maps, sdfFile, momentsDB, FPDB, order, nThreads, writeFiles, printDX, outPrefix);
            System.out.println("DONE PROCESSING SDF FILE " + sdfFile);
        }
        if (!writeFiles) {
            momentsDB.close();
            FPDB.close();
        }
    }

    @Command(name = "process_pharmscreen_POS", description = "Temporal routine to process again hydroele maps and only store the positive part.")
    public void process_pharmscreen_hydroelepos(@Parameters(description = "Folder to read DX files from") String folder,
            @Option(names = {
                    "--momentsdb" }, defaultValue = "moments.db", description = "RocksDB to store moments") String momentsDBPath,
            @Option(names = {
                    "--fpdb" }, defaultValue = "fp.db", description = "RocksDB to store invariant fingerprint") String fpDBPath,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order (0 to 20)") int order,
            @Option(names = {
                    "--rm" }, description = "Remove DX folder after processing") boolean remove)
            throws IOException {

        System.out.println("START PROCESSING PHARMSCREEN FOLDER " + folder);
        RocksDBInterface momentsDB = new RocksDBInterface(momentsDBPath);
        RocksDBInterface FPDB = new RocksDBInterface(fpDBPath);
        // Process hydroele and cav
        process_(folder, momentsDB, FPDB, order, 0.0, 20.0, false, true, false, "hydroele", 1.); // POSITIVE HYDROELE
        momentsDB.close();
        FPDB.close();
        if (remove) {
            // Remove dxf file after processing
            File f = new File(folder);
            FileUtils.deleteDirectory(f);
            System.out.println("PHARMSCREEN FOLDER DELETED: " + folder);
        }
        System.out.println("DONE PROCESSING PHARMSCREEN FOLDER " + folder);
    }

    @Command(name = "process_pharmscreen", description = "Compute 3D Zernike moments for all DX grids in folder generated by pharmscreen. Save moments to a rocksDB database.")
    public void process_pharmscreen(@Parameters(description = "Folder to read DX files from") String folder,
            @Option(names = {
                    "--momentsdb" }, defaultValue = "moments.db", description = "RocksDB to store moments") String momentsDBPath,
            @Option(names = {
                    "--fpdb" }, defaultValue = "fp.db", description = "RocksDB to store invariant fingerprint") String fpDBPath,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order (0 to 20)") int order,
            @Option(names = {
                    "--rm" }, description = "Remove DX folder after processing") boolean remove)
            throws IOException {

        System.out.println("START PROCESSING PHARMSCREEN FOLDER " + folder);
        RocksDBInterface momentsDB = new RocksDBInterface(momentsDBPath);
        RocksDBInterface FPDB = new RocksDBInterface(fpDBPath);
        // Process hydroele and cav
        process_(folder, momentsDB, FPDB, order, 0.0, 20.0, false, false, false, "hydro", 1.);
        process_(folder, momentsDB, FPDB, order, 0.0, 20.0, false, true, false, "hydroele", 1.); // POSITIVE HYDROELE
        process_(folder, momentsDB, FPDB, order, 0.0, 100.0, true, true, false, "hbond", 10.);
        momentsDB.close();
        FPDB.close();
        if (remove) {
            // Remove dxf file after processing
            File f = new File(folder);
            FileUtils.deleteDirectory(f);
            System.out.println("PHARMSCREEN FOLDER DELETED: " + folder);
        }
        System.out.println("DONE PROCESSING PHARMSCREEN FOLDER " + folder);
    }

    @Command(name = "process", description = "Compute 3D Zernike moments for all DX grids in current folder. Save moments to a rocksDB database.")
    public void process(@Parameters(description = "Folder to read DX files from") String folder,
            @Option(names = {
                    "--momentsdb" }, defaultValue = "moments.db", description = "RocksDB to store moments") String momentsDBPath,
            @Option(names = {
                    "--fpdb" }, defaultValue = "fp.db", description = "RocksDB to store invariant fingerprint") String fpDBPath,
            @Option(names = { "-N" }, required = true, description = "Zernike order (0 to 20)") int order,
            @Option(names = {
                    "--minCap" }, defaultValue = "0.0", description = "Volume preprocessing option: Discard values below this cutoff after flipping negatives to positives. Default 0.0.") double minCap,
            @Option(names = {
                    "--maxCap" }, defaultValue = "100.0", description = "Volume preprocessing option: Discard values above this cutoff after flipping negatives to positives. Default 100.0.") double maxCap,
            @Option(names = {
                    "--skipnormalize" }, description = "Volume preprocessing option: Do not normalize values after flipping negatives to positives. Default is to normalize.") boolean skipNormalize,
            @Option(names = {
                    "--skipflip" }, description = "Volume preprocessing option: Do not flip negative and positive values. Default is to flip.") boolean skipflip,
            @Option(names = {
                    "--rm" }, description = "Remove DX folder after processing") boolean remove,
            @Option(names = {
                    "--pattern" }, defaultValue = "", description = "Pattern to match DX files within folder. The filename should contain this substring") String pattern,
            @Option(names = {
                    "--multiplier" }, defaultValue = "1.0", description = "Volume preprocessing option: Multiplier to normalized volume. Used to scale the numbers. Default 1.") double multiplier)
            throws IOException {
        // find files matched `dx` file extension from folder
        List<String> result;
        try (Stream<Path> walk = Files.walk(Paths.get(folder))) {
            result = walk
                    .filter(p -> !Files.isDirectory(p)) // not a directory
                    .map(p -> p.toString()) // convert path to string
                    .filter(f -> f.endsWith("dx")) // check end with
                    .filter(f -> f.contains(pattern)) // meet pattern
                    .collect(Collectors.toList()); // collect all matched to a List
        }

        System.out.println("Processing folder: " + folder);
        System.out.println("Processing folder pattern: " + pattern);
        System.out.println("Moments DB: " + momentsDBPath);
        System.out.println("FP DB: " + fpDBPath);
        System.out.println("Processing folder: " + folder);
        System.out.println("ZERNIKE ORDER: " + order);
        System.out.println("Remove originals?: " + remove);
        if (result != null) {
            // If .dx files found, initialize rocksDB. If existing, it will just use the
            // existing one.
            RocksDBInterface momentsDB = new RocksDBInterface(momentsDBPath);
            RocksDBInterface FPDB = new RocksDBInterface(fpDBPath);
            // Process each DX file
            Volume v = new Volume();
            for (String dxf : result) {
                System.out.println("Processing volume file: " + dxf);
                try {
                    v = loadMoleculeVolume(dxf, minCap, maxCap, !skipflip, !skipNormalize, multiplier);
                } catch (Exception e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
                InvariantNorm n = new InvariantNorm(v, order);
                List<Double> fp = n.getFingerprint();
                List<Complex> moments = ZernikeMoments.flattenMomentsComplex(n.getMoments().getOriginalMoments());
                momentsDB.putMoments(dxf, moments);
                FPDB.putFP(dxf, fp);
            }
            if (remove) {
                // Remove dxf file after processing
                File f = new File(folder);
                FileUtils.deleteDirectory(f);
                System.out.println("Deleted folder: " + folder);
            }
            System.out.println("DONE PROCESSING FOLDER " + folder);
        }
    }

    public void process_(String folder, RocksDBInterface momentsDB, RocksDBInterface fpDB, int order, double minCap,
            double maxCap,
            boolean skipNormalize, boolean skipflip, boolean remove, String pattern, double multiplier)
            throws IOException {
        // find files matched `dx` file extension from folder
        List<String> result;
        try (Stream<Path> walk = Files.walk(Paths.get(folder))) {
            result = walk
                    .filter(p -> !Files.isDirectory(p)) // not a directory
                    .map(p -> p.toString()) // convert path to string
                    .filter(f -> f.endsWith("dx")) // check end with
                    .filter(f -> f.contains(pattern)) // meet pattern
                    .collect(Collectors.toList()); // collect all matched to a List
        }

        System.out.println("Processing folder: " + folder);
        System.out.println("Processing folder pattern: " + pattern);
        System.out.println("Processing folder: " + folder);
        System.out.println("ZERNIKE ORDER: " + order);
        System.out.println("Remove originals?: " + remove);
        if (result != null) {
            // If .dx files found, initialize rocksDB. If existing, it will just use the
            // existing one.
            // Process each DX file
            Volume v = new Volume();
            for (String dxf : result) {
                String key = "";
                if (dxf.contains("hydroele")) {
                    key = (skipflip) ? dxf + "_POS" : dxf;
                } else {
                    key = dxf;
                }
                System.out.println("Processing volume file: " + dxf);
                try {
                    v = loadMoleculeVolume(dxf, minCap, maxCap, !skipflip, !skipNormalize, multiplier);
                } catch (Exception e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
                InvariantNorm n = new InvariantNorm(v, order);
                List<Double> fp = n.getFingerprint();
                List<Complex> moments = ZernikeMoments.flattenMomentsComplex(n.getMoments().getOriginalMoments());
                momentsDB.putMoments(key, moments);
                fpDB.putFP(key, fp);
            }
            if (remove) {
                // Remove dxf file after processing
                File f = new File(folder);
                FileUtils.deleteDirectory(f);
                System.out.println("Deleted folder: " + folder);
            }
            System.out.println("DONE PROCESSING FOLDER " + folder);
        }
    }

    @Command(name = "moments", description = "Compute 3D Zernike moments")
    public void moments(
            @Option(names = { "-i" }, required = true, description = "Input Volume in DX format") String volumePath,
            @Option(names = {
                    "-o" }, required = true, description = "Destination file for full set of complex coefficients") String dest,
            @Option(names = {
                    "-f" }, description = "OPTIONAL Destination file for invariant coefficients") String fpdest,
            @Option(names = { "-N" }, required = true, description = "Zernike order (0 to 20)") int order,
            @Option(names = {
                    "--minCap" }, defaultValue = "0.0", description = "Volume preprocessing option: Discard values below this cutoff after flipping negatives to positives. Default 0.0.") double minCap,
            @Option(names = {
                    "--maxCap" }, defaultValue = "100.0", description = "Volume preprocessing option: Discard values above this cutoff after flipping negatives to positives. Default 100.0.") double maxCap,
            @Option(names = {
                    "--skipnormalize" }, description = "Volume preprocessing option: Do not normalize values after flipping negatives to positives. Default is to normalize.") boolean skipNormalize,
            @Option(names = {
                    "--skipflip" }, description = "Volume preprocessing option: Do not flip negative and positive values. Default is to flip.") boolean skipflip,
            @Option(names = {
                    "--multiplier" }, defaultValue = "1.0", description = "Volume preprocessing option: Multiplier to normalized volume. Used to scale the numbers. Default 1.") float multiplier)
            throws IOException {

        Volume v = new Volume();
        try {
            v = loadMoleculeVolume(volumePath, minCap, maxCap, !skipflip, !skipNormalize, multiplier);
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        InvariantNorm n = new InvariantNorm(v, order);
        n.getMoments().write(dest, false);
        System.out.println("Projecting volume file: " + volumePath + " and saving moments to: " + dest);

        if (fpdest != null) {
            List<Double> fp = n.getFingerprint();
            ZernikeMomentsIO.writeDouble(fpdest, fp);
        }
    }

    @Command(name = "reconstruct", description = "Reconstruct volume from 3D Zernike moments file")
    public void reconstruct(
            @Option(names = {
                    "-i" }, required = true, description = "Input text file with moments (full complex number set)") String coefPath,
            @Option(names = { "-o" }, required = true, description = "Destination DX file") String dest,
            @Option(names = { "-N" }, required = true, description = "Zernike order used") int order,
            @Option(names = {
                    "--dim" }, defaultValue = "32", description = "Reconstruction grid dimensions (n x n x n). Default is 32.") int dim,
            @Option(names = {
                    "--origin" }, description = "OPTIONAL Grid origin (x,y,z)", converter = DoubleListConverter.class) List<Double> origin,
            @Option(names = {
                    "--scale" }, defaultValue = "2.", description = "Reconstruction re-scaling factor") double scaling,
            @Option(names = {
                    "--width" }, defaultValue = "1.", description = "Reconstruction grid Width (cell size)") double width,
            @Option(names = {
                    "--useCache" }, description = "Reconstruction use saved cache info to speed up calculations. DEFAULT FALSE.") boolean useCache,
            @Option(names = {
                    "--saveCache" }, description = "Reconstruction SAVE cache info to speed up future calculations. DEFAULT FALSE.") boolean saveCache)
            throws IOException {

        ZernikeMoments m = null;
        try {
            m = ZernikeMoments.read(coefPath, order, false);
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        Volume reconstructVolume = new Volume();
        try {
            reconstructVolume = ZernikeMoments.reconstructVolume(m, dim, order, useCache, saveCache);
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        if (origin != null) {
            double[] origin_d = { origin.get(0), origin.get(1), origin.get(2) };
            reconstructVolume.setCorner(origin_d);
        }
        OpenDXIO.write(dest, reconstructVolume);
        System.out.println("Volume reconstructed " + dest);
    }

    @Command(name = "fp", description = "Print Invariant Fingerprint from full set of coefficients")
    public void fp(
            @Option(names = {
                    "-i" }, required = true, description = "Input text file with moments (full complex number set)") String coefPath,
            @Option(names = { "-N" }, required = true, description = "Zernike order used") int order,
            @Option(names = { "-o" }, description = "Optional output file. Else print to STDOUT.") String dest)
            throws IOException {

        ZernikeMoments m = null;
        try {
            m = ZernikeMoments.read(coefPath, order, false);
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        double[] center = { 0.0, 0.0, 0.0 };
        InvariantNorm n = new InvariantNorm(m.getOriginalMoments(), true, center);
        List<Double> fp = n.getFingerprint();

        if (dest != null) {
            ZernikeMomentsIO.writeDouble(dest, fp);
        } else {
            String fp_str = fp.stream().map(Object::toString).collect(Collectors.joining(",")).toString();
            // String fp_str = "";
            // for (int i=0; i < fp.size(); i++){
            // fp_str += fp.get(i);
            // if (i < fp.size()-1) fp_str += ",";
            // }
            System.out.println(fp_str);
        }
    }

    @Command(name = "distance", description = "Compute distance between 2 invariant fp files")
    public void distance(
            @Option(names = {
                    "-i1" }, required = true, description = "FIRST Input text file with invariant moments (fp)") String input1Path,
            @Option(names = {
                    "-i2" }, required = true, description = "SECOND Input text file with invariant moments (fp)") String input2Path)
            throws IOException, ClassNotFoundException {

        List<Double> invariantsThis = ZernikeMomentsIO.readDouble(input1Path);
        List<Double> invariantsOther = ZernikeMomentsIO.readDouble(input2Path);

        double sumDiffs = 0;
        for (int k = 0; k < invariantsThis.size(); k++) {
            double diffAbs = Math.abs(invariantsThis.get(k) - invariantsOther.get(k));
            double sumAbs = invariantsThis.get(k) + invariantsOther.get(k) + 1;
            sumDiffs += diffAbs / sumAbs;
        }
        System.out.println(sumDiffs);
    }

    @Command(name = "superpose", description = "Superpose two volumes and rotate a PDB file corresponding to second volume")
    public void superpose(
            @Option(names = {
                    "-v1" }, required = true, description = "Input Volume DX file for target") String input1Path,
            @Option(names = {
                    "-v2" }, required = true, description = "Input Volume DX file for moving part to superposte to target") String input2Path,
            @Option(names = {
                    "-m2" }, required = true, description = "Molecule to move according to volume superposition in PDB format") String pdbFile,
            @Option(names = {
                    "-o" }, required = true, description = "Output destination file. Second molecule fitted to first volume (pdb file) ") String dest,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order to use ") int order,
            @Option(names = {
                    "--minCap" }, defaultValue = "0.0", description = "Volume preprocessing option: Discard values below this cutoff after flipping negatives to positives. Default 0.0.") double minCap,
            @Option(names = {
                    "--maxCap" }, defaultValue = "100.0", description = "Volume preprocessing option: Discard values above this cutoff after flipping negatives to positives. Default 100.0.") double maxCap,
            @Option(names = {
                    "--skipnormalize" }, description = "Volume preprocessing option: Do not normalize values after flipping negatives to positives. Default is to normalize.") boolean skipNormalize)
            throws IOException, ClassNotFoundException {

        Volume v1 = new Volume();
        Volume v2 = new Volume();
        Structure mol = null;
        try {
            v1 = loadMoleculeVolume(input1Path, minCap, maxCap, true, !skipNormalize, 1);
            v2 = loadMoleculeVolume(input2Path, minCap, maxCap, true, !skipNormalize, 1);
            mol = loadPDB(pdbFile);
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        InvariantNorm n1 = new InvariantNorm(v1, order);
        InvariantNorm n2 = new InvariantNorm(v2, order);
        AlignmentResult alignment = n2.alignTo(n1);
        Matrix4d R = alignment.getTransforms().get(0);
        Calc.transform(mol, R);
        writeStructurePDB(mol, dest);
        System.out.println(getMatrixString(R, "transform"));
        System.out.println("SCORE: " + alignment.getScore());
    }

    @Command(name = "rotate", description = "Superpose two moments and output rotation matrix to STDOUT")
    public void rotate(
            @Option(names = {
                    "-m1" }, required = true, description = "Input moments file for target (static)") String input1Path,
            @Option(names = {
                    "-m2" }, required = true, description = "Input moments file for moving part to superposte to target") String input2Path,
            @Option(names = {
                    "-N" }, required = true, description = "Zernike order to use ") int order)
            throws IOException, ClassNotFoundException {

        ZernikeMoments m1 = null;
        ZernikeMoments m2 = null;
        try {
            m1 = ZernikeMoments.read(input1Path, order, false);
            m2 = ZernikeMoments.read(input2Path, order, false);
        } catch (ClassNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        double[] center = { 0, 0, 0 };
        InvariantNorm n1 = new InvariantNorm(m1, center);
        InvariantNorm n2 = new InvariantNorm(m2, center);
        AlignmentResult alignment = n2.alignTo(n1);
        Matrix4d R = alignment.getTransforms().get(0);
        // Calc.transform(mol, R);
        // writeStructurePDB(mol, dest);
        System.out.println(getMatrixString(R, "transform"));
        System.out.println("SCORE: " + alignment.getScore());
    }

    public static void main(String[] args) {
        CommandLine commandLine = new CommandLine(new LigZernike());
        int exitCode = commandLine.execute(args);
        System.exit(exitCode);

    }

}
