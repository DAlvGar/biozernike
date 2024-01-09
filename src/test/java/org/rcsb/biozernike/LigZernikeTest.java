package org.rcsb.biozernike;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Test;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.volume.VolumeIO;
import org.rcsb.biozernike.zernike.ZernikeMoments;
import org.rcsb.biozernike.zernike.ZernikeMomentsIO;

public class LigZernikeTest {

    public Volume loadMoleculeVolume(String filename) throws Exception {
        return loadMoleculeVolume(filename, 0.2, 20., true, 1.);
    }

    public Volume loadMoleculeVolume(String filename, double minCap, double maxCap, boolean flip, double multiplier)
            throws Exception {
        String path = LigZernikeTest.class.getResource("/" + filename).getPath();
        Volume v1 = OpenDXIO.read(path);
        if (flip)
            v1.flipValues(); // Turn negative to positive and viceversa
        v1.capMax(maxCap); // throw away positive values, cap to zero and keep negatives
        // v1.capMin(0);
        v1.applyContourAndNormalize(minCap, multiplier); // Eliminate values below 0.2 and normalize
        v1.updateCenter();
        return v1;
    }

    // @Test
    public void testMomentsWrite() throws Exception {
        String outPrefix = "test_norms";
        Volume vol1 = loadMoleculeVolume("mol1_HDON.dx");
        InvariantNorm n = new InvariantNorm(vol1, 15);
        ZernikeMomentsIO.writeDouble(outPrefix + "_fp.txt", n.getFingerprint());
        n.getMoments().write(outPrefix + "_complex.txt", false);
    }

    // @Test
    public void testMomentsRead() throws Exception {
        /*
         * object 1 class gridpositions counts 45 49 45
         * origin -23.410000 -29.040000 -18.110000
         * delta 1 0 0
         * delta 0 1 0
         * delta 0 0 1
         */
        Volume vol1 = loadMoleculeVolume("mol1_HDON.dx");
        System.out.println(vol1.getCenterVolume()[0]+" "+vol1.getCenterVolume()[1]+" "+vol1.getCenterVolume()[2]);
        int[] dimensions = { 45, 49, 45 };
        double[] origin = { -23.41, -29.04, -18.11 };
        double[] center = { dimensions[0] / 2., dimensions[1] / 2., dimensions[2] / 2. };
        //double[] center = {0,0,0};
        int _maxN = 15;
        double gridWidth = 1;
        double scale_factor = 2.;
        ZernikeMoments m = ZernikeMoments.read("test_norms_complex.txt", 15, false);
        Volume reconstructVolume = ZernikeMoments.reconstructVolume(m, dimensions, center, _maxN, gridWidth,
                scale_factor, false, false);
        reconstructVolume.setCorner(origin);
        OpenDXIO.write("reconstruct_mol1_HDON_fromfile.dx", reconstructVolume);
    }

    // @Test
    public void testMIP() throws Exception {
        int maxorder = 15;
        double minCap = 0.5;
        double maxCap = 20.;
        boolean flip = true;
        double multiplier = 1;
        double scaleFactor = 3.; // Reconstruction scaling factor
        String reconstructionFolder = "reconstruction_05_20_100";
        File directory = new File("Ligands/" + reconstructionFolder);
        if (!directory.exists()) {
            directory.mkdir();
            // If you require it to make the entire directory path including parents,
            // use directory.mkdirs(); here instead.
        }
        String[] names = { "mol1", "mol2", "mol3", "mol4", "mol5" };
        String probe = "NEG";
        String matrixFile = "Ligands/transformations.dat";
        String pymolFile = "Ligands/pymol.pml";
        int N = names.length;

        BufferedWriter mwriter = new BufferedWriter(new FileWriter(matrixFile));
        BufferedWriter pwriter = new BufferedWriter(new FileWriter(pymolFile));

        List<String> basenames = new ArrayList<>();
        List<Volume> volumeArray = new ArrayList<>();
        List<Structure> structureArray = new ArrayList<>();
        List<InvariantNorm> invariantList = new ArrayList<>();

        // Populate arrays

        for (int i = 0; i < N; i++) {
            basenames.add(names[i] + "_" + probe);
            // Volume v = loadMoleculeVolume(basenames.get(i) + ".dx");
            Volume v = loadMoleculeVolume(basenames.get(i) + ".dx", minCap, maxCap, flip, multiplier);
            volumeArray.add(v);
            Structure s = loadPDB(names[i] + ".pdb");
            structureArray.add(s);
            invariantList.add(new InvariantNorm(volumeArray.get(i), maxorder));
        }

        // Compute distances against first element
        System.out.println("COMPUTE INVARIANT DISTANCE ALL AGAINST FIRST");
        for (int i = 1; i < N; i++) {
            double d = invariantList.get(0).compareInvariants(invariantList.get(i), 0); // 0 is basic fingerprints
            System.out.println("Distance 0-" + i + ": " + d);
        }

        System.out.println("ALIGN INDEPENDENTLY ALL AGAINST FIRST");
        List<Matrix4d> matrixarray = new ArrayList<>();
        // Compute alignment scores against first
        for (int i = 1; i < N; i++) {
            AlignmentResult alignmentResult = align(invariantList.get(0), invariantList.get(i), maxorder);
            System.out.println("INDEPENDENT ALIGMENT SCORE 0-" + i + ": " + alignmentResult.getScore());
            Matrix4d rot1 = alignmentResult.getTransforms().get(0);
            Matrix4d rot2 = alignmentResult.getTransforms().get(1);
            rot1.invert();
            rot1.mul(rot2);
            // printMatrix(rot2, "indep_rot" + i);
            matrixarray.add(rot1);
            String ms = getMatrixString(rot1, "indep_rot" + i);
            mwriter.write(ms);
            pwriter.write(ms);
        }
        mwriter.close();

        // System.out.println("ALIGN ALL AT ONCE");
        // AlignmentResult alignmentResult =
        // RotationAlignment.alignMultiple(invariantList, maxorder, false);
        // List<Matrix4d> matrixarray = alignmentResult.getTransforms();
        // for (int i = 0; i < matrixarray.size(); i++) {
        // printMatrix(matrixarray.get(i), "all_rot" + i);
        // String ms = getMatrixString(matrixarray.get(i), "all_rot" + i);
        // mwriter.write(ms);
        // pwriter.write(ms);
        // }

        System.out.println("RECONSTRUCT AND WRITE ALL VOLUMES");
        for (int i = 0; i < N; i++) {
            System.out.println(" - Volume " + i);
            Volume volume_rec2 = ZernikeMoments.reconstVolumeLike(invariantList.get(i).getMoments(), volumeArray.get(i),
                    maxorder, scaleFactor, true, false);
            String volname = basenames.get(i) + "_" + maxorder;
            String volfile = reconstructionFolder + "/" + volname + ".dx";
            String originalVol = basenames.get(i) + ".dx";
            OpenDXIO.write("Ligands/" + volfile, volume_rec2);
            pwriter.write("load " + volfile + "," + volname + "\n");
            pwriter.write("load originals/" + originalVol + "\n");
            if (i > 0) {
                pwriter.write("cmd.transform_object(\"" + volname + "\", indep_rot" + i + ")\n");
                pwriter.write("cmd.transform_object(\"" + basenames.get(i) + "\", indep_rot" + i + ")\n");
            }
            pwriter.write("isomesh " + volname + "_M," + volname + ", 50" + "\n");
            pwriter.write("isomesh " + basenames.get(i) + "_M," + basenames.get(i) + ", -1" + "\n");
            // pwriter.write("isosurface "+volname+"_S,"+volname+", 2"+"\n");
        }

        System.out.println("TRANSFORM ALL AND WRITE (ALL AGAINST ALL)");
        for (int i = 0; i < N; i++) {
            String pdbname = "";
            if (i == 0) {
                pdbname = basenames.get(i) + "_original";
            } else {
                Calc.transform(structureArray.get(i), matrixarray.get(i - 1));
                pdbname = basenames.get(i) + "_fitted";
            }
            String pdbfile = pdbname + ".pdb";
            writeStructurePDB(structureArray.get(i), "Ligands/" + pdbfile);
            pwriter.write("load " + pdbfile + "," + pdbname + "\n");
        }

        pwriter.close();
        System.out.println("DONE");
    }

    // @Test
    public void testPhm() throws Exception {
        Volume vol1 = loadMoleculeVolume("phmscreen_ref_hydroele.dx");
        Volume vol2 = loadMoleculeVolume("phmscreen_good_hydroele.dx");
        Volume vol3 = loadMoleculeVolume("phmscreen_bad_hydroele.dx");
        Structure mol1 = loadPDB("phmscreen_ref.pdb");
        Structure mol2 = loadPDB("phmscreen_good.pdb");
        Structure mol3 = loadPDB("phmscreen_bad.pdb");
        /*
         * vol1.setRadiusVarMult(2);
         * vol2.setRadiusVarMult(2);
         * vol3.setRadiusVarMult(2);
         */

        /*
         * VolumeIO.write(vol1,"Ligands/mol1_original.mrc",MapFileType.MRC);
         * VolumeIO.write(vol2,"Ligands/mol2_original.mrc",MapFileType.MRC);
         * VolumeIO.write(vol3,"Ligands/mol3_original.mrc",MapFileType.MRC);
         */

        List<InvariantNorm> zc = new ArrayList<>();

        InvariantNorm normalization1 = new InvariantNorm(vol1, 15);
        InvariantNorm normalization2 = new InvariantNorm(vol2, 15);
        InvariantNorm normalization3 = new InvariantNorm(vol3, 15);

        double d = 100;
        d = normalization1.compareInvariants(normalization2, 0); // 0 is basic fingerprints
        System.out.println("Distance ref to good: " + d);

        d = normalization1.compareInvariants(normalization3, 0); // 0 is basic fingerprints
        System.out.println("Distance ref to bad: " + d);

        zc.add(normalization1);
        zc.add(normalization2);
        zc.add(normalization3);
        AlignmentResult alignmentResult = RotationAlignment.alignMultiple(zc, 15, false);

        System.out.println("ALIGNMENT SCORE (global): " + alignmentResult.getScore());

        Matrix4d rot1 = alignmentResult.getTransforms().get(0);
        Matrix4d rot2 = alignmentResult.getTransforms().get(1);
        Matrix4d rot3 = alignmentResult.getTransforms().get(2);
        rot1.invert(); // this is the way of mol1 to the center we need to undo this way by inverting
                       // the matrix
        Matrix4d rot21 = (Matrix4d) rot1.clone();
        rot21.mul(rot2); // concatenate transformation of molecule 2 to the center followed by the way of
                         // mol1 to the original position
        Matrix4d rot31 = (Matrix4d) rot1.clone();
        rot31.mul(rot3);

        // Calc.transform(mol1,rot1);
        Calc.transform(mol2, rot21); // superpose mol2 onto mol1
        Calc.transform(mol3, rot31); // superpose mol3 onto mol1

        // Create a new Structure object for the ligand
        writeStructurePDB(mol1, "Ligands/phmscreen_ref.original.pdb");
        writeStructurePDB(mol2, "Ligands/phmscreen_good.fitted.pdb");
        writeStructurePDB(mol3, "Ligands/phmscreen_bad.fitted.pdb");
    }

    // @Test
    public void testReconstruction() throws Exception {
        int maxOrder = 15;
        String mol = "mol5";
        String probe = "HDON";
        String basename = mol + "_" + probe;
        Volume volume1 = loadMoleculeVolume(basename + ".dx");
        ZernikeMoments zm = new ZernikeMoments(volume1, maxOrder);
        // describeVolume(volume1);
        // Volume volume_rec = ZernikeMoments.reconstructVolume(zm, 32, 15, false,
        // true);
        // VolumeIO.write(volume_rec, "D:/PT/Ligands/volume1_rec_orig.map",
        // MapFileType.MRC);

        for (int maxN = 5; maxN <= maxOrder; maxN += 5) {
            int[] nTerms = CalcNumZernikeTerms.calcNumTerms(maxN);
            System.out.println("ORDER " + maxN + " Number of coefficients: " + nTerms[0] + " total " + nTerms[1]
                    + " positive " + nTerms[2] + " invariant");
            Volume volume_rec2 = ZernikeMoments.reconstVolumeLike(zm, volume1, maxN, 1.6, true, true);
            VolumeIO.write(volume_rec2, "Ligands/orders/" + basename + "_" + maxN + ".mrc", MapFileType.MRC);
        }
    }

    @Test
    public void testLigands() throws Exception {
        int max_order = 15;
        String[] fileNames = { "mol1_HDON.dx", "mol2_HDON.dx", "mol3_HDON.dx", "mol4_HDON.dx",
                "phmscreen_ref_hydroele.dx" };
        int n_files = fileNames.length;
        Volume[] volumeArray = new Volume[n_files];
        List<InvariantNorm> invariantsArray = new ArrayList<>();
        // InvariantNorm[] invariantsArray = new InvariantNorm[n_files];

        // Load and process volumes
        for (var i = 0; i < n_files; i++) {
            volumeArray[i] = loadMoleculeVolume(fileNames[i]);
        }

        // Compute invariant norms from volumes
        for (var i = 0; i < n_files; i++) {
            invariantsArray.add(new InvariantNorm(volumeArray[i], max_order));
        }

        // Calculate distances
        for (var i = 0; i < n_files; i++) {
            InvariantNorm invariant1 = invariantsArray.get(i);
            for (var j = i + 1; j < n_files; j++) {
                double d = invariant1.compareInvariants(invariantsArray.get(j), 0); // 0 is basic fingerprints
                System.out.println("Distance " + i + " - " + j + " : " + d);
            }
        }

    }

    public Structure loadPDB(String filename) throws Exception {
        String path = LigZernikeTest.class.getResource("/" + filename).getPath();
        PDBFileReader pdbreader = new PDBFileReader();
        Structure structure = pdbreader.getStructure(path);
        return structure;
    }

    public void superpose(InvariantNorm m1, InvariantNorm m2) {
        // Superpose
        List<InvariantNorm> both = new ArrayList<>();
        both.add(m1);
        both.add(m2);
        AlignmentResult alignmentResult = RotationAlignment.alignMultiple(both);

        // alignment quality is reasonable
        System.out.println("Alignment score: " + alignmentResult.getScore());

        // IF A MOL2 PDB WAS GIVEN, ALIGN AND OUTPUT THE TRANSFORMED MOL
        Matrix4d rot1 = alignmentResult.getTransforms().get(0);
        Matrix4d rot2 = alignmentResult.getTransforms().get(1);
        rot1.invert();
        rot1.mul(rot2);
        PDBFileReader pdbreader = new PDBFileReader();

        try {
            URL url3 = OpenDXIOTest.class.getResource("/mol3.pdb");
            Structure structure = pdbreader.getStructure(url3.getPath());

            Calc.transform(structure, rot1);
            try (PrintWriter out = new PrintWriter("mol3_fitted.pdb")) {
                out.println(structure.toPDB());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void writeStructurePDB(Structure s, String fileName) throws FileNotFoundException {
        try (PrintWriter out = new PrintWriter(fileName)) {
            out.println(s.toPDB());
        }
    }

    public void describeVolume(Volume vol) {
        System.out.println("Volume Dimensions: " + vol.getDimensions()[0] + " " + vol.getDimensions()[1] + " "
                + vol.getDimensions()[2]);
        System.out.println("Volume Center: " + vol.getCenterVolume()[0] + " " + vol.getCenterVolume()[1] + " "
                + vol.getCenterVolume()[2]);
        System.out.println("Volume Real Center: " + vol.getCenterReal()[0] + " " + vol.getCenterReal()[1] + " "
                + vol.getCenterReal()[2]);
        System.out
                .println("Volume Origin: " + vol.getCorner()[0] + " " + vol.getCorner()[1] + " " + vol.getCorner()[2]);
        System.out.println("Volume spacing: " + vol.getGridWidth());
        System.out.println("Volume mass: " + vol.getVolumeMass());
        System.out.println("Volume original mass: " + vol.getOriginalVolumeMass());
        System.out.println("Volume getRadiusVarMult: " + vol.getRadiusVarMult());
        System.out.println("Volume stats: " + vol.getDescriptiveStatistics());

    }

    public AlignmentResult align(InvariantNorm m1, InvariantNorm m2, int maxOrder) {
        List<InvariantNorm> zc = new ArrayList<>();
        zc.add(m1);
        zc.add(m2);
        return RotationAlignment.alignMultiple(zc, maxOrder, false);
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

    public void printMatrix(Matrix4d m, String name) {
        System.out.println(getMatrixString(m, name));
    }

}
