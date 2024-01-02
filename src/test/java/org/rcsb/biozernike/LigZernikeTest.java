package org.rcsb.biozernike;

import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.apache.commons.lang.StringUtils;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.junit.Test;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;

public class LigZernikeTest {
    
    public Volume loadMolecule(String filename) throws Exception {
        String path = LigZernikeTest.class.getResource("/" + filename).getPath();
        System.err.println(path);
        Volume v1 = OpenDXIO.read(path);
        v1.flipValues(); // Turn negative to positive and viceversa
        //v1.capMax(20); // throw away positive values, cap to zero and keep negatives
        //v1.capMin(0);
        v1.applyContourAndNormalize(0, 10); // Eliminate values below zero and normalize
        v1.updateCenter();
        System.out.println(String.format("V1 center: %.3f %.3f %.3f", v1.getCenterReal()[0], v1.getCenterReal()[1],
                v1.getCenterReal()[2]));

        /*
         * BufferedWriter writer = new BufferedWriter(new FileWriter("v1_center.pdb"));
         * writer.write(String.
         * format("ATOM      1  N   CEN A   1    %8.3f%8.3f%8.3f  1.00139.11           N\n"
         * ,
         * v1.getCenterReal()[0], v1.getCenterReal()[1],
         * v1.getCenterReal()[2]));
         * writer.close();
         * OpenDXIO.write("v1_prepared.dx", v1);
         */
        return v1;
    }

    public void superpose(InvariantNorm m1, InvariantNorm m2) {
        // Superpose
        List<InvariantNorm> both = new ArrayList<>();
        both.add(m2);
        both.add(m1);
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

    @Test
    public void testLigands() throws Exception {
        int max_order = 10;
        String[] fileNames = { "mol1_HYD.dx", "mol2_HYD.dx", "mol3_HYD.dx", "mol4_HYD.dx" };
        int n_files = fileNames.length;
        Volume[] volumeArray = new Volume[n_files];
        InvariantNorm[] invariantsArray = new InvariantNorm[n_files];

        // Load and process volumes
        for (var i = 0; i < n_files; i++) {
            volumeArray[i] = loadMolecule(fileNames[i]);
        }

        // Compute invariant norms from volumes
        for (var i = 0; i < n_files; i++) {
            invariantsArray[i] = new InvariantNorm(volumeArray[i], max_order);
        }

        // Calculate distances
        for (var i = 0; i < n_files; i++) {
            InvariantNorm invariant1 = invariantsArray[i];
            for (var j = i + 1; j < n_files; j++) {
                double d = invariant1.compareInvariants(invariantsArray[j], 2); // 0 is basic fingerprints
                System.out.println("Distance " + i + " - " + j + " : " + d);
            }
        }
        /*
         * // Print 3DZD fingerprint for each molecule
         * List<Double> fp1 = m1.getFingerprint();
         * List<Double> fp2 = m2.getFingerprint();
         * 
         * try (PrintWriter out = new PrintWriter("mol1_fp.txt")) {
         * for (Double d : fp1) {
         * out.println(d);
         * }
         * }
         * 
         * try (PrintWriter out = new PrintWriter("mol2_fp.txt")) {
         * for (Double d : fp2) {
         * out.println(d);
         * }
         * }
         */


    }
}
