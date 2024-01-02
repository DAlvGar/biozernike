package org.rcsb.biozernike.ligzernike;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.apache.commons.lang3.StringUtils;
import org.rcsb.biozernike.AlignmentResult;
import org.rcsb.biozernike.InvariantNorm;
import org.rcsb.biozernike.RotationAlignment;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Hello world!
 *
 */
@Command(name = "LigZernike", mixinStandardHelpOptions = true, description = "LigZernike will extract fingerprints for the two volumes given, calcualte a distance and return the second moelcule superposed to the first", version = "1.0")
public class LigZernike {
    @Option(names = { "-v1", "--vol1" }, required = true, description = "Volume for molecule 1")
    private static String volumepath1;
    @Option(names = { "-v2", "--vol2" }, required = true, description = "Volume for molecule 2")
    private static String volumepath2;
    @Option(names = { "-m", "--mol" }, required = false, description = "PDB file for molecule 2 if you want it to be superposed to 1")
    private static String molPDB2;
    @Option(names = { "-o", "--out" }, required = false, description = "Output aligned pdb file with second molecule")
    private static String outPDB;

    public static void main(String[] args) throws FileNotFoundException {
        int exitCode = new CommandLine(new LigZernike()).execute(args);
        if (StringUtils.isNotEmpty(volumepath1) && StringUtils.isNotEmpty(volumepath2)) {
            Volume v1 = null;
            Volume v2 = null;
            try {
                v1 = OpenDXIO.read(volumepath1);
                v2 = OpenDXIO.read(volumepath2);
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }

            // Compute invariant norms from volumes
            v1.flipValues(); // Turn negative to positive and viceversa
            v1.capMax(10); // throw away positive values, cap to zero and keep negatives
            v1.applyContourAndNormalize(0, 1); // Eliminate values below zero and normalize
            v1.updateCenter();
            v2.flipValues();
            v2.capMax(10); // throw away positive values, cap to zero and keep negatives
            v2.applyContourAndNormalize(0, 1);
            v2.updateCenter(); //keep the same center
            
            System.out.println("Center v1 :"+v1.getCenterReal()[0]+" "+v1.getCenterReal()[1]+" "+v1.getCenterReal()[2]);
            System.out.println("Center v2 :"+v2.getCenterReal()[0]+" "+v2.getCenterReal()[1]+" "+v2.getCenterReal()[2]);
            
            InvariantNorm m1 = new InvariantNorm(v1, 10);
            InvariantNorm m2 = new InvariantNorm(v2, 10);

            // Print 3DZD fingerprint for each molecule
            /*
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

            // Calculate distance
            for (int i = 0; i < 3; i++) {
                double d = m1.compareInvariants(m2, i);
                System.out.println(String.format("%.2f",d));
            }

            if (molPDB2 != null) {

                // Superpose
                List<InvariantNorm> both = new ArrayList<>();
                both.add(m2);
                both.add(m1);
                AlignmentResult alignmentResult = RotationAlignment.alignMultiple(both);

                // alignment quality is reasonable
                System.out.println("Alignment score: " + alignmentResult.getScore());

                // IF A MOL2 PDB WAS GIVEN, ALIGN AND OUTPUT THE TRANSFORMED MOL
                if (StringUtils.isNotEmpty(molPDB2)) {
                    Matrix4d rot1 = alignmentResult.getTransforms().get(0);
                    Matrix4d rot2 = alignmentResult.getTransforms().get(1);
                    rot1.invert();
                    rot1.mul(rot2);
                    PDBFileReader pdbreader = new PDBFileReader();

                    try {
                        Structure structure = pdbreader.getStructure(molPDB2);
                        Calc.transform(structure, rot1);
                        try (PrintWriter out = new PrintWriter(outPDB)) {
                            out.println(structure.toPDB());
                        }
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }
        }

    }

}


