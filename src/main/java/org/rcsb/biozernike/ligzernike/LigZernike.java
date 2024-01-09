package org.rcsb.biozernike.ligzernike;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import javax.vecmath.Matrix4d;

import org.apache.commons.lang3.StringUtils;
import org.rcsb.biozernike.AlignmentResult;
import org.rcsb.biozernike.InvariantNorm;
import org.rcsb.biozernike.RotationAlignment;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.zernike.ZernikeMoments;
import org.rcsb.biozernike.zernike.ZernikeMomentsIO;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.io.PDBFileReader;

import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

/**
 * Hello world!
 *
 */
@Command(name = "LigZernike", mixinStandardHelpOptions = true, description = "Tool to work with volumetric files and 3D zernikes for ligands", version = "1.0")
public class LigZernike {

    public Volume loadMoleculeVolume(String filename, double minCap, double maxCap, boolean flip, boolean normalize,
            double multiplier)
            throws Exception {
        Volume v1 = OpenDXIO.read(filename);
        if (flip)
            v1.flipValues(); // Turn negative to positive and viceversa
        v1.capMax(maxCap); // throw away positive values, cap to zero and keep negatives
        v1.capMin(minCap);
        if (normalize)
            v1.applyContourAndNormalize(minCap, multiplier); // Eliminate values below 0.2 and normalize
        v1.updateCenter();
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

    public static void main(String[] args) {
        CommandLine commandLine = new CommandLine(new LigZernike());
        int exitCode = commandLine.execute(args);
        System.exit(exitCode);

    }

}
