package org.rcsb.biozernike;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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

public class LigZernikeTest {
    
    public Volume loadMoleculeVolume(String filename) throws Exception {
        String path = LigZernikeTest.class.getResource("/" + filename).getPath();
        Volume v1 = OpenDXIO.read(path);
        v1.flipValues(); // Turn negative to positive and viceversa
        v1.capMax(20); // throw away positive values, cap to zero and keep negatives
        //v1.capMin(0);
        v1.applyContourAndNormalize(0.2, 1); // Eliminate values below 0.2 and normalize
        v1.updateCenter();
        //System.out.println(String.format("V1 center: %.3f %.3f %.3f", v1.getCenterReal()[0], v1.getCenterReal()[1],
        //        v1.getCenterReal()[2]));

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

    public void writeStructurePDB(Structure s, String fileName) throws FileNotFoundException{
		try (PrintWriter out = new PrintWriter(fileName)) {
			out.println(s.toPDB());
		}
    }

    public void describeVolume(Volume vol){
        System.out.println("Volume Dimensions: "+vol.getDimensions()[0]+" "+vol.getDimensions()[1]+" "+vol.getDimensions()[2] );
        System.out.println("Volume Center: "+vol.getCenterVolume()[0]+" "+vol.getCenterVolume()[1]+" "+vol.getCenterVolume()[2] );
        System.out.println("Volume Real Center: "+vol.getCenterReal()[0]+" "+vol.getCenterReal()[1]+" "+vol.getCenterReal()[2] );
        System.out.println("Volume Origin: "+vol.getCorner()[0]+" "+vol.getCorner()[1]+" "+vol.getCorner()[2] );
        System.out.println("Volume spacing: "+vol.getGridWidth() );
        System.out.println("Volume mass: "+vol.getVolumeMass() );
        System.out.println("Volume original mass: "+vol.getOriginalVolumeMass() );
        System.out.println("Volume getRadiusVarMult: "+vol.getRadiusVarMult());
        System.out.println("Volume stats: "+vol.getDescriptiveStatistics());
        
    }

    public AlignmentResult align(InvariantNorm m1, InvariantNorm m2, int maxOrder){
        List<InvariantNorm> zc = new ArrayList<>();
        zc.add(m1);
		zc.add(m2);
        return RotationAlignment.alignMultiple(zc, maxOrder, false);
    }

	@Test
	public void testMIP() throws Exception {
        int maxorder = 10;
        Volume vol1 = loadMoleculeVolume("mol1_HDON.dx");
        Volume vol2 = loadMoleculeVolume("mol2_HDON.dx");
        Volume vol3 = loadMoleculeVolume("mol3_HDON.dx");
        Volume vol4 = loadMoleculeVolume("mol4_HDON.dx");
        Volume vol5 = loadMoleculeVolume("mol5_HDON.dx");
        Structure mol1 = loadPDB("mol1.pdb");
        Structure mol2 = loadPDB("mol2.pdb");
        Structure mol3 = loadPDB("mol3.pdb");
        Structure mol4 = loadPDB("mol4.pdb");
        Structure mol5 = loadPDB("mol5.pdb");
		/* vol1.setRadiusVarMult(2);
		vol2.setRadiusVarMult(2);
		vol3.setRadiusVarMult(2); */
		
		/* VolumeIO.write(vol1,"Ligands/mol1_original.mrc",MapFileType.MRC);
		VolumeIO.write(vol2,"Ligands/mol2_original.mrc",MapFileType.MRC);
		VolumeIO.write(vol3,"Ligands/mol3_original.mrc",MapFileType.MRC); */

		//List<InvariantNorm> zc = new ArrayList<>();

		InvariantNorm normalization1 = new InvariantNorm(vol1, maxorder);
		InvariantNorm normalization2 = new InvariantNorm(vol2, maxorder);
		InvariantNorm normalization3 = new InvariantNorm(vol3, maxorder);
		InvariantNorm normalization4 = new InvariantNorm(vol4, maxorder);
		InvariantNorm normalization5 = new InvariantNorm(vol5, maxorder);


        double d = 100;
        d = normalization1.compareInvariants(normalization2, 0); // 0 is basic fingerprints
        System.out.println("Distance 1 to 2: " + d);
                        
        d = normalization1.compareInvariants(normalization3, 0); // 0 is basic fingerprints
        System.out.println("Distance 1 to 3: " + d);
        
        d = normalization1.compareInvariants(normalization4, 0); // 0 is basic fingerprints
        System.out.println("Distance 1 to 4: " + d);

        d = normalization1.compareInvariants(normalization5, 0); // 0 is basic fingerprints
        System.out.println("Distance 1 to 5: " + d);

        d = normalization4.compareInvariants(normalization5, 0); // 0 is basic fingerprints
        System.out.println("Distance 4 to 5: " + d);
                        
        AlignmentResult alignmentResult =  align(normalization1, normalization2, maxorder);
        System.out.println("ALIGNMENT SCORE (1 - 2): "+alignmentResult.getScore());
		Matrix4d rot12 = alignmentResult.getTransforms().get(0);
		Matrix4d rot2 = alignmentResult.getTransforms().get(1);
        
        alignmentResult =  align(normalization1, normalization3, maxorder);
        System.out.println("ALIGNMENT SCORE (1 - 3): "+alignmentResult.getScore());
		Matrix4d rot13 = alignmentResult.getTransforms().get(0);
		Matrix4d rot3 = alignmentResult.getTransforms().get(1);
        
        alignmentResult =  align(normalization1, normalization4, maxorder);
        System.out.println("ALIGNMENT SCORE (1 - 4): "+alignmentResult.getScore());
		Matrix4d rot14 = alignmentResult.getTransforms().get(0);
		Matrix4d rot4 = alignmentResult.getTransforms().get(1);
        
        alignmentResult =  align(normalization1, normalization5, maxorder);
        System.out.println("ALIGNMENT SCORE (1 - 5): "+alignmentResult.getScore());
		Matrix4d rot15 = alignmentResult.getTransforms().get(0);
		Matrix4d rot5 = alignmentResult.getTransforms().get(1);
        
        /* rot1.invert(); // this is the way of mol1 to the center we need to undo this way by inverting the matrix
        rot2.mul(rot1);
        rot3.mul(rot1);
        rot4.mul(rot1);
        rot5.mul(rot1); */
        /* Matrix4d rot21 = (Matrix4d) rot1.clone();
        rot21.mul(rot2); // concatenate transformation of molecule 2 to the center followed by the way of mol1 to the original position
        Matrix4d rot31 = (Matrix4d) rot1.clone();
        rot31.mul(rot3);
        Matrix4d rot41 = (Matrix4d) rot1.clone();
        rot41.mul(rot4);
        Matrix4d rot51 = (Matrix4d) rot1.clone();
        rot51.mul(rot5); */

		//Calc.transform(mol1,rot1);
		Calc.transform(mol1,rot13); // superpose mol2 onto mol1
		Calc.transform(mol2,rot2); // superpose mol2 onto mol1
		Calc.transform(mol3,rot3); // superpose mol3 onto mol1
		Calc.transform(mol4,rot4); // superpose mol4 onto mol1
		Calc.transform(mol5,rot5); // superpose mol5 onto mol1
		/* Calc.transform(mol2,rot21); // superpose mol2 onto mol1
		Calc.transform(mol3,rot31); // superpose mol3 onto mol1
		Calc.transform(mol4,rot41); // superpose mol4 onto mol1
		Calc.transform(mol5,rot51); // superpose mol5 onto mol1 */

		// Create a new Structure object for the ligand
		writeStructurePDB(mol1, "Ligands/mol1.original.pdb");
		writeStructurePDB(mol1, "Ligands/mol1.fitted3.pdb");
		writeStructurePDB(mol2, "Ligands/mol2.fitted.pdb");
		writeStructurePDB(mol3, "Ligands/mol3.fitted.pdb");
		writeStructurePDB(mol4, "Ligands/mol4.fitted.pdb");
		writeStructurePDB(mol5, "Ligands/mol5.fitted.pdb");
	}


	@Test
	public void testPhm() throws Exception {
        Volume vol1 = loadMoleculeVolume("phmscreen_ref_hydroele.dx");
        Volume vol2 = loadMoleculeVolume("phmscreen_good_hydroele.dx");
        Volume vol3 = loadMoleculeVolume("phmscreen_bad_hydroele.dx");
        Structure mol1 = loadPDB("phmscreen_ref.pdb");
        Structure mol2 = loadPDB("phmscreen_good.pdb");
        Structure mol3 = loadPDB("phmscreen_bad.pdb");
		/* vol1.setRadiusVarMult(2);
		vol2.setRadiusVarMult(2);
		vol3.setRadiusVarMult(2); */
		
		/* VolumeIO.write(vol1,"Ligands/mol1_original.mrc",MapFileType.MRC);
		VolumeIO.write(vol2,"Ligands/mol2_original.mrc",MapFileType.MRC);
		VolumeIO.write(vol3,"Ligands/mol3_original.mrc",MapFileType.MRC); */

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
        
        System.out.println("ALIGNMENT SCORE (global): "+alignmentResult.getScore());
        
		Matrix4d rot1 = alignmentResult.getTransforms().get(0);
		Matrix4d rot2 = alignmentResult.getTransforms().get(1);
		Matrix4d rot3 = alignmentResult.getTransforms().get(2);
        rot1.invert(); // this is the way of mol1 to the center we need to undo this way by inverting the matrix
        Matrix4d rot21 = (Matrix4d) rot1.clone();
        rot21.mul(rot2); // concatenate transformation of molecule 2 to the center followed by the way of mol1 to the original position
        Matrix4d rot31 = (Matrix4d) rot1.clone();
        rot31.mul(rot3);

		//Calc.transform(mol1,rot1);
		Calc.transform(mol2,rot21); // superpose mol2 onto mol1
		Calc.transform(mol3,rot31); // superpose mol3 onto mol1

		// Create a new Structure object for the ligand
		writeStructurePDB(mol1, "Ligands/phmscreen_ref.original.pdb");
		writeStructurePDB(mol2, "Ligands/phmscreen_good.fitted.pdb");
		writeStructurePDB(mol3, "Ligands/phmscreen_bad.fitted.pdb");
	}


	@Test
	public void testReconstruction() throws Exception {
        int maxOrder = 15;
        String mol = "mol5";
        String probe = "HDON";
        String basename = mol+"_"+probe;
		Volume volume1 = loadMoleculeVolume(basename+".dx");
		ZernikeMoments zm = new ZernikeMoments(volume1, maxOrder);
        //describeVolume(volume1);
        //		Volume volume_rec = ZernikeMoments.reconstructVolume(zm, 32, 15, false, true);
        //		VolumeIO.write(volume_rec, "D:/PT/Ligands/volume1_rec_orig.map", MapFileType.MRC);
        
		for (int maxN=5;maxN<=maxOrder;maxN+=5) {
            int[] nTerms = CalcNumZernikeTerms.calcNumTerms(maxN);
            System.out.println("ORDER "+maxN+" Number of coefficients: "+ nTerms[0] +" total "+ nTerms[1] + " positive "+ nTerms[2] + " invariant");
			Volume volume_rec2 = ZernikeMoments.reconstructVolume(zm, volume1.getDimensions(), volume1.getCenterVolume(), maxN, volume1.getGridWidth(),false, false);
            volume_rec2.setCorner(volume1.getCorner());
            volume_rec2.setCenterReal(volume1.getCenterReal());
            //describeVolume(volume_rec2);
			VolumeIO.write(volume_rec2, "Ligands/orders/"+basename+"_"+maxN+".mrc", MapFileType.MRC);
		}
	}

    @Test
    public void testLigands() throws Exception {
        int max_order = 15;
        String[] fileNames = { "mol1_HDON.dx", "mol2_HDON.dx", "mol3_HDON.dx", "mol4_HDON.dx", "phmscreen_ref_hydroele.dx" };
        int n_files = fileNames.length;
        Volume[] volumeArray = new Volume[n_files];
        List<InvariantNorm> invariantsArray = new ArrayList<>();
        //InvariantNorm[] invariantsArray = new InvariantNorm[n_files];

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

}
