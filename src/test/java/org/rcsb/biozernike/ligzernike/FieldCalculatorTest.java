package org.rcsb.biozernike.ligzernike;

import java.io.IOException;

import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.rcsb.biozernike.molecules.SDFReader;
import org.rcsb.biozernike.molecules.SDFReaderTest;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;

public class FieldCalculatorTest {

    @Test
    public void testFieldCalculation() {
        String path = SDFReaderTest.class.getResource("/test_translated.sdf").getPath();
        int[] fields = { 2, 3, 5, 6 }; // Set the desired field ID
        String[] names = { "hydroele", "hydrocav", "hbond_Donors", "hbond_Acceptors" };
        double[] spacing = {0.5, 0.5, 0.5, 0.5};
        FieldCalculator calculator = new FieldCalculator();
        
        try (IteratingSDFReader reader = SDFReader.read(path)) {
            IAtomContainer molecule;
            while (reader.hasNext()) {
                molecule = reader.next();
                // Process coordinates of atoms in the molecule
                for (int i = 0; i < fields.length; i++) {
                    calculator.setSpacing(spacing[i]);
                    Volume volume = calculator.projectField(molecule, fields[i]);
                    OpenDXIO.write("Ligands/projection_" + names[i] + ".dx", volume);
                }
            }
            // Close the SDF reader
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
