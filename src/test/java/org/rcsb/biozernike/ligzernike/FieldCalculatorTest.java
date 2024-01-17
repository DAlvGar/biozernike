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
        String path = SDFReaderTest.class.getResource("/test.sdf").getPath();
        int fieldID = 6; // Set the desired field ID

        try (IteratingSDFReader reader = SDFReader.read(path)) {
            IAtomContainer molecule;
            while (reader.hasNext()) {
                molecule = reader.next();
                // Process coordinates of atoms in the molecule
                Volume volume = FieldCalculator.projectField(molecule, fieldID);
                OpenDXIO.write("project_acceptors.dx", volume);
            }
            // Close the SDF reader
            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
