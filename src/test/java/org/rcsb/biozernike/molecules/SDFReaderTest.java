package org.rcsb.biozernike.molecules;

import java.util.List;

import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

public class SDFReaderTest {

    @Test
    public void testReadSDF() throws Exception {
        String path = SDFReaderTest.class.getResource("/test.sdf").getPath();
        IteratingSDFReader reader = SDFReader.read(path);
        // Iterate through the molecules in the SDF file
        while (reader.hasNext()) {
            IAtomContainer molecule = reader.next();
            // Process coordinates of atoms in the molecule
            List<List<Double>> f = SDFReader.getAtomFields(molecule);
            System.out.println(f);
        }

        // Close the SDF reader
        reader.close();
    }
}
