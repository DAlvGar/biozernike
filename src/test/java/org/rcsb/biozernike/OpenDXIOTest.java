package org.rcsb.biozernike;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.net.URL;
import org.junit.Test;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.OpenDXIO;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.volume.VolumeIO;


public class OpenDXIOTest {
    
    ClassLoader classLoader = getClass().getClassLoader();

	@Test
	public void testDXReading() throws Exception {

		// some hardcoded scaling coefficients for the EM volume (as we do not control the density values)
		URL url = OpenDXIOTest.class.getResource("/phmscreen_ref_hydroele.dx");
		Volume volumeDX = OpenDXIO.read(url.getPath());
        assertEquals(0.5, volumeDX.getGridWidth(), 0);
    }

	@Test
	public void testDXWriting() throws Exception {
		// some hardcoded scaling coefficients for the EM volume (as we do not control the density values)
		URL url = OpenDXIOTest.class.getResource("/mol1_HDON.dx");
		// URL url = OpenDXIOTest.class.getResource("/phmscreen_ref_hydroele.dx");
		Volume volumeDX = OpenDXIO.read(url.getPath());
		OpenDXIO.write("reconstructed.dx", volumeDX);
		volumeDX.flipValues();
		//OpenDXIO.write("reconstructed_flip.dx", volumeDX);
		volumeDX.applyContourAndNormalize(1, 1);
		OpenDXIO.write("reconstructed_flip_contour.dx", volumeDX);
    }

	@Test
	public void testConversion() throws Exception {
		// some hardcoded scaling coefficients for the EM volume (as we do not control the density values)
		// URL url = OpenDXIOTest.class.getResource("/phmscreen_ref_hydroele.dx");
		URL url = OpenDXIOTest.class.getResource("/mol1_HDON.dx");
		Volume volumeDX = OpenDXIO.read(url.getPath());
		OpenDXIO.write("reconstructed.dx", volumeDX);
		VolumeIO.write(volumeDX, new File("convert.ccp4"), MapFileType.CCP4);
		VolumeIO.write(volumeDX, new File("convert.mrc"), MapFileType.MRC);
		File infile = new File("convert.mrc");
		//Volume ccp4 = VolumeIO.read(infile, MapFileType.MRC);
		//OpenDXIO.write("ccp4todx.dx", ccp4);
    }
}
